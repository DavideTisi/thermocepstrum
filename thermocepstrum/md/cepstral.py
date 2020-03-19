import numpy as np
from scipy.special import polygamma
from scipy.sparse import diags
from scipy.fftpack import dct,rfft,irfft, fft, ifft
from .tools import logtau_to_tau
from .aic import *
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()

EULER_GAMMA = 0.57721566490153286060651209008240243104215933593992   # Euler-Mascheroni constant

################################################################################


def multicomp_cepstral_parameters(NF, N_COMPONENTS):
    """
    Returns the theoretical variance of the cepstral coefficients and the mean of the log(PSD) distribution,
    generated from a periodogram that is the average of N_COMPONENTS.
    """

    N = 2 * (NF - 1)

    # variance of cepstral coefficients
    trigamma = polygamma(1, N_COMPONENTS)
    ck_THEORY_var = 1. / N * np.concatenate(([2 * trigamma], [trigamma] * (NF - 2), [2 * trigamma]))

    # bias of log(PSD)
    psd_THEORY_mean = (polygamma(0, N_COMPONENTS) - np.log(N_COMPONENTS)) * np.ones(NF)
    psd_THEORY_mean[0] = polygamma(0, 0.5 * N_COMPONENTS) - np.log(0.5 * N_COMPONENTS)
    psd_THEORY_mean[-1] = psd_THEORY_mean[0]

    return ck_THEORY_var, psd_THEORY_mean

def mel_multicomp_cepstral_parameters( N_COMPONENTS,bins):
    """
    Returns the theoretical variance of the cepstral coefficients and the mean of the log(PSD) distribution,
    generated from a periodogram that is the average of N_COMPONENTS.

    :param N_COMPONENTS : number of degrees of freedom
    :param bins : filterbank bins
    """
    NF = bins.shape[0]
    N = 2 * (NF - 1)
    Nbins = np.zeros(NF-2)
    #Nbins[0] = Nbins[-1] = 1
    #Nbins[1:-1] = bins[2:] - bins[:-2]
    Nbins = bins[2:] - bins[:-2] + 1
    T = 1/Nbins
    # variance of cepstral coefficients
    trigamma = polygamma(1, N_COMPONENTS)
    ck_THEORY_var = 1. / N * np.concatenate(([2 * trigamma], [trigamma] * (NF - 2), [2 * trigamma]))

    ### compute the covariance matrix for Xi_j Xi_i
    var_diag = T * trigamma
    #var_sdiag = np.zeros(NF-1)
    var_sdiag = np.zeros(NF-3)
    #var_sdiag [1:-2]=(bins[3:-1] - bins[2:-2])*T[1:-3]*T[2:-2]*trigamma
    var_sdiag = (bins[2:-1] - bins[1:-2] + 1) * T[:-1] * T[1:] * trigamma
    
    
    # bias of log(PSD)
    #psd_THEORY_mean = (polygamma(0, N_COMPONENTS) - np.log(N_COMPONENTS)) * np.ones(NF)
    #psd_THEORY_mean[0] = polygamma(0, 0.5 * N_COMPONENTS) - np.log(0.5 * N_COMPONENTS)
    #psd_THEORY_mean[-1] = psd_THEORY_mean[0]
    psd_THEORY_mean = (polygamma(0, N_COMPONENTS) - np.log(N_COMPONENTS)) * np.ones(NF-2)

    return ck_THEORY_var, psd_THEORY_mean, [var_diag, var_sdiag]

def dct_coefficients(y):
    """Compute the normalized Discrete Cosine Transform coefficients of y.
        yk = 0.5 * DCT(y) / (N-1)"""
    yk = dct(y, type=1) / (y.size - 1) * 0.5   # normalization
    return yk


def dct_filter_psd(y, K=None):
    # K=P*-1 is the maximum coefficient summed (c_k = 0 for k > K)
    if (K >= y.size):
        log.write_log('! Warning:  dct_filter_psd K value ({:}) out of range.'.format(K))
        return np.full(y.size, np.NaN)
    yk = dct(y, type=1)
    if K is not None:
        yk[K + 1:] = 0.
    ynew = dct(yk, type=1) / (y.size - 1) * 0.5
    return ynew


def dct_filter_tau(y):
    # K=P*-1 is the maximum coefficient summed (c_k = 0 for k > K)
    yk = dct(y, type=1) / (y.size - 1)
    ftau = np.zeros(y.size)
    ftau[0] = yk[0]
    ftau[1] = 0.5 * yk[0] + yk[1]
    for i in range(2, yk.size - 1):
        ftau[i] = ftau[i - 1] + yk[i]
    ftau[-1] = ftau[-2] + 0.5 * yk[-1]
    return ftau


################################################################################


class CosFilter(object):
    """
    CEPSTRAL ANALYSIS based filtering.

    ** INPUT VARIABLES:
    samplelogpsd    = the original sample log-PSD, \hat{L}_k
    ck_theory_var   = the theoretical variance of cepstral coefficients, \sigma*^2(P*,N)
    psd_theory_mean = the theoretical bias of log-PSD, \lambda_l
    aic_type        = type of AIC to use ('aic' (default), 'aicc')
    Kmin_corrfactor = cutoff correction factor (default: 1.0)
    K_PSD           = cutoff used to compute logpsd (default: K_PSD = aic_Kmin)

    ** INTERNAL VARIABLES:
    samplelogpsd  = the original sample log-PSD - logpsd_THEORY_mean

    logpsdK  = the cepstrum of the data, \hat{C}_n (i.e. the DCT of samplelogpsd)
    aic_min  = minimum value of the AIC
    aic_Kmin = cutoff K that minimizes the AIC, K = P*-1

    logtau          = filtered log(tau) as a function of K, L_0(P*-1)
    logtau_Kmin     = filtered log(tau) at the aic_Kmin, L*_0
    logtau_var_Kmin = theoretical L*_0 variance
    logtau_std_Kmin = theoretical L*_0 standard deviation
    logpsd          = filtered log-PSD at the specified cutoff K_PSD

    tau          = filtered tau as a function of K, S_0(P*-1)
    tau_Kmin     = filtered tau at the aic_Kmin, S*_0
    tau_var_Kmin = theoretical S*_0 variance
    tau_std_Kmin = theoretical S*_0 standard deviation
    psd          = filtered PSD at the specified cutoff K_PSD

    p_aic... = Bayesian AIC weighting stuff
    """

    def __init__(self, samplelogpsd, ck_theory_var=None, psd_theory_mean=None, aic_type='aic', Kmin_corrfactor=1.0):

        NF = samplelogpsd.size
        N = 2 * (NF - 1)

        if psd_theory_mean is None:
            # by default the THEORETICAL means are the one component ones:
            # ck THEORY mean:
            #    - EULER_GAMMA - log(2)   for k = {0, N/2}
            #    - EULER_GAMMA            otherwise
            self.logpsd_THEORY_mean = -EULER_GAMMA * np.ones(NF)
            self.logpsd_THEORY_mean[0] = -EULER_GAMMA - np.log(2)
            self.logpsd_THEORY_mean[-1] = -EULER_GAMMA - np.log(2)
        else:
            self.logpsd_THEORY_mean = psd_theory_mean

        # subtract the mean of the distribution
        self.samplelogpsd = samplelogpsd - self.logpsd_THEORY_mean

        # compute cepstral coefficients
        self.logpsdK = dct_coefficients(self.samplelogpsd)

        # estimate AIC
        if (aic_type == 'aic'):
            self.aic = dct_AIC(self.logpsdK, ck_theory_var)
        elif (aic_type == 'aicc'):
            self.aic = dct_AICc(self.logpsdK, ck_theory_var)
        else:
            raise ValueError('AIC type not valid.')
        self.aic_type = aic_type
        self.aic_min = np.min(self.aic)
        self.Kmin_corrfactor = Kmin_corrfactor
        self.aic_Kmin = int(round(np.argmin(self.aic) * Kmin_corrfactor))
        if (self.aic_Kmin >= NF):
            log.write_log('! Warning:  aic_Kmin ({:}) is out of range.'.format(self.aic_Kmin))

        # set theoretical errors
        if ck_theory_var is None:
            # by default the THEORETICAL variances are the one component ones:
            # ck THEORY variances:
            #    (pi^2)/3/N   for k = {0, N/2}
            #    (pi^2)/6/N   otherwise
            #self.logpsdK_THEORY_var = 1. / N * np.concatenate(
            #    ([np.pi**2 / 3], [np.pi**2 / 6.] * (NF - 2), [np.pi**2 / 3]))
            self.logpsdK_THEORY_var = 1. / N * np.concatenate(
                ([np.pi**2 / 3], [np.pi**2 / 6.] * (NF - 3), [np.pi**2 / 3]))
            self.logpsdK_THEORY_std = np.sqrt(self.logpsdK_THEORY_var)
            # logtau THEORY variances:  (we assume to be summing ck up to K, included)
            #    (pi^2)/3/N*(2*K+1)   for K = {0, N/2-1}
            #    (pi^2)/3             for K = N/2
            self.logtau_THEORY_var = 1. / N * np.concatenate(
                (np.pi**2 / 3. * (2 * np.arange(NF - 1) + 1), [np.pi**2 / 3. * N]))
            self.logtau_THEORY_std = np.sqrt(self.logtau_THEORY_var)
        else:
            self.logpsdK_THEORY_var = ck_theory_var
            self.logpsdK_THEORY_std = np.sqrt(self.logpsdK_THEORY_var)
            self.logtau_THEORY_var = np.zeros(NF)
            self.logtau_THEORY_var[0] = self.logpsdK_THEORY_var[0]
            for K in range(1, NF - 1):
                self.logtau_THEORY_var[K] = self.logtau_THEORY_var[K - 1] + 4. * self.logpsdK_THEORY_var[K]
            self.logtau_THEORY_var[-1] = self.logtau_THEORY_var[-2] + self.logpsdK_THEORY_var[-1]
            self.logtau_THEORY_std = np.sqrt(self.logtau_THEORY_var)
        return

    def __repr__(self):
        msg = 'CosFilter:\n' + \
              '  AIC type  = {:}\n'.format(self.aic_type) + \
              '  AIC min   = {:f}\n'.format(self.aic_min) + \
              '  AIC_Kmin  = {:d}  (P* = {:d})\n'.format(self.aic_Kmin, self.aic_Kmin + 1) + \
              '  L_0*   = {:15f} +/- {:10f}\n'.format(self.logtau_Kmin, self.logtau_std_Kmin) + \
              '  S_0*   = {:15f} +/- {:10f}\n'.format(self.tau_Kmin, self.tau_std_Kmin) + \
              '  K_PSD  = {:d}\n'.format(self.K_PSD)
        return msg

    def scan_filter_tau(self, K_PSD=None, correct_mean=True):
        """Computes the tau as a function of the cutoff K for the CosFilter.
        Also computes psd and logpsd for the given K_PSD cutoff (or otherwise the aic_Kmin is used)."""
        if K_PSD is None:
            self.K_PSD = self.aic_Kmin
        else:
            self.K_PSD = K_PSD
            self.aic_Kmin = self.K_PSD

        if (self.aic_Kmin >= self.samplelogpsd.size):
            log.write_log('! Warning:  aic_Kmin ({:}) is out of range.'.format(self.aic_Kmin))

        # COS-filter analysis with frequency cutoff K
        self.logtau = dct_filter_tau(self.samplelogpsd)
        self.logpsd = dct_filter_psd(self.samplelogpsd, self.K_PSD)   # usually is log(psd)@aic_Kmin
        self.psd = np.exp(self.logpsd)
        self.tau = np.exp(self.logtau)
        self.tau_THEORY_std = self.tau * self.logtau_THEORY_std

        if (self.aic_Kmin < self.samplelogpsd.size):
            self.logtau_Kmin = self.logtau[self.aic_Kmin]
            self.logtau_var_Kmin = self.logtau_THEORY_var[self.aic_Kmin]
            self.logtau_std_Kmin = self.logtau_THEORY_std[self.aic_Kmin]
            self.tau_Kmin = self.tau[self.aic_Kmin]
            self.tau_std_Kmin = self.tau_THEORY_std[self.aic_Kmin]
            self.tau_var_Kmin = self.tau_std_Kmin**2
        else:
            self.logtau_Kmin = np.NaN
            self.logtau_var_Kmin = np.NaN
            self.logtau_std_Kmin = np.NaN
            self.tau_Kmin = np.NaN
            self.tau_var_Kmin = np.NaN
            self.tau_std_Kmin = np.NaN

        if correct_mean:
            self.logpsd = self.logpsd + self.logpsd_THEORY_mean
            self.logtau = self.logtau + self.logpsd_THEORY_mean[0]
            self.logtau_Kmin = self.logtau_Kmin + self.logpsd_THEORY_mean[0]
        return

    def scan_filter_psd(self, K_LIST, correct_mean=True):
        """Computes the psd as a function of the cutoff K for the CosFilter.
        Repeats the procedure for all the cutoffs in the K_LIST."""
        self.K_LIST = K_LIST
        self.logpsd_K_LIST = np.zeros((self.samplelogpsd.size, len(self.K_LIST)))
        self.psd_K_LIST = np.zeros((self.samplelogpsd.size, len(self.K_LIST)))
        self.logtau_K_LIST = np.zeros(len(self.K_LIST))   # DEFINED AS log(PSD[0]), no factor 0.5 or 0.25
        self.tau_K_LIST = np.zeros(len(self.K_LIST))

        for k, K in enumerate(self.K_LIST):
            # COS-filter analysis with frequency cutoff K
            self.logpsd_K_LIST[:, k] = dct_filter_psd(self.samplelogpsd, K)
            self.logtau_K_LIST[k] = self.logpsd_K_LIST[0, k]

            self.psd_K_LIST[:, k] = np.exp(self.logpsd_K_LIST[:, k])
            self.tau_K_LIST[k] = np.exp(self.logtau_K_LIST[k])

            if correct_mean:
                self.logpsd_K_LIST[:, k] = self.logpsd_K_LIST[:, k] + self.logpsd_THEORY_mean
                self.logtau_K_LIST[k] = self.logtau_K_LIST[k] + self.logpsd_THEORY_mean[0]
        return

    #############################
    ####  Bayesian method
    #############################
    def compute_p_aic(self, method='ba'):
        """Define a weight distribution from the AIC, according to a method."""
        NF = self.samplelogpsd.size
        self.p_aic = produce_p(self.aic, method)
        self.p_aic_Kave, self.p_aic_Kstd = grid_statistics(np.arange(NF), self.p_aic)
        return

    def compute_logtau_density(self, method='ba', only_stats=False, density_grid=None, grid_size=1000,
                               correct_mean=True):
        if self.p_aic is None:
            raise ValueError('No P_AIC defined.')

        # compute statistics
        self.p_logtau_density_xave, self.p_logtau_density_xstd = \
                        grid_statistics(self.logtau, self.p_aic, self.logtau_THEORY_var + self.logtau**2)
        self.p_logtau_density_xstd2 = np.dot(
            self.p_aic, np.sqrt(self.logtau_THEORY_var + (self.logtau - self.p_logtau_density_xave)**2))
        ##self.p_logtau_density_xave, self.p_logtau_density_xstd = \
        ##                ta.grid_statistics(self.p_logtau_grid, self.p_logtau_density)

        # compute distribution
        if not only_stats:
            if density_grid is None:
                self.p_logtau_density, self.p_logtau_grid = produce_p_density(self.p_aic, \
                                        self.logtau_THEORY_std, self.logtau, grid_size=grid_size)
            else:
                self.p_logtau_grid = density_grid
                self.p_logtau_density = produce_p_density(self.p_aic, self.logtau_THEORY_std, \
                                            self.logtau, grid=self.p_logtau_grid)

        # tau distribution
        self.p_tau_density_xave, self.p_tau_density_xstd = logtau_to_tau(self.p_logtau_density_xave,
                                                                         self.logpsd_THEORY_mean[0],
                                                                         self.p_logtau_density_xstd)
        return

    def mel_compute_variance(self,mel_var_list,debug=False):
        '''

        :param mel_var_list: list with diagonal and convariance matrix of Xi (vedi Notre Mel #TODO spiegare)
        :param debug: debug flag if True the code return also covariance matrix
        :return: variance on the mel-filtered cepstrum, civariance matrix (only if debug=True)
        '''
        if debug: cov = diags([mel_var_list[0], mel_var_list[1], mel_var_list[1]],[0,1,-1]).toarray() #cov= covariance Xi
        

        
        ### Vecchia formula approssimata # Per velocizzare i conti: calcolo solo l'errore in zero, ovvero quello che diventerà l'errore su \kappa
        ###cov_cc = ifft(fft(cov, axis = 1), axis = 0)  #/cov.shape[0]
        ###cov_cc[self.aic_Kmin + 1:,:] = 0.
        ###cov_cc[:,self.aic_Kmin + 1:] = 0.
        ###tmp1 = np.sum(cov_cc).real/cov.shape[0]*np.ones(cov.shape)
        

        if(debug):
              (n1, n2) = cov.shape
              r1 = np.arange(n1)
              r2 = np.arange(n2)

              eps = np.exp(2*np.pi*1j*np.outer(r1, r2)/n1)
              ems = 1/eps # np.exp(-2*np.pi*1j*np.outer(r1, r2)/n1)
              ep  = np.copy(eps) #np.exp(2*np.pi*1j*np.outer(r1, r2)/n1)
              em  = 1/ep #np.exp(-2*np.pi*1j*np.outer(r1, r2)/n1)
              #eps[self.aic_Kmin + 1:,self.aic_Kmin + 1:] = 0.
              #ems[self.aic_Kmin + 1:,self.aic_Kmin + 1:] = 0.
              eps[self.aic_Kmin + 1:,:] = 0.
              eps[:,self.aic_Kmin + 1:] = 0.
              ems[self.aic_Kmin + 1:,:] = 0.
              ems[:,self.aic_Kmin + 1:] = 0.
              tmp = np.einsum('am,bn,jn,mi,ji->ab', eps, ems, ep, em, cov, optimize='greedy').real

        ### Correct formula (fft makes the calculation efficient; possible improvements in the calculation of s5)
        nfilt = len(mel_var_list[0])
        nj = mel_var_list[0]
        alphaj = np.zeros(nfilt)
        alphaj[:-1] = mel_var_list[1]
        twopinf = 2*np.pi/nfilt

        fftnj = fft(nj)
        fftalpha = fft(alphaj)
        pstar = self.aic_Kmin + 1
        range0 = range(1, pstar)
        range2 = range(-(pstar-1), 0)
        ###
        s1 = np.sum(nj + 2*alphaj)
        ###
        tmp = fftnj.real + fftalpha.real
        s2 = np.sum(tmp[1:pstar])
        ###
        tmp = fftnj.real
        s3 = np.sum([tmp.take(n1, mode='wrap') for n in range2 for n1 in range(n+1, n+pstar)])
        ###
        s4 = np.sum([fftalpha[n]*np.exp(-1j*twopinf*n) for n in range0]).real
        ###
        pinf = np.pi/nfilt
        s5 = np.sum([np.cos(pinf*(n+n1)) * np.real(np.exp(-1j*pinf*(n-n1)) * \
                                                   fftalpha.take(n-n1, mode='wrap')) \
                      for n in range0 \
                    for n1 in range0]).real
    
        v =  s1 + 4*(s2 + s3 + s4) + 8*s5

        if debug : return np.sqrt(v)/nfilt, cov, np.sqrt(np.diag(tmp)/n1/n2)
        return np.sqrt(v)/nfilt

#    def optimize_cos_filter(self, thr=0.05, K_LIST=None, logtauref=None):
#        if K_LIST is not None:
#            self.K_LIST = K_LIST
#        self.scan_cos_filter_K()
#        ## find minimum cutoff K that satisfies  |log(tau) - tauref| < thr
#        if logtauref is not None:
#            self.logtauref = logtauref
#        else:
#            self.logtauref = self.logtau[-1]  # if tauref is not given, use logtau with max cutoff
#        self.optimalK_idx = len(self.K_LIST) - np.argmin(np.abs(self.logtau - self.logtauref)[::-1] <= thr)
#        if (self.optimalK_idx < len(self.K_LIST)):
#            self.optimalK = self.K_LIST[self.optimalK_idx]
#        else:
#            self.optimalK_idx = np.NaN
#            self.optimalK = np.NaN
#            log.write_log(Warning: optimal cutoff K NOT FOUND.')
#        return
