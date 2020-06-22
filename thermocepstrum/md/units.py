def scale_kappa_REALtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units.
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    return (4184. / NA / temp)**2 / kB / volume * timestep * 100.

def scale_kappa_viscosity_REALtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units (1 Poise=0.1 Pa*s).
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    atmtoPa=1.01325
   # print('kappascale= ',volume * atmtoPa**2  / kB / temp  * timestep * 10**(-12))
   # print('vol= ',volume)
   # print('atmtoPa= ',atmtoPa)
   # print('kb= ',kB)
   # print('temp= ',temp)
   # print('timestep= ',timestep)
    return volume * atmtoPa**2  / kB / temp  * timestep * 1.0e-11

def scale_kappa_viscosity_CPtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units (1 Poise=0.1 Pa*s).
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ a_bohr^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    atmtoPa=1.01325
    a_bohr = 0.529177
    # print('kappascale= ',volume * atmtoPa**2  / kB / temp  * timestep * 10**(-12))
    # print('vol= ',volume)
    # print('atmtoPa= ',atmtoPa)
    # print('kb= ',kB)
    # print('temp= ',temp)
    # print('timestep= ',timestep)
    return volume *a_bohr**3 / kB / temp  * timestep * 1.0e-03

def scale_kappa_viscosity_METALtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units (1 Poise=0.1 Pa*s).
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    atmtoPa=1.01325
    a_bohr = 0.529177
    # print('kappascale= ',volume * atmtoPa**2  / kB / temp  * timestep * 10**(-12))
    # print('vol= ',volume)
    # print('atmtoPa= ',atmtoPa)
    # print('kb= ',kB)
    # print('temp= ',temp)
    # print('timestep= ',timestep)
    return volume / kB / temp * timestep * 1.0e-11

def scale_kappa_METALtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from METAL LAMMPS units to SI units.
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    return (charge / temp)**2 / kB / volume * timestep * 10000.


def scale_kappa_QEPWtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from Quantum Espresso PW/HARTREE-HEATCURRENT units to SI units.
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    J_PWtoMETAL = 1.0 / 6.719329152e-6   # [Ry*Bohr/tau_PW] --> [ev*A/ps]  (1tau_PW = 4.8378e-5 ps)
    return (charge / temp)**2 / kB / volume * timestep * 10000. * J_PWtoMETAL**2


def scale_kappa_GPUMDtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from GPUMD units to SI units.
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.02214
    massunit = 1.660538921
    charge = 1.6021765
    return (charge)**3 / (temp)**2 / massunit / kB / volume * timestep * 1.0e8


def scale_kappa_DLPOLYtoSI(temp, volume, timestep):
    """
    Conversion factor for the thermal conductivity from DL_POLY units to SI units.
    INPUT:
    temp      =  temperature [ K ]
    volume    =  cell volume [ A^3 ]
    timestep  =  integration time step [ fs ]
    """

    kB = 1.3806504
    NA = 6.022140857
    return (1.0 / NA / temp)**2 / kB / volume * timestep * 1e10

def scale_kappa_METAL_ELtoSI(temp, volume, timestep):
    """
    Conversion factor for the electrical conductivity from METAL LAMMPS units to [S/cm] units.
    INPUT:
    : param temp :  temperature [ K ]
    : param volume :  cell volume [ A^3 ]
    :param timestep  :  integration time step [ fs ]
    """

    kB = 1.3806504
    charge = 1.6021765
    return charge**2 / temp / kB / volume * timestep * 100
