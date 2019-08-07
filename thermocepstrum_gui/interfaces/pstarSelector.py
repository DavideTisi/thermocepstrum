from thermocepstrum_gui.utils.custom_widgets import *
from thermocepstrum_gui.core import control_unit as cu
from thermocepstrum_gui.core.control_unit import log


class PStarSelector(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main

        self.next_frame = None
        self.prev_frame = None

        self.parent = parent
        self.main_frame = self

        self.main_frame.grid(column=0, row=0, sticky='nsew', padx=20)

        sections = Frame(self.main_frame)
        sections.grid(row=0, column=0, sticky='nsew', padx=20, columnspan=2)

        self.graph = GraphWidget(sections, sections, size=(7, 4), toolbar=True)

        self.container_frame = Frame(self.main_frame)
        self.container_frame.grid(row=1, column=0, sticky='nsew')

        variable_frame = Frame(self.container_frame, bd=1, relief=SOLID)
        variable_frame.pack(side=TOP, anchor='w', padx=20, fill='x', expand=1, pady=20)

        Label(variable_frame, text='f*: ', font='Arial 12 bold').grid(row=0, column=0, sticky='we')
        self.fstar_label = Label(variable_frame, text='')
        self.fstar_label.grid(row=0, column=1, sticky='we')
        Label(variable_frame, text='P*: ', font='Arial 12 bold').grid(row=0, column=2, sticky='we')
        self.pstar_label = Label(variable_frame, text='')
        self.pstar_label.grid(row=0, column=3, sticky='we')
        Label(variable_frame, text='\u03f0:', font='Arial 12 bold').grid(row=0, column=4, sticky='we')
        self.kmin_label = Label(variable_frame, text='')
        self.kmin_label.grid(row=0, column=5, sticky='we')

        value_frame = Frame(self.container_frame)
        value_frame.pack(side=TOP, anchor='w', padx=20, fill=BOTH, expand=1)

        Label(value_frame, text='P* selector', font='Arial 12 bold').grid(row=0, column=0, sticky='w')
        ttk.Separator(value_frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we', pady=10, columnspan=5)

        Label(value_frame, text='P*: ', font='Arial 12').grid(row=2, column=0, sticky='w')
        self.value_entry = Spinbox(value_frame, bd=1, relief=SOLID, increment=1)
        self.value_entry.grid(row=2, column=1, sticky='w')

        self.increment = IntVar()
        Label(value_frame, text='Increment by: ', font='Arial 12').grid(row=3, column=0, sticky='w', pady=10)

        rdbt_frame = Frame(value_frame)
        rdbt_frame.grid(row=3, column=1, sticky='w')

        Radiobutton(rdbt_frame, text='1', font='Arial 11 bold', variable=self.increment, value=1,
                    command=self._change_increment).pack(side=LEFT)
        Radiobutton(rdbt_frame, text='10', font='Arial 11 bold', variable=self.increment, value=10,
                    command=self._change_increment).pack(side=LEFT)
        Radiobutton(rdbt_frame, text='100', font='Arial 11 bold', variable=self.increment, value=100,
                    command=self._change_increment).pack(side=LEFT)

        Button(value_frame, text='Recalculate', font='Arial 12 bold', bd=1, relief=SOLID,
               command=self._recalc, width=20).grid(row=2, column=2, sticky='wens', rowspan=2, padx=50)

        value_frame.columnconfigure(0, weight=1, minsize=110)
        value_frame.columnconfigure(1, weight=1, minsize=150)
        value_frame.columnconfigure(2, weight=1, minsize=1)

        button_frame = Frame(self.container_frame)
        button_frame.pack(side=TOP, anchor='w', fill=BOTH, expand=1, padx=15)

        back_button = Button(button_frame, text='Back', bd=1, relief=SOLID, command=lambda: self.back(), width=10)
        back_button.grid(row=0, column=0, sticky='we', padx=5)

        self.info_section = Frame(self.main_frame)
        self.info_section.grid(row=1, column=1, sticky='nswe')

        self.main_frame.columnconfigure(0, weight=1, minsize=500)
        self.main_frame.columnconfigure(1, weight=1, minsize=200)

        self.logs = None
        self.info = None

        self._init_output_frame()

        self.setted = False

    def _init_output_frame(self):
        if not TopBar.show_info.get() and not TopBar.show_logs.get():
            self.container_frame.grid(row=1, column=0, sticky='nsew', columnspan=2)
            self.info_section.grid_remove()
        else:
            self.container_frame.grid(row=1, column=0, sticky='nsew', columnspan=1)
            self.info_section.grid(row=1, column=1, sticky='nswe')

        if TopBar.show_logs.get():
            if not self.logs:
                self.logs = TextWidget(self.main_frame, self.info_section, 'Logs', 5, 45)
        else:
            if self.logs:
                self._del_out_frames()

        if TopBar.show_info.get():
            if not self.info:
                self.info = TextWidget(self.main_frame, self.info_section, 'Info', 5, 45)
        else:
            if self.info:
                self._del_out_frames()

    def _del_out_frames(self):
        for el in self.info_section.winfo_children():
            el.destroy()
        log.set_func(None)
        self.logs = None
        self.info = None
        self.update()

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def back(self):
        if self.prev_frame:
            self.main.show_frame(self.prev_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def _get_pstar(self, aic_type='aic', Kmin_corrfactor=1.0):
        cu.data.xf.cepstral_analysis(aic_type=aic_type, K_PSD=Kmin_corrfactor - 1)

    def _pstar(self):
        self.value_entry.config(from_=2, to=cu.data.xf.Nfreqs)
        self.value_entry.delete(0, END)
        self.value_entry.insert(0, (cu.data.xf.dct.aic_Kmin + 1))

        self.fstar_label.config(text='{:4f}'.format(cu.data.fstar))
        self.pstar_label.config(text=f'{cu.data.xf.dct.aic_Kmin + 1}')
        self.kmin_label.config(text='{:18f} +/- {:10f} W/mK'.format(cu.data.xf.kappa_Kmin, cu.data.xf.kappa_Kmin_std))

    def _change_increment(self):
        self.value_entry.config(increment=int(self.increment.get()))

    def _recalc(self):
        self._get_pstar(aic_type='aic', Kmin_corrfactor=int(self.value_entry.get()))
        self._get_pstar(aic_type='aic', Kmin_corrfactor=int(self.value_entry.get()))
        self.graph.add_graph(cu.gm.plot_cepstral_spectrum, 'cepstral', x=cu.data.xf)
        self.graph.update_cut()

    def _setup_pstar(self):
        cu.data.xf.cepstral_analysis(aic_type='aic', K_PSD=None)
        self._pstar()

    def update(self):
        super().update()

        if cu.data.fstar == cu.data.old_fstar:
            self.setted = True
        else:
            self.setted = False
            cu.data.old_fstar = cu.data.fstar

        if not self.setted:
            self.setted = True
            self._setup_pstar()

        self.graph.show(cu.gm.GUI_plot_periodogram, x=cu.data.j)
        self.graph.add_graph(cu.gm.resample_current, 'resample', x=cu.data.j, fstar_THz=cu.data.fstar,
                             PSD_FILTER_W=cu.data.psd_filter_width)

        self._init_output_frame()
        if self.info:
            cu.update_info(self.info)
        if self.logs:
            log.set_func(self.logs.write)
        self._recalc()
        self.graph.update_cut()
