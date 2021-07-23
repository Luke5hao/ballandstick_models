from neuron import h 
#from neuron.units import ms, mV  
import matplotlib.pyplot as plt 
#from bokeh.io import output_notebook
#import bokeh.plotting as plt
#output_notebook()

ms = 1
mV = 1

h.load_file('stdrun.hoc')

class BallAndStick:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_biophysics()
    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.all = [self.soma, self.dend]
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        self.dend.L = 200
        self.dend.diam = 1
    def _setup_biophysics(self):
        for sec in self.all:
            sec.Ra = 100    # Axial resistance in Ohm * cm
            sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        self.soma.insert('hh')                                          
        for seg in self.soma:
            seg.hh.gnabar = 0.12  # Sodium conductance in S/cm2
            seg.hh.gkbar = 0.036  # Potassium conductance in S/cm2
            seg.hh.gl = 0.0003    # Leak conductance in S/cm2
            seg.hh.el = -54.3     # Reversal potential in mV
        # Insert passive current in the dendrite                  
        self.dend.insert('pas')                                  
        for seg in self.dend:                                       
            seg.pas.g = 0.001  # Passive conductance in S/cm2          
            seg.pas.e = -65    # Leak reversal potential mV            
    def __repr__(self):
        return 'BallAndStick[{}]'.format(self._gid)

my_cell = BallAndStick(0)

stim = h.IClamp(my_cell.dend(1))
stim.delay = 5
stim.dur = 1
stim.amp = 0.1

soma_v = h.Vector().record(my_cell.soma(0.5)._ref_v)
dend_v = h.Vector().record(my_cell.dend(0.5)._ref_v)
t = h.Vector().record(h._ref_t)

h.finitialize(-65 * mV)
h.continuerun(25 * ms)

plt.plot(t, soma_v, label = 'soma(0.5)')

amps = [0.075 * i for i in range(1, 5)]
colors = ['green', 'blue', 'red', 'black']
for amp, color in zip(amps, colors):
    stim.amp = amp
    for my_cell.dend.nseg, width in [(1, 2), (101, 1)]:
        h.finitialize(-65)
        h.continuerun(25)
        plt.plot(t, list(soma_v), linewidth = width, color=color)
        plt.plot(t, list(dend_v), linewidth = width, color=color, linestyle = 'dashed')
plt.show()

