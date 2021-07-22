"""from neuron import h 
#from neuron.units import ms, mV  

class BallAndStick:
    def __init__(self, gid):
        self.agid = gid
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.dend.connect(self.soma)
    def __repr__(self):
        return 'BallAndStick[{}]'.format(self.agid)
         
my_cell = BallAndStick(0)"""
