import IPython.display as ipd
import numpy
from IPython.core.display import display
import time

def beep(frq = 440, sr = 22050,T = 0.5,):
    ''' make a noise when code finishes running'''
    t = numpy.linspace(0, T, int(T*sr), endpoint=False) # time variable
    x = 0.5*numpy.sin(2*numpy.pi*frq*t)              # pure sine wave at 440 Hz
    display(ipd.Audio(x, rate=sr, autoplay=True)) # load a NumPy array
    
def beep_repeat(delay=2, **kwargs):
    while True:
        beep(**kwargs)
        time.sleep(delay)


def beep_on_error(*args, frq=880, beep_func=beep_repeat, **kwargs):
    def on_error(self, etype, value, tb, tb_offset=None):
        self.showtraceback((etype, value, tb), tb_offset=tb_offset)
        beep_func(frq=frq, *args,**kwargs)
    get_ipython().set_custom_exc((Exception,), on_error)