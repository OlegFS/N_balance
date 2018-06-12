import nest
import nest.raster_plot
import time
from numpy import exp
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as sps
from itertools import product
from pathlib import Path
from func.brunel import fixed_brunel 
from func.brunel import save_dict
from func.brunel_meta import meta_brunel,disconnected_brunel, stim_brunel
from func.helpers import return_sc

# Transition probability for N stim neurons 
directory   = 'test'
simulation = 'init'

j = 0.85*1.5
g = 4.0
eta = 0.5
d =  [3.5]
sim_time = 110
Tprob = np.zeros([50,5000])
for ind,i in enumerate(np.arange(1,50,1)):
    for tr in np.arange(5000):
        A = stim_brunel(directory= directory,
             simulation = simulation,
             g=np.round(g,decimals=3), # inhibitor0y strenght
             eta = np.round(eta,decimals=3),
             d=d, # synaptic delay
             J=j, #synaptic strength NE =800, # fraction of inh neurons
             NI= 200,
             N_rec = 1000,
             epsilon = 0.1,
             simtime=sim_time,
             verbose = False,
             chunk= False,
             chunk_size= 50000,
            voltage = False)
        A.build(amp =10e16,init_voltage = True)
        A.connect(n_ext=i,Poisson=True)
        A.run()
        sc = return_sc('test/','init',(3.5,103.5),N=1000,bin_size = 1)
        Tprob[ind,tr] = np.sum(sc)
    np.save('test/StimTransitionProb_40_50',Tprob)


