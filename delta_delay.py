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
from func.brunel import delta_brunel 
from func.brunel import alpha_brunel 

g_= [4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0]
eta_ = [0.85]
d =[2.5,9.5]
values =  list(product(g_,eta_))
repeat = False
sim_time = 60000

#g = [3.6,3.8,38]#[3.8]#[30,3.8,0.618,10e-10]#[2000,2,.5,.1]
#50
for i in range(len(values)):
    g = values[i][0]
    eta = values[i][1]
    print('delta %s'%(i))
    simulation =Path('sim/alpha_delta/brunel_esp_ex_%s_g_%s_test01dt_delta_delay_2595.npy'%(g,eta))
    if not  simulation.exists() or repeat ==True:
        ispikes,espikes= delta_brunel(g=np.round(g,decimals=3), # inhibitor0y strenght
                 eta = np.round(eta,decimals=3),
                 d=d, # synaptic delay
                 J=0.1, #synaptic strength
                 NE =8000, # fraction of inh neurons
                 NI= 2000,
                 N_rec = 8000,
                 epsilon = 0.1,
                 tauMem= 20.0,
                 CMem = 250.0,
                # tauSyn=0.5,
                 simtime=sim_time,
                 verbose = True)

        n_events = nest.GetStatus(espikes,"n_events")[0]
        n_events_i = nest.GetStatus(ispikes,"n_events")[0]
        ts, gids = nest.raster_plot._from_memory(espikes)
        dictionary = {'ts':ts,'gids':gids}
        np.save(simulation, dictionary) 
    else:
        print('exists')
