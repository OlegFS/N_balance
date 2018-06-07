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

g_=[7.5]# [4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0]#
eta_ = [0.50]
d =[3.5]
j = 0.9
values =  list(product(g_,eta_))
repeat = True
sim_time = 200000
for i in range(len(values)):
    print('alpha %s'%(i))
    g = values[i][0]
    eta = values[i][1]
    simulation=Path('sim/alpha_coupled/brunel_esp_ex_%s_g_%s_test01dt_alpha_coupled_%s_long.npy'%(g,eta,j))
    if not  simulation.exists() or repeat ==True:
        ispikes,espikes= alpha_brunel(g=np.round(g,decimals=3), # inhibitor0y strenght
                 eta = np.round(eta,decimals=3),
                 d=d[0], # synaptic delay
                 J=j, #synaptic strength
                 NE =8000, # fraction of inh neurons
                 NI= 2000,
                 N_rec = 8000,
                 epsilon = 0.1,
                 tauMem= 20.0,
                 CMem = 250.0,
                 tauSyn=0.5,
                 simtime=sim_time,
                 verbose = True)

        n_events = nest.GetStatus(espikes,"n_events")[0]
        n_events_i = nest.GetStatus(ispikes,"n_events")[0]
        ts, gids = nest.raster_plot._from_memory(espikes)
        dictionary = {'ts':ts,'gids':gids}
        np.save(simulation, dictionary) 
    else:
        print('exists')
