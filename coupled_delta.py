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
from func.brunel import save_dict

g_=[10.0,10.5,11,11.5,12.0,12.5,13.0]#[4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,9.5]# [4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0]#
eta_ = [0.50]
d =[3.5]
j_ =[0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9]
values =  list(product(g_,j_))
repeat  = False 
sim_time = 200000
for i in range(len(values)):
    g = values[i][0]
    j = values[i][1]
    eta = eta_[0]#values[0][1]
    print('delta %s, %s'%(values[i]))
    simulation=Path('sim/delta_coupled/tau40/brunel_esp_ex_%s_g_%s_test01dt_delta_%s.hdf5'%(g,eta,j))
    simulation_=Path('sim/delta_coupled/tau40/brunel_esp_ex_%s_g_%s_test01dt_delta_%s.npy'%(g,eta,j))
    if (not simulation.exists()) and (not simulation_.exists()) or repeat ==True:
        print(not simulation.exists()) 
        print(simulation)
        print('Simulation delta %s, %s'%(values[i]))
        ispikes,espikes= delta_brunel(g=np.round(g,decimals=3), # inhibitor0y strenght
                 eta = np.round(eta,decimals=3),
                 d=d, # synaptic delay
                 J=j, #synaptic strength
                 NE =8000, # fraction of inh neurons
                 NI= 2000,
                 N_rec = 8000,
                 epsilon = 0.1,
                 tauMem= 40.0,
                 CMem = 250.0,
                # tauSyn=0.5,
                 simtime=sim_time,
                 verbose = True,
                 chunk= True,
                 chunk_size= 50000)
        save_dict(espikes, ispikes,simulation)
        #n_events = nest.GetStatus(espikes,"n_events")[0]
        #n_events_i = nest.GetStatus(ispikes,"n_events")[0]
        #ts, gids = nest.raster_plot._from_memory(espikes)
        #ts_i, gids_i = nest.raster_plot._from_memory(ispikes)
        #dictionary = {'ts':ts,
        #              'gids':gids,
        #              'ts_i':ts_i,
        #              'gids_i':gids_i}
        #np.save(simulation, dictionary) 
    else:
        print('exists')

