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
from func.brunel_meta import meta_brunel,fixedOutDeg_brunel,fixed_brunel,disconnected_brunel,ExcitatoryBrunel
#Load edges
#edges = np.load('fixed_1000.npy')
g_=[5.0]#[10.0,10.5,11,11.5,12.0,12.5,13.0]#[4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,9.5]# [4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0]#
eta_ = [0.50]
d =[3.5]
j_ =[0.85*1.5]#*1.5]#[0.81,0.82,0.8300,0.84,0.85,0.86,0.87,0.88,0.89,0.9]
values =  list(product(g_,j_))
repeat  =True 
sim_time = 400000
for i in range(len(values)):
    g = values[i][0]
    j = values[i][1]
    eta = eta_[0]#values[0][1]
    print('delta %s, %s'%(values[i]))
    simulation='brunel_esp_ex_%s_g_%s_test01dt_delta_%s-1000-'%(g,eta,j)
    directory = 'out_deg'
    s = Path(directory+'/'+simulation)
    if (not s.exists()) and (not s.exists()) or repeat ==True:
        A =meta_brunel(directory= directory,
             simulation = simulation,
             g=np.round(g,decimals=3), # inhibitor0y strenght
             eta = np.round(eta,decimals=3),
             d=d, # synaptic delay
             J=j, #synaptic strength
             NE =800, # fraction of inh neurons
             NI= 200,
             N_rec = 1000,
             epsilon = 0.1,
             simtime=sim_time,
             verbose = True,
             chunk= False,
             chunk_size= 10000,
             voltage = True)
        A.build()
        A.connect()#edge = edges)
        A.run()
        #_ = A.get_Nout_deg()
        #_ = A.get_Nin_deg()
        _ = A.get_connectivity()
        print(not s.exists()) 
        print(simulation)
        print('Simulation delta %s, %s'%(values[i]))
    else:
        print('exists')

