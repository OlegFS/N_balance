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
from func.brunel_meta import meta_brunel,fixedOutDeg_brunel
#Load edges
#edges = np.load('Fixed_in_out_edges.npy')
g_=[7]# np.round(np.arange(2,4,1),2).tolist()#[10.0,10.5,11,11.5,12.0,12.5,13.0]#[4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,9.0,9.5]# [4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5,5.6,5.7,5.8,5.9,6.0]#
eta_ = [0.50]
d =[3.5]
j_ =np.round(np.arange(0.8,1.5,0.01),2).tolist()#[0.81,0.82,0.83,0.84,0.85,0.86,0.87,0.88,0.89,0.9]
values =  list(product(g_,j_))
repeat  =False 
NE =8000#2400 
NI = 2000
sim_time = 400000
for i in range(len(values)):
    g = values[i][0]
    j = values[i][1]
    eta = eta_[0]#values[0][1]
    print('delta %s, %s'%(values[i]))
    simulation='brunel_esp_ex_%s_g_%s_delta_%s-10000'%(g,eta,j)
    directory = 'sim/Jscale/N10000'
    s = Path(directory+'/'+simulation)
    if (not s.exists()) and (not s.exists()) or repeat ==True:
        A = meta_brunel(directory= directory,
             simulation = simulation,
             g=np.round(g,decimals=3), # inhibitor0y strenght
             eta = np.round(eta,decimals=3),
             d=d, # synaptic delay
             J=j, #synaptic strength
             NE =NE, # fraction of inh neurons
             NI= NI,
             N_rec = NE+NI,
             epsilon = 0.1,
             simtime=sim_time,
             verbose = True,
             chunk= False,
             chunk_size= 50000,
             voltage = False)
        A.build()
        A.connect()
        A.run()
        #_ = A.get_Nout_deg()
        print(not s.exists()) 
        print(simulation)
        print('Simulation delta %s, %s'%(values[i]))
    else:
        print('exists')

