# Script for matching IBI with real data 
import nest
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set_context('talk')
from func.brunel_meta import *
from func.helpers import *


directory   = 'sim/match/'

j =0.1#1.27
g = 4
eta = 1.0#0.5
d =  [3.5]
sim_time =2000#110
simulation = 'match_ex1_N= %s_g=%s_j=%s'%(1000,g,j)
A = meta_brunel(directory= directory,
                             simulation = simulation,
                             g=np.round(g,decimals=3), # inhibitor0y strenght
                             eta = np.round(eta,decimals=3),
                             d=d, # synaptic delay
                             J=j, #synaptic strengthd
                             NE =800, # fraction of inh neurons
                             NI= 200,
                             N_rec = 1000,
                             epsilon = 0.1,
                             tauMem =40.,
                             simtime=sim_time,
                             master_seed = 1000,
                             verbose = True,
                             chunk= False,
                             chunk_size= 50000,
                             voltage = False)
A.build()#rate = 180.)#196
A.connect()
A.run()
conn = A.get_connectivity()
