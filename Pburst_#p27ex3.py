# Script for the second sampling experiment
import nest
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set_context('talk')
from func.brunel_meta import *
from func.helpers import *

#
#directory   = 'test/'
#simulation = 'monstrous'
#
#j =0.85*1.5
#g = 4
#eta = 0.5#0.5
#d =  [3.5]
#sim_time =5000000#110
#A = meta_brunel(directory= directory,
#                             simulation = simulation,
#                             g=np.round(g,decimals=3), # inhibitor0y strenght
#                             eta = np.round(eta,decimals=3),
#                             d=d, # synaptic delay
#                             J=j, #synaptic strengthd
#                             NE =800, # fraction of inh neurons
#                             NI= 200,
#                             N_rec = 1000,
#                             epsilon = 0.1,
#                             tauMem =40.,
#                             simtime=sim_time,
#                             master_seed = 2000,
#                             verbose = True,
#                             chunk= False,
#                             chunk_size= 50000,
#                             voltage = False)
#A.build()#rate = 180.)#196
#A.connect()
#A.run()
#conn = A.get_connectivity()
## Transition probability for N stim neurons 
directory   = 'test'
simulation = 'init'
#
j = 0.85*1.5
g = 4.0
eta = 0.5
d =  [3.5]
sim_time = 110
repeat = 5000
x = 40
Tprob = np.zeros([x,repeat])
RealTsp_e = np.zeros([x,repeat])
RealTsp_i = np.zeros([x,repeat])
for ind,i in enumerate(np.arange(1,x,1)):
    for tr in np.arange(repeat):
        A = stim_brunel(directory= directory,
             simulation = simulation,
             g=np.round(g,decimals=3), # inhibitor0y strenght
             eta = np.round(eta,decimals=3),
             d=d, # synaptic delay
             J=j, #synaptic strength NE =800, # fraction of inh neurons
             NE = 800,
             NI= 200,
             N_rec = 1000,
             epsilon = 0.1,
             simtime=sim_time,
             master_seed = 2000,
             verbose = False,
             chunk= False,
             chunk_size= 50000,
            voltage = False)
        A.build(amp =10e16,init_voltage = True,nu = 0.137)
        A.connect(n_ext=i,Poisson=True)
        A.run()
        sc = return_sc('test/','init',(3.5,103.5),N=1000,bin_size = 1)
        Tprob[ind,tr] = np.sum(sc)
        sc = return_sc('test/','init',(0,100),bin_size=0.5,cells = 'ex')
        RealTsp_e[ind,tr] = sc[6]
        sc = return_sc('test/','init',(0,100),bin_size=0.5,cells = 'in')
        RealTsp_i[ind,tr] = sc[6]
    np.save('test/StimTransitionProb_ex3_fixedNu',Tprob)
    np.save('test/StimTransitionProb_ex3_1_40_realScE',Tprob)
    np.save('test/StimTransitionProb_ex3_1_40_realScI',Tprob)

