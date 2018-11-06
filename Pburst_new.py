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
##
#directory   = 'test/'
#simulation = 'hyper'
#
#j =0.85*1.5
#g = 7
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
### Transition probability for N stim neurons 
directory   = 'P_burst'
#
j = 1.8
g = 4.1
eta = 0.4
d =  [3.5]
sim_time = 140
repeat =1000 
n_st1 = 9
n_st2 = 10
Tprob = np.zeros([n_st2,repeat])
RealTsp_e = np.zeros([n_st2,repeat])
RealTsp_i = np.zeros([n_st2,repeat])
simulation ='brunel_j=%s_g=%s_eta=%s_NI=%s_Nt=%s'%(j,g,eta,200,1000)
for ind,i in enumerate(np.arange(n_st1,n_st2,1)):
#for tr in np.arange(repeat):
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
    A.build(amp =0.,c_start = 0.,init_voltage = True,nu = 0.137)#, n_st =i)
    A.connect(n_ext=0,Poisson=True)
    A.run(repeat = 15,n_st = i,prime = True)
#    sc = return_sc(directory,simulation,(3.5,103.5),bin_size = 1)
#    Tprob[ind,tr] = np.sum(sc)
#    sc = return_sc(directory,simulation,(0,10),bin_size=0.5,cells = 'ex')
#    RealTsp_e[ind,tr] = sc[0]
#    sc = return_sc(directory,simulation(0,10),bin_size=0.5,cells = 'in')
#    RealTsp_i[ind,tr] = sc[0]
#    np.save(directory+'/'+'Tprob',Tprob)
#    np.save(directory+'/'+'Tprob_realScE_random_sl',RealTsp_e)
#    np.save(directory+'/'+'Tprob_realScI_random_sl',RealTsp_i)
#
