import nest
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
from func.brunel_meta import *
from func.helpers import *

# Check size effects in Brunlel' analytics

## Test if I get the Nu right 
J = [1.275]#[0.1,0.5,1.,1.275]
g = 4.0
eta = 0.5
d =  [3.5]
meanV = np.zeros(len(J))
stdV = np.zeros(len(J))
rate = np.zeros(len(J))
NE = 800
NI = 200
#
for ind,j in enumerate(J):
#    sim_time =5000#110
#    # Tprob = np.zeros([1,10000])
#    directory = 'sim/voltage_test'
    simulation = 'rate'
#    # for ind,i in #enumerate(np.arange(1,50,1)):
#    #     for tr in np.arange(10000):
#    A = meta_brunel(directory= directory,
#         simulation = simulation,
#         g=np.round(g,decimals=3), # inhibitor0y strenght
#         eta = np.round(eta,decimals=3),
#         d=d, # synaptic delay
#         J=j, #synaptic strength NE =800, # fraction of inh neurons
#         NE = NE,
#         NI= NI,
#         N_rec = NE+NI,
#         epsilon = 0.1,
#         simtime=sim_time,
#         master_seed = 1000,
#         verbose = False,
#         chunk= False,
#         chunk_size= 50000,
#        voltage = False)
#    A.build()
#    A.connect()
#    A.run() 
#    try:
#        print('no spikes')
#        sc = return_sc('sim/voltage_test/','rate-1002',(0,50000),N=NE+NI,bin_size = 1)
#        rate[ind] = np.sum(sc)/(NE+NI)/(sim_time/1000)
#    except:
#        pass
#    np.save('sim/voltage_test/sampled_rate_1000',rate)
#
#    sim_time =500#110
#    
    directory = 'sim/voltage_test'
    simulation = 'SecVoltTest_j=%s_n=%s'%(j,NE+NI)
#    # for ind,i in #enumerate(np.arange(1,50,1)):
#    #     for tr in np.arange(10000):
#    A = meta_brunel(directory= directory,
#         simulation = simulation,
#         g=np.round(g,decimals=3), # inhibitor0y strenght
#         eta = np.round(eta,decimals=3),
#         d=d, # synaptic delay
#         J=j, #synaptic strength NE =800, # fraction of inh neurons
#         NE = NE,
#         NI= NI,
#         N_rec = 1000,
#         epsilon = 0.1,
#         simtime=sim_time,
#         master_seed = 1000,
#         verbose = False,
#         chunk= False,
#         chunk_size= 50000,
#        voltage = True)
#    A.build()
#    A.connect()
#    A.run()
#
    s_voltage = read_voltage('sim/voltage_test/',
                         simulation,
                         1000,
                         (0,2180),
                         nice_format=True)

    #np.save('sim/voltage_test/sampled_voltage',s_voltage)
    mm = []
    sstd = []
    for i in range(0,s_voltage.shape[1]):
        mm.append(np.mean(s_voltage[:,i]))
        sstd.append(np.std(s_voltage[:,i]))
    meanV[ind] = np.mean(mm[5000:])
    stdV[ind] = np.mean(sstd[5000:])
    np.save('sim/voltage_test/sampled_meanV1',meanV)
    np.save('sim/voltage_test/sampled_stdV1',stdV)
