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
from func.brunel_meta import meta_brunel,disconnected_brunel,stim_brunel,PoisBrunel
from func.helpers import return_sc

# Transition probability for N stim neurons 
directory   = 'test'
simulation = 'init'


j = 0.5#*1.5
g = 5
#eta = 0.5
eta =  np.linspace(0.2,0.8,100)#[3.5]
d = [3.5]
sim_time = 100000#110
rate_ex = np.zeros(len(eta))
for ind,i in enumerate(eta):
	A = PoisBrunel(directory= directory,
				 simulation = simulation,
				 g=np.round(g,decimals=3), # inhibitor0y strenght
				 eta = np.round(i,decimals=3),
				 d=d, # synaptic delay
				 J=j, #synaptic strengthd
				 NE =1, # fraction of inh neurons
				 NI= 1,
				 N_rec = 2,
				 epsilon = 0.1,
				 simtime=sim_time,
				 verbose = True,
				 chunk= False,
				 chunk_size= 50000,
				voltage = False)
	A.build()#amp =0.,init_voltage = True)
	A.connect()#n_ext=1,Poisson=True)
	A.run()
	rate_ex[ind] = A.rate_ex
np.save('sim/rate_test/slow_rate_100_J0.5_sim.npy',rate_ex)
