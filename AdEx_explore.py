# Analyze AdEx Behaviour 
from func.AdEx_meta import *
from func.helpers import *

%load_ext autoreload
import nest
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
#sns.set_context('talk')
from func.AdEx_meta import *
from func.helpers import *
from scipy.signal import find_peaks
na = np.array
plt.style.use('dark_background')

# %%

# plt.plot(SC[0])
J = 0.1

NI =2500#200#200
NE = int(12500-NI)#int(1000-NI)#
sim_time=500

params = [(3.,2), (6.0, 4),(5.0,2),(4.5,0.9)]#[(5.,2),(4.5,0.9), (1.0,2)]
titles = ['SR', 'SI_fast','AI','SI_slow']

bin_size = 0.1
plot_count = 1
plt.figure(figsize=(20,15))
path = 'sim/AdEx_test/'
for i,p in enumerate(params):

    # AI
    g,eta= p

    plt.subplot(4,2,plot_count)
    plt.title('AdEx')
    simulation ='AdEx_J=%s_g=%s_eta=%s_NI=%s_N_tot=%s_sim_time=%s'%(J,g,eta,NI,NI+NE,sim_time)
    plot_raster(path,simulation,t=(100,250),N=100,u_id=(0,100),marker=1)
    sc = return_sc(path,simulation,(0,sim_time),bin_size=bin_size)
    plt.plot(np.linspace(0,sim_time,len(sc)), 50*sc/np.max(sc),linewidth = 2,alpha =0.5)
    plot_count += 1

    plt.subplot(4,2,plot_count)
    plt.title('LIF')
    simulation ='LIF_J=%s_g=%s_eta=%s_NI=%s_N_tot=%s_sim_time=%s'%(J,g,eta,NI,NI+NE,sim_time)
    plot_raster(path,simulation,t=(100,250),N=100,u_id=(0,100),marker=1)
    sc = return_sc(path,simulation,(0,sim_time),bin_size=bin_size)
    plt.plot(np.linspace(0,sim_time,len(sc)),50*sc/np.max(sc),linewidth = 2,alpha =0.5)
    plot_count += 1
# plt.
plt.tight_layout()
# plt.savefig('figs/AdEx_LIF_comp_diagnostic.eps')