
from func.tsodyks_meta import *
from func.helpers import *

J =[1.27]#400.,450.,450.,450.,600.]#[3822.]#[392.,427.,540.,854.,1709.,3822.]#[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]

eta =0.5#0.5
d =  [3.5]
NI_ =[200]#,200,500,800,990]#[40,200,800,3200,16000]#[50,200,500,800,950,990]#[200,500,800,950]#[200,800,2000,3200,3800]
U = 0.035
G = [4.0]#0.001,4,0.8,0.18,0.00101]
sim_time= 100000#110
tau_rec = 130000.
tau_fac = 0.0
directory   = 'sim/brunel_bic'
for bic in [0.9]:#@0.8,0.7,0.5,0.3,0.2,0.01]:#[0.8,0.5,0.2]: 
    for i,NI in enumerate(NI_):
        NE =int(1000-NI)#int(5000-NI)
        jA = J[i]
        j = jA
        g =(NE/NI)*bic#(NE/NI)*bic
        simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,g,eta,NI,U,tau_rec,tau_fac,NI+NE)
        A = meta_brunel(directory= directory,
                                    simulation = simulation,
                                    g=np.round(g,decimals=3), # inhibitor0y strenght
                                    eta = np.round(eta,decimals=3),
                                    d=d, # synaptic delay
                                    J=jA, #synaptic strengthd
                                    NE =NE, # fraction of inh neurons
                                    NI= NI,
                                    N_rec = NE+NI,
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
