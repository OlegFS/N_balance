from func.brunel_meta import *
from func.helpers import *

J =1.9#1.45#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]
G =[4.1]# 5,2.5,1.5625, 1.316]#[0.001,4,0.8,0.18,0.00101]#[4.,4.,0.875,0.195,0.04]#[4.,4.,1.0,0.2,0.052] 
eta =0.4#0.5
d =  [3.5]
NI_ =[0.25]#, 0.25, 0.5,0.8,0.95]#, 0.25,0.5,0.8,0.99]#[200,800,2000,3200,3800]

sim_time =35000#110
directory   = 'sim/burts_exponents/'
for i,NI_f in enumerate(NI_):
    
    NE =800#int(1000-NI)
    j = J
    g = G[i]
    NI  =int(NE*NI_f) 
    #simulation ='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s'%(j,g,eta,NI,U,tau_rec,tau_fac)

    simulation ='brunel_j=%s_g=%s_eta=%s_NI=%s_Nt=%s'%(j,g,eta,NI,NI+NE) 
    A = meta_brunel(directory= directory,
                                simulation = simulation,
                                g=np.round(g,decimals=3), # inhibitor0y strenght
                                eta = np.round(eta,decimals=3),
                                d=d, # synaptic delay
                                J=j, #synaptic strengthd
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
    A.build()
    A.connect()
    A.run()
