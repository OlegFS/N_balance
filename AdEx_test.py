from func.AdEx_meta import *
from func.helpers import *
import sys
J =[1.4]#,700.,700.,700.,700.,500.,510.,510.,450.]#400.,450.,450.,450.,600.]#[3822.]#[392.,427.,540.,854.,1709.,3822.]#[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]

g =5.#np.float(sys.argv[1])#4.0#(NE/NI)
eta =0.5#np.float(sys.argv[2])#0.5

d =  [3.5]
NI =200#,100,200,250,300,500,700,800,950]#,200,500,800,990]#[40,200,800,3200,16000]#[50,200,500,800,950,990]#[200,500,800,950]#[200,800,2000,3200,3800]
sim_time= 10000#110
directory   = 'sim/AdEx_test'
NE =int(1000-NI)#int(5000-NI)

jA = J[0]
j = jA
simulation ='AdEx_J=%s_g=%s_eta=%s_NI=%s_N_tot=%s_sim_time=%s'%(J[0],g,eta,NI,NI+NE,sim_time)
#simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,g,eta,NI,U,tau_rec,tau_fac,NI+NE)
A = AdExModel(directory= directory,
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
                        master_seed = 2000,
                        verbose = True,
                        chunk= False,
                        chunk_size= 50000,
                        voltage = False)
A.build(tau_w = 100000., a=0.,b = 0., constant_k=[0.1,0.1])#rate = 180.)#196
A.connect(j_cnqx=False)
A.run()
#simulation = 'comp'
simulation ='LIF_J=%s_g=%s_eta=%s_NI=%s_N_tot=%s_sim_time=%s'%(J[0],g,eta,NI,NI+NE,sim_time)
B = meta_brunel(directory= directory,
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
B.build()#tau_w = 100., a=0.,b = 0., constant_k=[0.1,0.1])#rate = 180.)#196
B.connect()#(j_cnqx=False)
B.run()
