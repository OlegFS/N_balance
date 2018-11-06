

from func.tsodyks_meta import *
from func.helpers import *

J =[250.]#[10.,100.,400.,800.,1000.]#[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]
G =[19]#[4,1,0.25,0.053]# [0.001,4,0.8,0.18,0.00101]#[4.,4.,0.875,0.195,0.04]#[4.,4.,1.0,0.2,0.052] 
eta =0.1#0.5
d =  [3.5]
NI_ =[5]#[200,500,800,950]#[200,800,2000,3200,3800]
U = 0.035

sim_time =300000#110
tau_rec =180000.
tau_fac = 0.0
directory   = 'sim/tsodyks_exc_singleTau/'
for i in range(len(J)):
    NE =4000#int(5000-NI)
    jA = J[i]
    j = jA
    g = 19.0
    NI = 10
    #simulation ='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s'%(j,g,eta,NI,U,tau_rec,tau_fac)

    simulation='tsodyks_j=%s_eta=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s_constantK'%(j,eta,U,tau_rec,tau_fac,NE)
    A = ExNetDeprS(directory= directory,
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
    A.build(U,tauSyn = 3.,constant_k=False)#rate = 180.)#196
    A.connect(U=U,u=1.,x=1., tau_rec=tau_rec ,tau_fac=tau_fac,tau_psc= 3.)
    A.run()
