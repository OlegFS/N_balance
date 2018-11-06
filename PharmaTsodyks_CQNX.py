from func.tsodyks_meta import *
from func.helpers import *

#G =[19,4,1,0.25,0.05263157894736842]#[19,4,1,0.25,0.05]#[19,1,1.,0.3125,0.066,0.013]# [0.001,4,0.8,0.18,0.00101]#[4.,4.,0.875,0.195,0.04]#[4.,4.,1.0,0.2,0.052] 
eta =0.2#0.5
d =  [3.5]
NI_ =[10,200,500,800,990]#[40,200,800,3200,16000]#[50,200,500,800,950,990]#[200,500,800,950]#[200,800,2000,3200,3800]
J =[400.,450.,450.,450.,600.]#[3822.]#[392.,427.,540.,854.,1709.,3822.]#[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]
U = 0.035
G = [0.001,4,0.8,0.18,0.00101]
sim_time= 300000#110
tau_rec = 120000.
tau_fac = 0.0
directory   = 'sim/tsodyks_bic'
for bic in [0.9,0.7,0.5,0.3,0.2]: 
    for i,NI in enumerate(NI_):
        NE =int(1000-NI)#int(5000-NI)
        jA = J[i]
        j = jA
        g = G[i]#NI/NE
        #g = (G[i]*1.25)#/bic
        #simulation ='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s'%(j,g,eta,NI,U,tau_rec,tau_fac)

        simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j*bic,g,eta,NI,U,tau_rec,tau_fac,NI+NE)
        A = DeprSynapses(directory= directory,
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
        A.connect(U=U,j_cnqx=j*bic,u=1.,x=1., tau_rec=tau_rec ,tau_fac=tau_fac,tau_psc= 3.)
        A.run()
