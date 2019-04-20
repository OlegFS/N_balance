from func.brunel_meta import *
from func.helpers import *

#G =[19,4,1,0.25,0.05263157894736842]#[19,4,1,0.25,0.05]#[19,1,1.,0.3125,0.066,0.013]# [0.001,4,0.8,0.18,0.00101]#[4.,4.,0.875,0.195,0.04]#[4.,4.,1.0,0.2,0.052] 
eta =0.5#0.5
d   =[3.5]# ,100,200,300,500,700,800
NI_=[100]#,200,250,300,500,700,800]#[100,200,250,300,500,700,800,950]#,800,950]
J   =[1.19]#,1.27,1.3,1.33,1.56,1.9,2.0]#,2.1,2.5]#[100.,200.,200.,200.,200.,200.,240.,260.,250.]#[3822.]#[392.,427.,540.,854.,1709.,3822.]#[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]
U   =0.035

sim_time= 1000000#110
directory   = 'sim/brunel_balaned_match_cosyne_dt01/'#'sim/brunel_balaned/'
for i,NI in enumerate(NI_):
    NE =int(1000-NI)#int(5000-NI)
    jA = J[i]
    j = jA
    g =(NE/NI)#*20
    #if NI==10:
    #    g=g*0.001
    #print(int(NI*1.1))
    #g = (G[i]*1.25)#/bic
    #simulation ='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s'%(j,g,eta,NI,U,tau_rec,tau_fac)

    simulation='brunel_j=%s_g=%s_eta=%s_NI=%s_Nt=%s'%(j,g,eta,NI,NI+NE)
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
    A.build()
    A.connect()
    A.run()
    print(simulation)
    #A.get_connectivity()
