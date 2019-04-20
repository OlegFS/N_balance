

from func.tsodyks_meta import *
from func.helpers import *
def EIratio(J,Ne,Ni,pe,pi,g):
#     print('Ai',(Ki*p)*(g*J))
#     print('Ae',(Ke*p)*(J))
    print('E/I',((Ne*pe)*(J))/((Ni*pi)*(g*J)))
    return ((Ne*pe)*(J))/((Ni*pi)*(g*J))

#G =[19,4,1,0.25,0.05263157894736842]#[19,4,1,0.25,0.05]#[19,1,1.,0.3125,0.066,0.013]# [0.001,4,0.8,0.18,0.00101]#[4.,4.,0.875,0.195,0.04]#[4.,4.,1.0,0.2,0.052] 
eta =0.2#0.5
d   =[3.5]# ,100,200,300,500,700,800
NI_=[100,200,250,300,500,700,800,950]
J   =[200.]*8#,200.,200.,200.,200.,200.,240.,260.,250.]#[3822.]#[392.,427.,540.,854.,1709.,3822.]#[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]
U   =0.035

sim_time= 500000#110
tau_rec = 30000.
tau_fac = 0.0
directory   = 'sim/tsodyks_Kbalaned_match_cosyne_dt01'
for i,NI in enumerate(NI_):
    NE =int(1000-NI)#int(5000-NI)
    jA = J[i]
    j = jA
    g = 4
    #if NI==10:
    #    g=g*0.001
    #print(int(NI*1.1))
    #g = (G[i]*1.25)#/bic
    #simulation ='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s'%(j,g,eta,NI,U,tau_rec,tau_fac)
    K = [80/NE , 20/NI] 

    if EIratio(1,NE,NI,K[0],K[1],4)!=1:
        print('Warning')
        break
    if NI==10:
        K[1]=0.001
    elif NI==950:
        K[0] = 0.8
    simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Ke=%s_Nt=%s'%(j,g,eta,NI,U,tau_rec,tau_fac,K[0],NI+NE)
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
    A.build(U=U, tauSyn = 3.,constant_k=K)#rate = 180.)#196
    A.connect(U=U,j_cnqx=False,u=1.,x=1., tau_rec=tau_rec ,tau_fac=tau_fac,tau_psc= 3.)
    A.run()
