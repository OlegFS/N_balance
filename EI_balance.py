from func.tsodyks_meta import *
from func.helpers import *

J = [200.]#[300.,700.,700.,700.,700.,500.,510.,510.,450.]
eta =[0.2]#0.5
d =  [3.5]
NI =2000#[10,100,200,250,300,500,700,800,950]#[100,200,500,800,990]#[200,500,800,950]#[200,800,2000,3200,3800]
U = 0.035
G = [5.0]#,4,1,0.25]
sim_time =40000#10
tau_rec =130000.#260000.
tau_fac = 0.
directory   = 'sim/tsodyks_EIbalance/Reyes_protocol/tau/'
for n,g in enumerate(G):
    #for i in range(len(J)):
    for scale in [1]:#[0.2,0.5,0.75,1,1.25,1.5,2,4]:
        NE =int(10000-NI)
        jA = J[0]
        j = jA
        simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,g,eta[0],NI,U,tau_rec,tau_fac,NI+NE)
        #simulation='tsodyks_j=%s_eta=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,eta,U,tau_rec,tau_fac,NE)
        A = DeprSynapses(directory= directory,
                                    simulation = simulation,
                                    g=np.round(g,decimals=3), # inhibitor0y strenght
                                    eta = np.round(eta[0],decimals=3),
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
                                    voltage = True)

        A.build(U,tauSyn = 3.,constant_k=[0.1,0.1])#rate = 181.)#196
        A.connect(U=U,u=1.,x=1.,j_cnqx=False, tau_rec=tau_rec ,tau_fac=tau_fac,tau_psc= 3.)
        A.run(stim=1,stim_eta = 10.5,n_trials =1 )
        conn = A.get_connectivity()
