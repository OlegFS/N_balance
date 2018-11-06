

from func.tsodyks_meta import *
from func.helpers import *

J =np.double(np.arange(10,1500,25)) #[ 450., 450., 450.,450.]#[400.,450.,450.,450.,600.]#[400.,400.,400.,400.,400.]
eta =0.1#0.5
d =  [3.5]
NE_ =[100,1000,4000,10000]
U = 0.035

sim_time= 150000#110
tau_rec = 100000.
tau_fac = 0.0
directory   = 'sim/tsodyks_Jscaling/'
for j in J:
    for i,NE in enumerate(NE_):

        g= 4
        #NE =#iint(1000)#int(5000-NI)
        NI = int(0.2*NE)

        simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,g,eta,NI,U,tau_rec,tau_fac,NI+NE)
        A = DeprSynapses(directory= directory,
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
                                    master_seed = 2000,
                                    verbose = True,
                                    chunk= False,
                                    chunk_size= 50000,
                                    voltage = False)
        A.build(tauSyn = 3.,constant_k=False)#rate = 180.)#196
        A.connect(U=U,u=1.,x=1., tau_rec=tau_rec ,tau_fac=tau_fac,tau_psc= 3.)
        A.run()
