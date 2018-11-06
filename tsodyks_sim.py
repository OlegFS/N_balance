from func.tsodyks_meta import *
from func.helpers import *

j =450.

jA = j#j/calJ#1.4
g =4.2#0.25
eta =0.2#0.5
d =  [3.5]
sim_time =3500000#110

NI =200
NE =int(1000-NI)
U = 0.035

tau_rec = 120000.
tau_fac = 0.0
directory   = 'sim/tsodyks/'
#simulation = 'disconnected'
simulation ='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,g,eta,NI,U,tau_rec,tau_fac,NI+NE)

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
                            master_seed = 1000,
                            verbose = True,
                            chunk= False,
                            chunk_size= 50000,
                            voltage = False)
A.build(tauSyn = 3.)#rate = 180.)#196
A.connect(U=U,u=1.,x=1., tau_rec=tau_rec ,tau_fac=tau_fac,tau_psc= 3.)
A.run()
print('Unit Amp', jA*U)
