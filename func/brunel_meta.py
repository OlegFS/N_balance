import sys
sys.path.append('/home/nest/nest-install/lib/python3.6/site-packages')

import nest
import nest.raster_plot
import time
from numpy import exp
import numpy as np
import h5py as h5py
import os 

class ChunkError(Exception): pass
class BuildError(Exception): pass

class meta_brunel(object):
    """ Network of LIF neurons. 
    inspired by Gewaltig et al. 2012"""
    threads = 4
    built = False
    connected = False
    dt      = float(0.1)    # the resolution in ms for the clock-based
    CMem =1.0
    V_m = 0.0
    V_res = 10.0#10
    t_ref = 2.0
    theta = 20.0

    def __init__(self,
                simulation = 'test',
                directory = 'test',
                g=7, # inhibitory strenght
                eta = 3.5, #Poisson rate ratio
                d=[3.5], # synaptic delay
                J=0.1, #synaptic strength
                NE =8000, # fraction of inh neurons
                NI =2000,
                N_rec = 50,
                epsilon = 0.1,
                tauMem= 20.0,
                simtime=1000,
                master_seed = 1000,
                verbose = True,
                chunk = False,
                chunk_size = 10,
                voltage= False):
        """Initialize the simulation , setup data directory"""
        self.name=self.__class__.__name__
        self.data_path=directory
        self.simulation = simulation
        nest.ResetKernel()
        if verbose ==False:
            nest.set_verbosity("M_ALL")
        if not os.path.exists(self.data_path):
            os.makedirs(self.data_path)
            print('”Writing data to :'+self.data_path)
        
        nest.SetKernelStatus({'resolution': self.dt, 'print_time': True,
                              'overwrite_files':True,
                              'local_num_threads':self.threads,
                              'data_path': self.data_path})
        self.g = g
        self.eta = eta
        self.delay = d
        self.J = J
        self.NE =int(NE)
        self.NI = int(NI)
        self.N_neurons = NE + NI
        self.epsilon = epsilon
        self.CE    = int(epsilon*NE) # number of excitatory synapses per neuron
        self.CI    = int(epsilon*NI) # number of inhibitory synapses per neuron  
        self.C_tot = int(self.CI+self.CE)
        self.N_rec = N_rec
        self.tauMem = tauMem
        self.simtime = simtime
        self.verbose = verbose
        self.chunk = chunk
        self.chunk_size = chunk_size
        self.record_vol = voltage
        #Some statistics
        self.in_deg_saved = False
        self.out_deg_saved = False
        if verbose:
            #neurons output
            print('Number of neurons (E/I):')
            print(NE,NI)
            print('Connection Probability: %s'%(epsilon))

        #Set the RNG seeds for each thread
        self.master_seed = master_seed#2000#1000
        self.n_vp = nest.GetKernelStatus('total_num_virtual_procs')
        print(self.n_vp)
        self.msdrange1 = range(self.master_seed, self.master_seed+self.n_vp)
        self.pyrngs = [ np.random.RandomState(s) for s in self.msdrange1]
        self.msdrange2 = range(self.master_seed+self.n_vp+1, self.master_seed+1+2*self.n_vp)
        nest.SetKernelStatus({'grng_seed': self.master_seed+self.n_vp,
                              'rng_seeds' :self.msdrange2}) 
        
    def build(self):
        
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                         "tau_m":self.tauMem,
                         "t_ref":self.t_ref,
                         "E_L": 0.0,
                         "V_reset": self.V_res,
                         "V_m": self.V_m,
                         "V_th":self.theta,
                         "refractory_input": False}
        if self.verbose:
            print(self.neuron_params)


        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        self.nu_th =self.theta / (self.J * self.CE * self.tauMem)# (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex * self.CE
        print(self.p_rate, self.nu_ex)        
        nest.SetDefaults("iaf_psc_delta", self.neuron_params)
        nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("iaf_psc_delta",self.NE)
        self.nodes_in = nest.Create("iaf_psc_delta",self.NI)
        self.noise    = nest.Create("poisson_generator")
        self.espikes  = nest.Create("spike_detector")
        if self.record_vol:
            self.voltmeter = nest.Create("multimeter")
        #self.ispikes  = nest.Create("spike_detector")
        self.nodes_al = self.nodes_ex +self.nodes_in
        nest.SetStatus(self.espikes,[{"label": self.simulation,
                                 "withtime": True,
                                 "withgid": True,#%(n_i_),
                                 "to_file": True}])
        if self.record_vol:
            nest.SetStatus(self.voltmeter,[{"label": self.simulation+'_voltage',
                                    'record_from': ['V_m'],
                                    'interval':self.dt,
                                     "withtime": True,
                                     "withgid": True,#%(n_i_),
                                     "to_file": True}])
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True
        self.built = True
    def connect(self):
        print(nest.GetStatus([1],'C_m'))
        """Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})

        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                   }
            self.syn_dict_in = {'model': 'inhibitory',
                    'delay':self.delay,
                   }
        else:
            self.d_min = self.delay[0]#3.5
            self.d_max = self.delay[1]#5.5
            self.syn_dict = {'model': 'excitatory',
                    'delay':{'distribution': 'uniform', 'low': self.d_min,
                             'high': self.d_max},
                   }
            print('uniform delay')

            self.syn_dict_in = {'model': 'inhibitory',
                    'delay': {'distribution': 'uniform', 'low': self.d_min,
                              'high': self.d_max},
                   }
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict)
        nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_al[:self.N_rec])
        print("Connecting network")

        print("Excitatory connections")


        self.conn_params_ex = {'rule': 'fixed_indegree', 'indegree': self.CE}
        nest.Connect(self.nodes_ex, self.nodes_ex+self.nodes_in,
                     self.conn_params_ex, syn_spec =self.syn_dict)

        print("Inhibitory connections")

        self.conn_params_in = {'rule': 'fixed_indegree', 'indegree': self.CI}
        nest.Connect(self.nodes_in, self.nodes_ex+self.nodes_in,
                     self.conn_params_in, syn_spec =self.syn_dict_in)

    def run(self,stim=0,stim_eta = 0.5,stim_tim = 5000,n_trials = 1):
        """run the simulation
        Chunk - to prevent from running unnecessary
        stim - with extra Poiss inoput"""
        print("Simulating")
        # New option to simulate in chunk and check a special condition in between
        # Thus we avoid unwanted condition
        if self.chunk:
            self.chunk_times = np.arange(0,self.simtime+1,self.chunk_size)[1::]
            if len(self.chunk_times) != int(self.simtime/self.chunk_size):
               raise ChunkError('check the simulation length') 
            for ch in self.chunk_times: 
                nest.Simulate(self.chunk_size)
                self.events_ex = nest.GetStatus(self.espikes,"n_events")[0]
                self.rate_ex   = self.events_ex/self.ch*1000.0/self.N_rec
                if self.rate_ex>60.0:
                    print('Rate is too high. Finishing simulation at %s ms'%(self.ch))
                    self.simtime = self.ch
                    break
        elif stim==1:
            self.nu_th= self.theta / (self.J * self.CE * self.tauMem)## (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
            self.new_nu_ex = stim_eta * self.nu_th
            self.new_p_rate = 1000.0 * self.new_nu_ex  * self.CE
            
            
            self.add_noise    = nest.Create("poisson_generator")
            nest.Simulate(self.simtime)
            nest.SetStatus(self.add_noise,{"rate":self.new_p_rate})
            nest.Connect(self.add_noise,self.nodes_ex[:100], syn_spec=self.syn_dict)
            nest.Connect(self.add_noise,self.nodes_in[:100], syn_spec=self.syn_dict)
            for trial in range(n_trials):   
                nest.SetStatus(self.add_noise,{"rate":self.new_p_rate})
                nest.Simulate(stim_tim)
                nest.SetStatus(self.add_noise,{"rate":0.0})
                nest.Simulate(10000)
        else:
            nest.Simulate(self.simtime)
        self.events_ex = nest.GetStatus(self.espikes,"n_events")[0]
        #self.events_in = nest.GetStatus(self.ispikes,"n_events")[0]

        self.rate_ex   = self.events_ex/self.simtime*1000.0/self.N_rec
        #self.rate_in   = self.events_in/self.simtime*1000.0/self.N_rec


        self.num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
                       nest.GetDefaults("inhibitory")["num_connections"])

        if self.verbose:
            print("Brunel network simulation (Python)")
            print("Number of neurons : {0}".format(self.N_neurons))
            print("Number of synapses: {0}".format(self.num_synapses))
            print("       Exitatory  : {0}".format(int(self.CE *self.N_neurons)
                                                   +self.N_neurons))
            print("       Inhibitory : {0}".format(int(self.CI * self.N_neurons)))
            print("Rate   : %.2f Hz" % self.rate_ex)
            #print("Inhibitory rate   : %.2f Hz" % self.rate_in)
            print("Simulation time   : %.2f s" % self.simtime)

    def get_Nin_deg(self):
        self.in_deg = np.zeros(self.N_rec)
        #out_deg = []#
        if self.in_deg_saved == False:
            for indx,i in enumerate(np.arange(1,self.N_rec+1)):
                conn = nest.GetConnections(target=[i])
                self.out_deg[indx] = len(nest.GetStatus(conn))
                #out_deg.append(nest.GetStatus(conn))
            with h5py.File(self.data_path+'/'+self.simulation+'-in_deg','w') as f:
                f['in_deg'] = self.in_deg
            print('saved as hdf5')
            self.in_deg_saved = True
        else:
            with h5py.File(self.data_path+'/'+self.simulation+'-in_deg','r') as f:
                self.in_deg=np.array(f['in_deg']) 
        return self.in_deg 

    def get_Nout_deg(self):
        self.out_deg = np.zeros(self.N_rec)
        #out_deg = []#
        if self.out_deg_saved == False:
            for indx,i in enumerate(np.arange(1,self.N_rec+1)):
                conn = nest.GetConnections(source=[i])
                self.out_deg[indx] = len(nest.GetStatus(conn))
                #out_deg.append(nest.GetStatus(conn))
            with h5py.File(self.data_path+'/'+self.simulation+'-out_deg','w') as f:
                f['out_deg'] = self.out_deg
            print('saved as hdf5')
            self.out_deg_saved = True
        else:
            with h5py.File(self.data_path+'/'+simulation+'-out_deg','r') as f:
                self.out_deg=np.array(f['out_deg']) 
        return self.out_deg 
    def get_connectivity(self):
        self.conn = nest.GetConnections()#source = range(1,self.N_rec+1)
        self.conn = np.array(self.conn)[:,0:2]
        with h5py.File(self.data_path+'/'+self.simulation+'-conn','w') as f:
            f['conn'] = self.conn
        print('saved as hdf5')
        return self.conn

class fixed_brunel(meta_brunel):
    """ Add fixed indegree-outdegree connectivity"""
    def connect(self,edge=None):
        """Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')
        if edge is None:
            raise BuildError('Need to provide the edges,(TODO: automatic build)')
        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})
        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                   }
            self.syn_dict_in = {'model': 'inhibitory',
                    'delay':self.delay,
                   }
        else:
            self.d_min = self.delay[0]#3.5
            self.d_max = self.delay[1]#5.5
            self.syn_dict = {'model': 'excitatory',
                    'delay':{'distribution': 'uniform', 'low': self.d_min,
                             'high': self.d_max},
                   }
            print('uniform delay')

            self.syn_dict_in = {'model': 'inhibitory',
                    'delay': {'distribution': 'uniform', 'low': self.d_min,
                              'high': self.d_max},
                   }
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict)
        nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)
        nest.Connect(self.nodes_al[:self.N_rec], self.espikes, syn_spec=self.syn_dict)

        print("Connecting network")

        print("Excitatory connections: building from the edges.Might take a while...")
        # TODO: compute the number of connections per populations
        self.ex_syn =  int(self.CE *self.N_neurons)
        print('number of excitatory syn (in=out): %s '%self.ex_syn)
        #nest.Connect(edge[:self.ex_syn][:,0], edge[:self.ex_syn][:,1], syn_spec =self.syn_dict)
        #ind_nodes = np.arange(0,self.NE*self.CE,self.CE)
        for n in edge[0:80000]:
            nest.Connect([n[0]],
                         [n[1]],syn_spec =self.syn_dict)
        #for n in range(800):
        #    nest.Connect(edge[ind_nodes[n]:ind_nodes[n+1],0].tolist(),
        #                 edge[ind_nodes[n]:ind_nodes[n+1],1].tolist(), syn_spec =self.syn_dict)
        print("Inhibitory connections")

        for n in edge[80000:]:
            nest.Connect([n[0]],
                         [n[1]],syn_spec =self.syn_dict_in)
        #nest.Connect(edges[self.ex_syn:][:,0], edges[self.ex_syn:][:,1], syn_spec =self.syn_dict)

        #ind_nodes = np.arange(self.NE*self.CE,self.NE*self.CE+self.NI*self.CI,self.CI)
        #for n in range(2):
        #    nest.Connect(edge[ind_nodes[n]:ind_nodes[n+1],0].tolist(),
        #                 edge[ind_nodes[n]:ind_nodes[n+1],1].tolist(), syn_spec=self.syn_dict_in)
        #for edge in edges[self.ex_syn:]:
        #    nest.Connect([edge[0]], [edge[1]], syn_spec =self.syn_dict_in)


class fixedOutDeg_brunel(meta_brunel):
    def connect(self):
        """Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})
        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                   }
            self.syn_dict_in = {'model': 'inhibitory',
                    'delay':self.delay,
                   }
        else:
            self.d_min = self.delay[0]#3.5
            self.d_max = self.delay[1]#5.5
            self.syn_dict = {'model': 'excitatory',
                    'delay':{'distribution': 'uniform', 'low': self.d_min,
                             'high': self.d_max},
                   }
            print('uniform delay')

            self.syn_dict_in = {'model': 'inhibitory',
                    'delay': {'distribution': 'uniform', 'low': self.d_min,
                              'high': self.d_max},
                   }
        nest.Connect(self.noise,self.nodes_al, syn_spec=self.syn_dict)
        #nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict)
        #nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)

        print("Connecting network")

        print("Excitatory connections")


        self.conn_params_ex = {'rule': 'fixed_outdegree', 'outdegree': self.CE}
        nest.Connect(self.nodes_ex, self.nodes_ex+self.nodes_in,
                     self.conn_params_ex, syn_spec =self.syn_dict)

        print("Inhibitory connections")

        self.conn_params_in = {'rule': 'fixed_outdegree', 'outdegree': self.CI}
        nest.Connect(self.nodes_in, self.nodes_ex+self.nodes_in, self.conn_params_in, syn_spec =self.syn_dict_in)

class disconnected_brunel(meta_brunel):
    """ Only Poisson Drive """
    def connect(self):
        """Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})
        self.delay = self.delay[0]
        print('sinlge delay')
        self.syn_dict = {'model': 'excitatory',
                'delay':self.delay,
               }
        self.syn_dict_in = {'model': 'inhibitory',
                'delay':self.delay,
               }

        nest.Connect(self.noise,self.nodes_al, syn_spec=self.syn_dict)
        #nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)

        print("Connecting network")




class ExcitatoryBrunel(meta_brunel):
    """ Only Excitatory Neurons """
        
    def build(self):
        
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                         "tau_m":self.tauMem,
                         "t_ref":self.t_ref,
                         "E_L": 0.0,
                         "V_reset": self.V_res,
                         "V_m": self.V_m,
                         "V_th":self.theta}
        if self.verbose:
            print(self.neuron_params)

        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        self.nu_th =self.theta / (self.J * self.CE * self.tauMem)# (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex * self.CE
        
        nest.SetDefaults("iaf_psc_delta", self.neuron_params)
        nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("iaf_psc_delta",self.NE)
        self.noise    = nest.Create("poisson_generator")
        self.espikes  = nest.Create("spike_detector")
        #self.ispikes  = nest.Create("spike_detector")
        self.nodes_al = self.nodes_ex# +self.nodes_in
        nest.SetStatus(self.espikes,[{"label": self.simulation,
                                 "withtime": True,
                                 "withgid": True,#%(n_i_),
                                 "to_file": True}])

        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True
        self.built = True
    def connect(self):
        """Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                    }
        else:
            self.d_min = self.delay[0]#3.5
            self.d_max = self.delay[1]#5.5
            self.syn_dict = {'model': 'excitatory',
                    'delay':{'distribution': 'uniform', 'low': self.d_min,
                             'high': self.d_max},
                   }
            print('uniform delay')
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)

        print("Connecting network")

        print("Excitatory connections")


        self.conn_params_ex = {'rule': 'fixed_indegree', 'indegree': self.CE}
        nest.Connect(self.nodes_ex, self.nodes_ex,
                     self.conn_params_ex, syn_spec =self.syn_dict)


    def run(self):
        """run the simulation"""
        print("Simulating")
        # New option to simulate in chunk and check a special condition in between
        # Thus we avoid unwanted condition
        if self.chunk:
            self.chunk_times = np.arange(0,self.simtime+1,self.chunk_size)[1::]
            if len(self.chunk_times) != int(self.simtime/self.chunk_size):
               raise ChunkError('check the simulation length') 
            for ch in self.chunk_times: 
                nest.Simulate(self.chunk_size)
                self.events_ex = nest.GetStatus(self.espikes,"n_events")[0]
                self.rate_ex   = self.events_ex/self.ch*1000.0/self.N_rec
                if self.rate_ex>60.0:
                    print('Rate is too high. Finishing simulation at %s ms'%(self.ch))
                    self.simtime = self.ch
                    break
        else:
            nest.Simulate(self.simtime)

        self.events_ex = nest.GetStatus(self.espikes,"n_events")[0]
        #self.events_in = nest.GetStatus(self.ispikes,"n_events")[0]

        self.rate_ex   = self.events_ex/self.simtime*1000.0/self.N_rec
        #self.rate_in   = self.events_in/self.simtime*1000.0/self.N_rec


        self.num_synapses = (nest.GetDefaults("excitatory")["num_connections"])

        if self.verbose:
            print("Brunel network simulation (Python)")
            print("Number of neurons : {0}".format(self.N_neurons))
            print("Number of synapses: {0}".format(self.num_synapses))
            print("       Exitatory  : {0}".format(int(self.CE *self.N_neurons)
                                                   +self.N_neurons))
            print("Rate   : %.2f Hz" % self.rate_ex)
            #print("Inhibitory rate   : %.2f Hz" % self.rate_in)
            print("Simulation time   : %.2f s" % self.simtime)
            
            
class stim_brunel(meta_brunel):
    dt      = float(0.005)    # the resolution in ms for the clock-based
    def build(self,amp = 0,c_start = 0,init_voltage = False,nu = 0.15):#n_st = 0):
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                         "tau_m":self.tauMem,
                         "t_ref":self.t_ref,
                         "E_L": 0.0,
                         "V_reset": self.V_res,
                         "V_m": self.V_m,
                         "V_th":self.theta}
        if self.verbose:
            print(self.neuron_params)

        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        self.nu_th =self.theta / (self.J * self.CE * self.tauMem)# (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex * self.CE

        nest.SetDefaults("iaf_psc_delta_canon", self.neuron_params)
        nest.SetDefaults("poisson_generator_ps",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("iaf_psc_delta_canon",self.NE)
        self.nodes_in = nest.Create("iaf_psc_delta_canon",self.NI)
        self.noise    = nest.Create("poisson_generator_ps")
        self.espikes  = nest.Create("spike_detector")
        if self.record_vol:
            self.voltmeter = nest.Create("multimeter")
        #self.ispikes  = nest.Create("spike_detector")
        self.nodes_al = self.nodes_ex +self.nodes_in
        nest.SetStatus(self.espikes,[{"label": self.simulation,
                                 "withtime": True,
                                 "precise_times":True,
                                 "precision":50,
                                 "withgid": True,#%(n_i_),
                                 "to_file": True}])
        if self.record_vol:
            nest.SetStatus(self.voltmeter,[{"label": self.simulation+'_voltage',
                                    'record_from': ['V_m'],
                                    'interval':0.05,
                                     "withtime": True,
                                     "withgid": True,#%(n_i_),
                                     "to_file": True}])
            
        if init_voltage:
        #    #_ = self.get_analytical_V(nu)
            self.mu_s = 9.76
            self.sigma_s = 2.701
        #    self.vinit = np.random.normal(self.mu_s,self.sigma_s,[self.NE+self.NI]) #for FR = 0.1
        #    #self.vinit = np.random.normal(10,3.7,[self.NE+self.NI]) #for FR = 0.1
        #    self.stim_n = np.random.randint(1,self.NE, n_st)
        #    print(self.stim_n)
        #    nest.SetStatus(self.nodes_al, "V_m", self.vinit)
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True

        # DC stimulation
        #t_stim = 0.             # perturbation time (time of the extra spike)
        #fade_out = 0.5*self.delay[0]      # fade out time (ms) *0.5
        #c_start =c_start# 2.4
        #c_stop = c_start#3.5
        #self.amp = amp
        #self.current = nest.Create("dc_generator",
        #                           params={'amplitude': self.amp,
        #                                   'start': c_start,
        #                                   'stop': c_stop+t_stim+fade_out})#10e4
        self.built = True
    def connect(self,n_ext=1,Poisson=True):
        """Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})
        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                   }
            self.syn_dict_in = {'model': 'inhibitory',
                    'delay':self.delay,
                   }
        else:
            self.d_min = self.delay[0]#3.5
            self.d_max = self.delay[1]#5.5
            self.syn_dict = {'model': 'excitatory',
                    'delay':{'distribution': 'uniform', 'low': self.d_min,
                             'high': self.d_max},
                   }
            print('uniform delay')

            self.syn_dict_in = {'model': 'inhibitory',
                    'delay': {'distribution': 'uniform', 'low': self.d_min,
                              'high': self.d_max},
                   }
        if Poisson:
            nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict)
            nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)

        #Stimulation
        #self.random_slice_start = np.random.randint(1,self.NE) 
        #self.random_slice_stop =self.random_slice_start +n_ext
        #print(np.arange(self.random_slice_start,self.random_slice_stop))
        #nest.Connect(self.current,self.nodes_al[self.random_slice_start:self.random_slice_stop])#

        #print(self.random_slice_start)
        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_al[:self.N_rec])
        print("Connecting network")

        print("Excitatory connections")


        self.conn_params_ex = {'rule': 'fixed_indegree', 'indegree': self.CE}
        nest.Connect(self.nodes_ex, self.nodes_ex+self.nodes_in,
                     self.conn_params_ex, syn_spec =self.syn_dict)

        print("Inhibitory connections")

        self.conn_params_in = {'rule': 'fixed_indegree', 'indegree': self.CI}
        nest.Connect(self.nodes_in, self.nodes_ex+self.nodes_in,
                     self.conn_params_in, syn_spec =self.syn_dict_in)

    def run(self,repeat,n_st,prime = False):
        """run the simulation"""
        print("Simulating")
        # New option to simulate in chunk and check a special condition in between
        # Thus we avoid unwanted condition
        for rep in range(repeat):
            if prime==True:
                nest.SetStatus(self.nodes_al,"V_m", np.zeros([self.NE+self.NI]))
                nest.Simulate(350.0)
                self.vinit = np.array(nest.GetStatus(self.nodes_al,"V_m"))
                print(np.mean(self.vinit),np.std(self.vinit))
            else:
                self.mu_s =0# 7.5
                self.sigma =0.5# 3.3
                self.vinit = np.random.normal(self.mu_s,self.sigma_s,[self.NE+self.NI]) #for FR = 0.1
            if n_st>0:
                self.stim_n = np.random.randint(1,self.NE, n_st)
                self.vinit[self.stim_n] =20.05# self.theta+0.05#20.05
                print(self.stim_n)
            print(np.sum(self.vinit>20.0),'spikes')
            nest.SetStatus(self.nodes_al, "V_m", self.vinit)
            nest.Simulate(self.simtime)

    def get_analytical_V(self,nu=0.15):
        #Analytical mean and std
	#J = self.j#(0.85*1.5)
	#eta = 0.5
	#Ce = 800*0.1
	#tau = 0.04 # membrane timescale in s
        tau = self.tauMem/1000
        self.nu_thr = self.theta/(self.J*self.CE*tau)
	#nu = 0.15
	#g= 4
        self.gamma = self.CI/self.CE
        mu_int = self.CE*self.J*(1-(self.gamma*self.g))*nu*tau
        mu_ext = (self.CE*self.J*(self.eta*self.nu_thr)*tau)
        self.mu_an = mu_int+mu_ext
        sigma_int = self.J*np.sqrt(self.CE*(1+(self.gamma*(self.g**2)))*nu*tau)
        sigma_ext = self.J*np.sqrt((self.CE*(self.eta*self.nu_thr)*tau))
        self.sigma_an = np.sqrt((sigma_int**2)+(sigma_ext**2))
        print(self.mu_an,self.sigma_an)
        return self.mu_an, self.sigma_an
class PoisBrunel(meta_brunel):
    def build(self):
        
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                         "tau_m":self.tauMem,
                         "t_ref":self.t_ref,
                         "E_L": 0.0,
                         "V_reset": self.V_res,
                         "V_m": self.V_m,
                         "V_th":self.theta}
        if self.verbose:
            print(self.neuron_params)

        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        self.nu_th =self.theta / (self.J * self.CE * self.tauMem)# (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex * self.CE
        nest.SetDefaults("iaf_psc_delta", self.neuron_params)
        nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("iaf_psc_delta",self.NE)
        self.nodes_in = nest.Create("iaf_psc_delta",self.NI)
        self.noise    = nest.Create("poisson_generator")
        self.espikes  = nest.Create("spike_detector")
        if self.record_vol:
            self.voltmeter = nest.Create("multimeter")
        #self.ispikes  = nest.Create("spike_detector")
        self.nodes_al = self.nodes_ex +self.nodes_in
        nest.SetStatus(self.espikes,[{"label": self.simulation,
                                 "withtime": True,
                                 "withgid": True,#%(n_i_),
                                 "to_file": True}])
        if self.record_vol:
            nest.SetStatus(self.voltmeter,[{"label": self.simulation+'_voltage',
                                    'record_from': ['V_m'],
                                    'interval':0.05,
                                     "withtime": True,
                                     "withgid": True,#%(n_i_),
                                     "to_file": True}])
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True
        self.vinit =np.random.normal(10,3.5,[self.NE+self.NI]) #for FR = 0.1
        #self.vinit[:]  = -25.
        nest.SetStatus(self.nodes_al, "V_m", self.vinit)

        self.built = True
    def connect(self):
        """ Only Poisson Drive Connect nodes"""
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})
        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                   }
            self.syn_dict_in = {'model': 'inhibitory',
                    'delay':self.delay,
                   }
        else:
            self.d_min = self.delay[0]#3.5
            self.d_max = self.delay[1]#5.5
            self.syn_dict = {'model': 'excitatory',
                    'delay':{'distribution': 'uniform', 'low': self.d_min,
                             'high': self.d_max},
                   }
            print('uniform delay')

            self.syn_dict_in = {'model': 'inhibitory',
                    'delay': {'distribution': 'uniform', 'low': self.d_min,
                              'high': self.d_max},
                   }
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict)
        nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)

        print("Connecting network")
