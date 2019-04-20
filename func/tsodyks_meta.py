from func.brunel_meta import *
from func.helpers import *

class DeprSynapses(meta_brunel):

    threads = 4
    built = False
    connected = False
    dt      = float(1.)    # the resolution in ms for the clock-based
    CMem = 1.0
    V_m = 10.0
    V_res = 10.0#10
    t_ref = 2.0
    theta = 20.0
    def build(self,U,tauSyn=3,constant_k = False):
        self.U = U
        if constant_k:
            self.CE = int(self.NE* constant_k[0])
            self.CI = int(self.NI*constant_k[1])
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                         "tau_m":self.tauMem,
#                         "t_ref_abs":self.t_ref,
#                         "t_ref_tot":self.t_ref,# set abs and total equall
                         #"t_ref":self.t_ref,
                         "E_L": 0.0,
                         "V_reset": self.V_res,
                         "V_m": self.V_m,
                         "V_th":self.theta}
#                         'tau_syn_ex': tauSyn,
#                         "tau_syn_in":tauSyn }
        if self.verbose:
            print(self.neuron_params)

        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        #self.nu_th =(self.theta * self.CMem) / (self.J_ex * self.CE * np.exp(1) * self.tauMem * tauSyn)
        self.nu_th= self.theta / (self.J * self.CE * self.tauMem)## (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex  * self.CE
        
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
                                    'interval':1.,
                                     "withtime": True,
                                     "withgid": True,#%(n_i_),
                                     "to_file": True}])
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True
        self.built = True
        
        
    def connect(self,j_cnqx=False, U=1., x=1., tau_rec=1000. , tau_fac= 0.0,tau_psc= 1.,u=1.):
        
#        self.J_ext = self.J
        self.dep_params = {"U":0.035,
                           # Fasccilitation increment!
                           # Carefull between tsodyks1 and tsodyks2 models
                           "u": 1.0,
                           'x': x,
                           "tau_rec": tau_rec,
                           "tau_fac":tau_fac,
                            #"tau_psc":tau_psc,
                           "weight": 1.}

        nest.SetDefaults("tsodyks2_synapse", self.dep_params)
        
        
        if self.built ==False:
            raise BuildError('Build the network first')
        if j_cnqx !=False:
            self.J_ext = self.J_ex
            self.J_ex = j_cnqx
            print('APPLYING CNQX')
        else:
            self.J_ext = self.J_ex
        nest.CopyModel("tsodyks2_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","external",{"weight":self.J_ext})
        nest.CopyModel("tsodyks2_synapse","inhibitory",{"weight":self.J_in})
        self.delay = self.delay[0]
        print('sinlge delay')
        self.syn_dict = {'model': 'excitatory',
                'delay':self.delay,
               }
        self.syn_dict_in = {'model': 'inhibitory',
                'delay':self.delay,
               }
        self.syn_dict_ext = {'model': 'external',
                'delay':self.delay,
               }
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict_ext)
#        nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict_ext)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.noise, self.rec_noise)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_ex[211:212])
            #nest.Connect(self.voltmeter, self.nodes_in[1:2])
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
        """run the simulation"""
        print("Simulating")
        # New option to simulate in chunk and check a special condition in between
        # Thus we avoid unwanted condition
        if stim==1:
            self.nu_th= self.theta / (self.J * self.CE * self.tauMem)## (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
            self.new_nu_ex = stim_eta * self.nu_th
            self.new_p_rate = 1000.0 * self.new_nu_ex  * self.CE
            
            
            self.add_noise    = nest.Create("poisson_generator")
            nest.Simulate(self.simtime)
            nest.SetStatus(self.add_noise,{"rate":self.new_p_rate})
            nest.Connect(self.add_noise,self.nodes_ex[:800], syn_spec=self.syn_dict_ext)
            nest.Connect(self.add_noise,self.nodes_in[:200], syn_spec=self.syn_dict_ext)
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

class ExNetDeprS(DeprSynapses):
    """ Network of excitatory neurons with adaptation """
    def connect(self, U=1., x=1., tau_rec=1000. , tau_fac= 0.0,tau_psc= 1.,u=1.):
        self.J_ext = 20.
        self.dep_params = {"U": U,
                           "u": 1.,
                           'x': x,
                           "tau_rec": tau_rec,
                           "tau_fac":tau_fac,
                            #"tau_psc":tau_psc,
                           "weight": 1.}
        nest.SetDefaults("tsodyks2_synapse", self.dep_params)
        if self.built ==False:
            raise BuildError('Build the network first')
        nest.CopyModel("tsodyks2_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","external",{"weight":self.J_ext})
        #nest.CopyModel("tsodyks2_synapse","inhibitory",{"weight":self.J_in})
        self.delay = self.delay[0]
        print('sinlge delay')
        self.syn_dict = {'model': 'excitatory',
                'delay':self.delay,
               }
        self.syn_dict_ext = {'model': 'external',
                'delay':self.delay,
               }
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict_ext)
        #nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict_ext)
        nest.Connect(self.nodes_ex[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_ex[:self.N_rec])
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

class CurrentDepr(DeprSynapses):
    theta = 20.0
    def build(self, init_E_L, mu,std,vinit =[0],slow_pm=False, taupm=100 ):
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
        #self.p_rate = 1000.0 * self.nu_ex * self.CE

        nest.SetDefaults("iaf_psc_delta", self.neuron_params)
        #nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("iaf_psc_delta",self.NE)
        self.nodes_in = nest.Create("iaf_psc_delta",self.NI)
        #self.noise    = nest.Create("poisson_generator")
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
            
        if init_E_L:
            #_ = self.get_analytical_V(nu)
            self.mu_s = mu#20.1
            self.sigma_s = std#10.7
            if len(vinit)>1:
                print('user specified')
                self.vinit =vinit
            else:
                self.vinit = np.random.uniform(self.mu_s,self.sigma_s,[self.NE+self.NI]) #for FR = 0.1
            self.Npacemakers = np.sum(self.vinit>20)
            print(self.Npacemakers)
            self.pacemakers = np.where(self.vinit>20)[0]
            if slow_pm:
                self.tauM = np.zeros([self.NE+self.NI])
                self.tauM[:] = self.tauMem
                self.tauM[self.pacemakers] = taupm
                nest.SetStatus(self.nodes_al, "tau_m", self.tauM)
            #self.vinit = np.random.normal(10,3.7,[self.NE+self.NI]) #for FR = 0.1
            nest.SetStatus(self.nodes_al, "E_L", self.vinit)
            
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True


        self.built = True

    def connect(self, U=1., x=1., tau_rec=1000. , tau_fac= 0.0,tau_psc= 1.,u=1.):
        
        self.dep_params = {"U": U,
                           "u": 1.,
                           'x': x,
                           "tau_rec": tau_rec,
                           "tau_fac":tau_fac,
                            #"tau_psc":tau_psc,
                           "weight": 1.}

        nest.SetDefaults("tsodyks2_synapse", self.dep_params)
        
        
        if self.built ==False:
            raise BuildError('Build the network first')

        nest.CopyModel("tsodyks2_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","external",{"weight":self.J_ex})
        nest.CopyModel("tsodyks2_synapse","inhibitory",{"weight":self.J_in})
        self.delay = self.delay[0]
        print('sinlge delay')
        self.syn_dict = {'model': 'excitatory',
                'delay':self.delay,
               }
        self.syn_dict_in = {'model': 'inhibitory',
                'delay':self.delay,
               }
        self.syn_dict_ext = {'model': 'external',
                'delay':self.delay,
               }

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_al[:self.N_rec])
        print("Connecting network")

        print("Excitatory connections")


        self.conn_params_ex = {'rule': 'fixed_indegree', 'indegree': self.CE}
       # nest.Connect(self.nodes_ex, self.nodes_ex+self.nodes_in,
        #             self.conn_params_ex, syn_spec =self.syn_dict)

        print("Inhibitory connections")

        self.conn_params_in = {'rule': 'fixed_indegree', 'indegree': self.CI}
       # nest.Connect(self.nodes_in, self.nodes_ex+self.nodes_in,
       #              self.conn_params_in, syn_spec =self.syn_dict_in)

def LambertWm1(x):
    nest.sli_push(x)
    nest.sli_run('LambertWm1')
    y = nest.sli_pop()
    return y


def ComputePSPnorm(tauMem, CMem, tauSyn):
    a = (tauMem / tauSyn)
    b = (1.0 / tauSyn - 1.0 / tauMem)

    # time of maximum
    t_max = 1.0 / b * (-LambertWm1(-exp(-1.0 / a) / a) - 1.0 / a)

    # maximum of PSP for current of unit amplitude
    return (exp(1.0) / (tauSyn * CMem * b) *
            ((exp(-t_max / tauMem) - exp(-t_max / tauSyn)) / b -
             t_max * exp(-t_max / tauSyn)))



class DeprSyn_tauPSC(DeprSynapses ):

    threads = 4
    built = False
    connected = False
    dt      = float(0.1)    # the resolution in ms for the clock-based
    CMem = 1.0
    V_m = 10.0
    V_res = 10.0#10
    t_ref = 2.0
    theta = 20.0
    def build(self,U,tauSyn=3,constant_k = False):
        self.U = U
        self.tauSyn = tauSyn
        if constant_k:
            self.CE = int(self.NE* constant_k[0])
            self.CI = int(self.NI*constant_k[1])
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                              "tau_m":self.tauMem,
#                         "t_ref_abs":self.t_ref,
#                         "t_ref_tot":self.t_ref,# set abs and total equall
                         #"t_ref":self.t_ref,
                              "E_L": 0.0,
                              "V_reset": self.V_res,
                              "V_m": self.V_m,
                              "V_th":self.theta,
                              'tau_syn_ex': tauSyn,
							  'tau_syn_in':tauSyn }
        if self.verbose:
            print(self.neuron_params)
        J_unit = ComputePSPnorm(self.tauMem, self.CMem, self.tauSyn)
        J_ex = self.J / J_unit  # amplitude of excitatory postsynaptic current
        J_in = -self.g * J_ex    # amplitude of inhibitory postsynaptic current
        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        self.nu_th =(self.theta * self.CMem) / (self.J_ex * self.CE * np.exp(1) * self.tauMem * tauSyn)
        #self.nu_th= self.theta / (self.J * self.CE * self.tauMem)## (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex  * self.CE
        
        nest.SetDefaults("iaf_psc_alpha", self.neuron_params)
        
        nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("iaf_psc_alpha",self.NE)
        self.nodes_in = nest.Create("iaf_psc_alpha",self.NI)
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
                                    'interval':1.,
                                     "withtime": True,
                                     "withgid": True,#%(n_i_),
                                     "to_file": True}])
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True
        self.built = True
        
        
    def connect(self,j_cnqx=False, U=1., x=1., tau_rec=1000. , tau_fac= 0.0,tau_psc= 1.,u=1.):
        
#        self.J_ext = self.J
        self.dep_params = {"U":0.035,
                           # Fasccilitation increment!
                           # Carefull between tsodyks1 and tsodyks2 models
                           "u": 1.0,
                           'x': x,
                           "tau_rec": tau_rec,
                           "tau_fac":tau_fac,
                            "tau_psc":tau_psc,
                           "weight": 1.}

        nest.SetDefaults("tsodyks_synapse", self.dep_params)
        
        
        if self.built ==False:
            raise BuildError('Build the network first')
        if j_cnqx !=False:
            self.J_ext = self.J_ex
            self.J_ex = j_cnqx
            print('APPLYING CNQX')
        else:
            self.J_ext = self.J_ex
        nest.CopyModel("tsodyks_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","external",{"weight":self.J_ext})
        nest.CopyModel("tsodyks_synapse","inhibitory",{"weight":self.J_in})
        self.delay = self.delay[0]
        print('sinlge delay')
        self.syn_dict = {'model': 'excitatory',
                'delay':self.delay,
               }
        self.syn_dict_in = {'model': 'inhibitory',
                'delay':self.delay,
               }
        self.syn_dict_ext = {'model': 'external',
                'delay':self.delay,
               }
        nest.Connect(self.noise,self.nodes_ex, syn_spec=self.syn_dict_ext)
#        nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict_ext)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.noise, self.rec_noise)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_ex[211:212])
            #nest.Connect(self.voltmeter, self.nodes_in[1:2])
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
        """run the simulation"""
        print("Simulating")
        # New option to simulate in chunk and check a special condition in between
        # Thus we avoid unwanted condition
        if stim==1:
            self.nu_th= self.theta / (self.J * self.CE * self.tauMem)## (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
            self.new_nu_ex = stim_eta * self.nu_th
            self.new_p_rate = 1000.0 * self.new_nu_ex  * self.CE
            
            
            self.add_noise    = nest.Create("poisson_generator")
            nest.Simulate(self.simtime)
            nest.SetStatus(self.add_noise,{"rate":self.new_p_rate})
            nest.Connect(self.add_noise,self.nodes_ex[:800], syn_spec=self.syn_dict_ext)
            nest.Connect(self.add_noise,self.nodes_in[:200], syn_spec=self.syn_dict_ext)
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
