from func.brunel_meta import *
from func.helpers import *
#
class tauMem_error(Exception): pass
class AdExModel(meta_brunel):

    threads = 4
    built = False
    connected = False
    dt      = float(0.001)    # the resolution in ms for the clock-based
    CMem = 400.0
    V_m = 0.0
    V_res = 10.0#10
    t_ref = 2.0
    theta = 20.0
    def build(self,
              Cmem = 4.,
              g_L = 0.1,
              Delta_T = 0.,
              tau_w=1000.,
              a=0.,
              b=100., constant_k = False):
        self.CMem = Cmem

        if self.tauMem != self.CMem/g_L:
            raise tauMem_error('C/g_L is not equal specified tau_mem. \
                               AdEx uses conductance(g_L) to calculate tau_mem')
        if constant_k:
            self.CE = int(self.NE* constant_k[0])
            self.CI = int(self.NI*constant_k[1])
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.neuron_params = {"C_m":self.CMem,
                              "V_peak":20.,
                              "I_e":0.,
                              "E_L": 0.0,
                              "g_L":g_L,
                              "t_ref":self.t_ref,
                              'Delta_T':Delta_T,
                              "V_reset": self.V_res,
                              "V_m": self.V_m,
                              "V_th":self.theta,
                              "tau_w" : tau_w,
                              'gsl_error_tol':1e-6,
                              "a":a,
                              "b":b}
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
        
        nest.SetDefaults("aeif_psc_delta", self.neuron_params)
        
        nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("aeif_psc_delta",self.NE)
        self.nodes_in = nest.Create("aeif_psc_delta",self.NI)
        self.noise    = nest.Create("poisson_generator")
        self.espikes  = nest.Create("spike_detector")
        if self.record_vol:
            self.voltmeter = nest.Create("multimeter")
        #self.ispikes  = nest.Create("spike_detector")
        self.nodes_al = self.nodes_ex +self.nodes_in
        nest.SetStatus(self.nodes_ex,[{ 'V_m':self.V_res}])
        nest.SetStatus(self.espikes,[{"label": self.simulation,
                                 "withtime": True,
                                 "to_memory": True,
                                 "withgid": True,#%(n_i_),
                                 "to_file": True}])
        if self.record_vol:
            nest.SetStatus(self.voltmeter,[{"label": self.simulation+'_voltage',
                                    'record_from': ['V_m'],
                                    'interval':1.,
                                     "withtime": True,
                                     "to_memory":False,
                                     "withgid": True,#%(n_i_),
                                     "to_file": True}])
        #nest.SetStatus(self.ispikes,[{"label": "brunel-py-in",
        #                         "withtime": True,
        #                         "withgid": True}])#%(n_i_)"to_file": True
        self.built = True
        
        
    def connect(self,j_cnqx=False):
        print(nest.GetStatus([1],'C_m'))

        if self.built ==False:

            raise BuildError('Build the network first')

        if j_cnqx !=False:
            self.J_ext = self.J_ex
            self.J_ex = j_cnqx
            print('APPLYING CNQX')
        else:
            self.J_ext = self.J_ex
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
#        nest.Connect(self.noise,self.nodes_in, syn_spec=self.syn_dict)

        nest.Connect(self.nodes_al[:self.N_rec],self.espikes, syn_spec=self.syn_dict)
        #nest.Connect(self.nodes_in[:self.N_rec] ,self.ispikes, syn_spec=self.syn_dict)
        if self.record_vol:
            nest.Connect(self.voltmeter, self.nodes_al[:self.N_rec])
        print("Connecting network")

        print("Excitatory connections")


        self.conn_params_ex = {'rule': 'fixed_indegree', 'indegree': self.CE}
       # nest.Connect(self.nodes_ex, self.nodes_ex+self.nodes_in,
       #              self.conn_params_ex, syn_spec =self.syn_dict)

        print("Inhibitory connections")

        self.conn_params_in = {'rule': 'fixed_indegree', 'indegree': self.CI}
        #nest.Connect(self.nodes_in, self.nodes_ex+self.nodes_in,
        #             self.conn_params_in, syn_spec =self.syn_dict_in)
        
    
class ExAdModel(AdExModel):
    """ Model with asynchronous adaptation """
    threads = 4
    built = False
    connected = False
    dt      = float(1.)    # the resolution in ms for the clock-based
    CMem = 200.0
    V_m = 0.0
    V_res = 10.0#10
    t_ref = 2.0
    theta = 20.0
    def build(self,
              Cmem = 400.,
              g_L = 10.,
              Delta_T = 0.,
              tau_w=1000.,
              a=20.,
              b=100., constant_k = False):
        self.CMem = Cmem

        if self.tauMem != self.CMem/g_L:
            raise tauMem_error('C/g_L is not equal specified tau_mem. \
                               AdEx uses conductance(g_L) to calculate tau_mem')
        if constant_k:
            self.CE = int(self.NE* constant_k[0])
            self.CI = int(self.NI*constant_k[1])
        #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
        self.E_neuron_params = {"C_m":self.CMem,
                              "V_peak":20.001,
                              "I_e":0.,
                              "E_L": 0.0,
                              "g_L":g_L,
                              "t_ref":self.t_ref,
                              'Delta_T':Delta_T,
                              "V_reset": self.V_res,
                              "V_m": self.V_m,
                              "V_th":self.theta,
                              "tau_w" : tau_w,
                              'gsl_error_tol':1e-6,
                              "a":a,
                              "b":b}
#                         "tau_syn_in":tauSyn }

        self.I_neuron_params = {"C_m":self.CMem,
                         "tau_m":self.tauMem,
                         "t_ref":self.t_ref,
                         "E_L": 0.0,
                         "V_reset": self.V_res,
                         "V_m": self.V_m,
                         "V_th":self.theta}
        if self.verbose:
            print(self.E_neuron_params)

        self.J_ex  =self.J       # amplitude of excitatory postsynaptic potential
        self.J_in  = -self.g*self.J_ex # amplitude of inhibitory postsynaptic potential
        if self.verbose:
            print('J_i = %s'%(self.J_in))

        #Poisson Rate
        #self.nu_th =(self.theta * self.CMem) / (self.J_ex * self.CE * np.exp(1) * self.tauMem * tauSyn)
        self.nu_th= self.theta / (self.J * self.CE * self.tauMem)## (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
        self.nu_ex = self.eta * self.nu_th
        self.p_rate = 1000.0 * self.nu_ex  * self.CE
        
        nest.SetDefaults("aeif_psc_delta", self.E_neuron_params)
        nest.SetDefaults("iaf_psc_delta", self.I_neuron_params)
        
        nest.SetDefaults("poisson_generator",{"rate":self.p_rate})
        self.nodes_ex = nest.Create("aeif_psc_delta",self.NE)
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
        
        
    def connect(self,j_cnqx=False):
        print(nest.GetStatus([1],'C_m'))

        if self.built ==False:

            raise BuildError('Build the network first')

        if j_cnqx !=False:
            self.J_ext = self.J_ex
            self.J_ex = j_cnqx
            print('APPLYING CNQX')
        else:
            self.J_ext = self.J_ex
        nest.CopyModel("static_synapse","excitatory",{"weight":self.J_ex})
        nest.CopyModel("static_synapse","inhibitory",{"weight":self.J_in})
        #nest.CopyModel("static_synapse","I_inhibitory",{"weight":20.})
       # nest.CopyModel("static_synapse","I_excitatory",{"weight":20.})
        if len(self.delay)<2:
            self.delay = self.delay[0]
            print('sinlge delay')
            self.syn_dict = {'model': 'excitatory',
                    'delay':self.delay,
                   }
            self.syn_dict_in = {'model': 'inhibitory',
                    'delay':self.delay,
                   }
            self.I_syn_dict = {'model': 'I_excitatory',
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
        
