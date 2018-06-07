import nest
import nest.raster_plot
import time
from numpy import exp
import numpy as np
import h5py as h5py

class ChunkError(Exception): pass

def save_dict(espikes, ispikes,simulation,conn_meta = None):
        # save current simulation
        ts, gids = nest.raster_plot._from_memory(espikes)
        ts_i, gids_i = nest.raster_plot._from_memory(ispikes)
        if 'npy' in simulation.name:
            dictionary = {'ts':ts,
                          'gids':gids,
                          'ts_i':ts_i,
                          'gids_i':gids_i,
                          'conn':conn_meta}
            np.save(simulation, dictionary) 
            print('saved as npy')
        else:
            with h5py.File(simulation,'w') as f:
                f['ts'] = ts
                f['gids'] = gids
                f['ts_i'] = ts_i
                f['gids_i'] = gids_i
                f['conn'] = conn_meta
            print('saved as hdf5')
            
def get_Noutdeg(N_rec):
    out_deg = np.zeros(N_rec)
    for indx,i in enumerate(np.arange(1,N_rec+1)):
        conn = nest.GetConnections(source=[i])
        out_deg[indx] = len(nest.GetStatus(conn))
    return out_deg

def alpha_brunel(g=7, # inhibitory strenght
                 eta = 3.5, #Poisson rate ratio
                 d=3.5, # synaptic delay
                 J=0.1, #synaptic strength
                 NE =10000, # fraction of inh neurons
                 NI =2500,
                 N_rec = 50,
                 epsilon = 0.1,
                 tauMem= 20.0,
                 CMem = 250.0,
                 tauSyn=0.5,
                 V_m = 0.0,
                 V_res = 10.0,
                 simtime=1000,
                 verbose = True):
    
    
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
    
    
    nest.ResetKernel()
    startbuild = time.time()
    dt      = float(0.05)    # the resolution in ms for the clock-based

    #Definition of the number of neurons in the network and the number of neuron recorded from
    NE =int(NE) #total - NI
    NI  =int(NI)# int(n_i_*total)
    epsilon = 0.1  # connection probability
    N_neurons = NE+NI   # number of neurons in total
    CE    = int(epsilon*NE) # number of excitatory synapses per neuron
    CI    = int(epsilon*NI) # number of inhibitory synapses per neuron  
    C_tot = int(CI+CE)      # total number of synapses per neuron
    if verbose:
        #neurons output
        print('Number of neurons (E/I):')
        print(NE,NI)
        print('Connection Probability: %s'%(epsilon))

    #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
    delay   = d    # synaptic delay in ms
    tauSyn = 0.5  # synaptic time constant in ms
    tauMem = 20.0  # time constant of membrane potential in ms
    CMem = 250.0  # capacitance of membrane in in pF
    theta = 20.0  # membrane threshold potential in mV
    neuron_params = {"C_m": CMem,
                     "tau_m": tauMem,
                     "tau_syn_ex": tauSyn,
                     "tau_syn_in": tauSyn,
                     "t_ref": 2.0,
                     "E_L": 0.0,
                     "V_reset": V_res,
                     "V_m": V_m,
                     "V_th": theta}
    if verbose:
        print(neuron_params)
    #J     = 0.1     # postsynaptic amplitude in mV

    #Synapses   
    J_unit = ComputePSPnorm(tauMem, CMem, tauSyn)
    J_ex = J / J_unit
    #J_ex  = J       # amplitude of excitatory postsynaptic potential
    J_in  = -g*J_ex # amplitude of inhibitory postsynaptic potential
    if verbose:
        print('J_i = %s'%(J_in))

    #Poisson Rate
    nu_th = (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
    nu_ex = eta * nu_th
    p_rate = 1000.0 * nu_ex * CE


    nest.SetKernelStatus({"resolution": dt, "print_time": True, "overwrite_files": True})
    print("Building network")

    nest.SetDefaults("iaf_psc_alpha", neuron_params)
    nest.SetDefaults("poisson_generator",{"rate": p_rate})
    nodes_ex = nest.Create("iaf_psc_alpha",NE)
    nodes_in = nest.Create("iaf_psc_alpha",NI)
    noise    = nest.Create("poisson_generator")
    espikes  = nest.Create("spike_detector")
    ispikes  = nest.Create("spike_detector")

    nest.SetStatus(espikes,[{"label": "brunel-py-ex",
                             "withtime": True,
                             "withgid": True}])#%(n_i_),
                             #"to_file": True

    nest.SetStatus(ispikes,[{"label": "brunel-py-in",
                             "withtime": True,
                             "withgid": True}])#%(n_i_)"to_file": True

    print("Connecting devices")

    nest.CopyModel("static_synapse","excitatory",{"weight":J_ex, "delay":delay})
    nest.CopyModel("static_synapse","inhibitory",{"weight":J_in, "delay":delay})

    nest.Connect(noise,nodes_ex, syn_spec="excitatory")
    nest.Connect(noise,nodes_in, syn_spec="excitatory")

    nest.Connect(nodes_ex[:N_rec], espikes, syn_spec="excitatory")
    nest.Connect(nodes_in[:N_rec], ispikes, syn_spec="excitatory")

    print("Connecting network")

    print("Excitatory connections")


    conn_params_ex = {'rule': 'fixed_indegree', 'indegree': CE}
    nest.Connect(nodes_ex, nodes_ex+nodes_in, conn_params_ex, "excitatory")

    print("Inhibitory connections")

    conn_params_in = {'rule': 'fixed_indegree', 'indegree': CI}
    nest.Connect(nodes_in, nodes_ex+nodes_in, conn_params_in, "inhibitory")

    endbuild=time.time()

    print("Simulating")

    nest.Simulate(simtime)
    endsimulate= time.time()

    events_ex = nest.GetStatus(espikes,"n_events")[0]
    events_in = nest.GetStatus(ispikes,"n_events")[0]

    rate_ex   = events_ex/simtime*1000.0/N_rec
    rate_in   = events_in/simtime*1000.0/N_rec


    num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
                   nest.GetDefaults("inhibitory")["num_connections"])

    build_time = endbuild-startbuild
    sim_time   = endsimulate-endbuild
    print("Brunel network simulation (Python)")
    print("Number of neurons : {0}".format(N_neurons))
    print("Number of synapses: {0}".format(num_synapses))
    print("       Exitatory  : {0}".format(int(CE * N_neurons) + N_neurons))
    print("       Inhibitory : {0}".format(int(CI * N_neurons)))
    print("Excitatory rate   : %.2f Hz" % rate_ex)
    print("Inhibitory rate   : %.2f Hz" % rate_in)
    print("Building time     : %.2f s" % build_time)
    print("Simulation time   : %.2f s" % sim_time)
    
    return ispikes,espikes
    
def delta_brunel(g=7, # inhibitory strenght
                 eta = 3.5, #Poisson rate ratio
                 d=[3.5,5.5], # synaptic delay
                 J=0.1, #synaptic strength
                 NE =10000, # fraction of inh neurons
                 NI =2500,
                 N_rec = 50,
                 epsilon = 0.1,
                 theta = 20.0,
                 t_ref = 2.0, #refractory period
                 tauMem= 20.0,
                 CMem = 1.0,
                 V_m = 0.0,
                 V_res = 10.0,
                 simtime=1000,
                 verbose = True,
                 chunk = False,
                 chunk_size = 1000,
                 deg = False):
    
    nest.ResetKernel()
    startbuild = time.time()
    dt      = float(0.05)    # the resolution in ms for the clock-based

    #Definition of the number of neurons in the network and the number of neuron recorded from
    NE =int(NE) #total - NI
    NI  =int(NI)# int(n_i_*total)
    N_neurons = NE+NI   # number of neurons in total
    CE    = int(epsilon*NE) # number of excitatory synapses per neuron
    CI    = int(epsilon*NI) # number of inhibitory synapses per neuron  
    C_tot = int(CI+CE)      # total number of synapses per neuron
    if verbose:
        #neurons output
        print('Number of neurons (E/I):')
        print(NE,NI)
        print('Connection Probability: %s'%(epsilon))

    #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
    delay   = d    # synaptic delay in ms
    CMem = CMem  # capacitance of membrane in in pF
    neuron_params = {"C_m": CMem,
                     "tau_m": tauMem,
                     "t_ref":t_ref,
                     "E_L": 0.0,
                     "V_reset": V_res,
                     "V_m": V_m,
                     "V_th": theta}
    if verbose:
        print(neuron_params)
    #J     = 0.1     # postsynaptic amplitude in mV

    #Synapses   
    J_ex  = J       # amplitude of excitatory postsynaptic potential
    J_in  = -g*J_ex # amplitude of inhibitory postsynaptic potential
    if verbose:
        print('J_i = %s'%(J_in))

    #Poisson Rate
    nu_th =theta / (J * CE * tauMem)# (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
    nu_ex = eta * nu_th
    p_rate = 1000.0 * nu_ex * CE


    nest.SetKernelStatus({"resolution": dt, "print_time": True, "overwrite_files": True})
    print("Building network")

    nest.SetDefaults("iaf_psc_delta", neuron_params)
    nest.SetDefaults("poisson_generator",{"rate": p_rate})
    nodes_ex = nest.Create("iaf_psc_delta",NE)
    nodes_in = nest.Create("iaf_psc_delta",NI)
    noise    = nest.Create("poisson_generator")
    espikes  = nest.Create("spike_detector")
    ispikes  = nest.Create("spike_detector")

    nest.SetStatus(espikes,[{"label": "brunel-py-ex",
                             "withtime": True,
                             "withgid": True}])#%(n_i_),
                             #"to_file": True

    nest.SetStatus(ispikes,[{"label": "brunel-py-in",
                             "withtime": True,
                             "withgid": True}])#%(n_i_)"to_file": True

    print("Connecting devices")
    nest.CopyModel("static_synapse","excitatory",{"weight":J_ex})
    nest.CopyModel("static_synapse","inhibitory",{"weight":J_in})
    if len(d)<2:
        d = d[0]
        print('sinlge delay')
        syn_dict = {'model': 'excitatory',
                'delay':d,
               }
        syn_dict_in = {'model': 'inhibitory',
                'delay':d,
               }
    else:
        d_min = d[0]#3.5
        d_max = d[1]#5.5
        syn_dict = {'model': 'excitatory',
                'delay':{'distribution': 'uniform', 'low': d_min, 'high': d_max},
               }
        print('uniform delay')

        syn_dict_in = {'model': 'inhibitory',
                'delay': {'distribution': 'uniform', 'low': d_min, 'high': d_max},
               }
    nest.Connect(noise,nodes_ex, syn_spec=syn_dict)
    nest.Connect(noise,nodes_in, syn_spec=syn_dict)

    nest.Connect(nodes_ex[:N_rec], espikes, syn_spec=syn_dict)
    nest.Connect(nodes_in[:N_rec], ispikes, syn_spec=syn_dict)

    print("Connecting network")

    print("Excitatory connections")


    conn_params_ex = {'rule': 'fixed_indegree', 'indegree': CE}
    nest.Connect(nodes_ex, nodes_ex+nodes_in, conn_params_ex, syn_spec =syn_dict)

    print("Inhibitory connections")

    conn_params_in = {'rule': 'fixed_indegree', 'indegree': CI}
    nest.Connect(nodes_in, nodes_ex+nodes_in, conn_params_in, syn_spec =syn_dict_in)

    endbuild=time.time()

    print("Simulating")
    # New option to simulate in chunk and check a special condition in between
    # Thus we avoid unwanted condition
    if chunk:
        chunk_times = np.arange(0,simtime+1,chunk_size)[1::]
        if len(chunk_times) != int(simtime/chunk_size):
           raise ChunkError('check the simulation length') 
        for ch in chunk_times: 
            nest.Simulate(chunk_size)
            events_ex = nest.GetStatus(espikes,"n_events")[0]
            rate_ex   = events_ex/ch*1000.0/N_rec
            if rate_ex>60.0:
                print('Rate is too high. Finishing simulation at %s ms'%(ch))
                simtime = ch
                break
    else:
        nest.Simulate(simtime)
    endsimulate= time.time()

    events_ex = nest.GetStatus(espikes,"n_events")[0]
    events_in = nest.GetStatus(ispikes,"n_events")[0]

    rate_ex   = events_ex/simtime*1000.0/N_rec
    rate_in   = events_in/simtime*1000.0/N_rec


    num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
                   nest.GetDefaults("inhibitory")["num_connections"])

    
    build_time = endbuild-startbuild
    sim_time   = endsimulate-endbuild
    print("Brunel network simulation (Python)")
    print("Number of neurons : {0}".format(N_neurons))
    print("Number of synapses: {0}".format(num_synapses))
    print("       Exitatory  : {0}".format(int(CE * N_neurons) + N_neurons))
    print("       Inhibitory : {0}".format(int(CI * N_neurons)))
    print("Excitatory rate   : %.2f Hz" % rate_ex)
    print("Inhibitory rate   : %.2f Hz" % rate_in)
    print("Building time     : %.2f s" % build_time)
    print("Simulation time   : %.2f s" % sim_time)
    if deg:
        conn_meta = get_Noutdeg(N_rec)
        return ispikes,espikes,conn_meta
    
    return ispikes,espikes





def fixed_brunel(g=7, # inhibitory strenght
                 eta = 3.5, #Poisson rate ratio
                 d=[3.5,5.5], # synaptic delay
                 J=0.1, #synaptic strength
                 NE =10000, # fraction of inh neurons
                 NI =2500,
                 N_rec = 50,
                 epsilon = 0.1,
                 theta = 20.0,
                 t_ref = 2.0, #refractory period
                 tauMem= 20.0,
                 CMem = 1.0,
                 V_m = 0.0,
                 V_res = 10.0,
                 simtime=1000,
                 verbose = True,
                 chunk = False,
                 chunk_size = 1000,
                 deg = False,
                 edges = None):
    
    nest.ResetKernel()
    startbuild = time.time()
    dt      = float(0.05)    # the resolution in ms for the clock-based

    #Definition of the number of neurons in the network and the number of neuron recorded from
    NE =int(NE) #total - NI
    NI  =int(NI)# int(n_i_*total)
    N_neurons = NE+NI   # number of neurons in total
    CE    = int(epsilon*NE) # number of excitatory synapses per neuron
    CI    = int(epsilon*NI) # number of inhibitory synapses per neuron  
    C_tot = int(CI+CE)      # total number of synapses per neuron
    if verbose:
        #neurons output
        print('Number of neurons (E/I):')
        print(NE,NI)
        print('Connection Probability: %s'%(epsilon))

    #Initialization of the parameters of the integrate and fire neuron and the synapses. The parameter of the neuron are stored in a dictionary.
    delay   = d    # synaptic delay in ms
    CMem = CMem  # capacitance of membrane in in pF
    neuron_params = {"C_m": CMem,
                     "tau_m": tauMem,
                     "t_ref":t_ref,
                     "E_L": 0.0,
                     "V_reset": V_res,
                     "V_m": V_m,
                     "V_th": theta}
    if verbose:
        print(neuron_params)
    #J     = 0.1     # postsynaptic amplitude in mV

    #Synapses   
    J_ex  = J       # amplitude of excitatory postsynaptic potential
    J_in  = -g*J_ex # amplitude of inhibitory postsynaptic potential
    if verbose:
        print('J_i = %s'%(J_in))

    #Poisson Rate
    nu_th =theta / (J * CE * tauMem)# (theta * CMem) / (J_ex * CE * exp(1) * tauMem * tauSyn)
    nu_ex = eta * nu_th
    p_rate = 1000.0 * nu_ex * CE


    nest.SetKernelStatus({"resolution": dt, "print_time": True, "overwrite_files": True})
    print("Building network")

    nest.SetDefaults("iaf_psc_delta", neuron_params)
    nest.SetDefaults("poisson_generator",{"rate": p_rate})
    nodes_ex = nest.Create("iaf_psc_delta",NE)
    nodes_in = nest.Create("iaf_psc_delta",NI)
    noise    = nest.Create("poisson_generator")
    espikes  = nest.Create("spike_detector")
    ispikes  = nest.Create("spike_detector")

    nest.SetStatus(espikes,[{"label": "brunel-py-ex",
                             "withtime": True,
                             "withgid": True}])#%(n_i_),
                             #"to_file": True

    nest.SetStatus(ispikes,[{"label": "brunel-py-in",
                             "withtime": True,
                             "withgid": True}])#%(n_i_)"to_file": True

    print("Connecting devices")
    nest.CopyModel("static_synapse","excitatory",{"weight":J_ex})
    nest.CopyModel("static_synapse","inhibitory",{"weight":J_in})
    if len(d)<2:
        d = d[0]
        print('sinlge delay')
        syn_dict = {'model': 'excitatory',
                'delay':d,
               }
        syn_dict_in = {'model': 'inhibitory',
                'delay':d,
               }
    else:
        d_min = d[0]#3.5
        d_max = d[1]#5.5
        syn_dict = {'model': 'excitatory',
                'delay':{'distribution': 'uniform', 'low': d_min, 'high': d_max},
               }
        print('uniform delay')

        syn_dict_in = {'model': 'inhibitory',
                'delay': {'distribution': 'uniform', 'low': d_min, 'high': d_max},
               }
    nest.Connect(noise,nodes_ex, syn_spec=syn_dict)
    nest.Connect(noise,nodes_in, syn_spec=syn_dict)

    nest.Connect(nodes_ex[:N_rec], espikes, syn_spec=syn_dict)
    nest.Connect(nodes_in[:N_rec], ispikes, syn_spec=syn_dict)

    print("Connecting network")

    print("Excitatory connections")
 
    #conn_params_ex = {'rule': 'fixed_indegree', 'indegree': CE}
    #nest.Connect(nodes_ex, nodes_ex+nodes_in, conn_params_ex, syn_spec =syn_dict)
    all_nodes = nodes_ex+nodes_in
    #print(CE)
    #print(nodes_in)
    for edge in edges[:800000]:
        nest.Connect([edge[0]], [edge[1]], syn_spec =syn_dict)

    print("Inhibitory connections")

    #conn_params_in = {'rule': 'fixed_indegree', 'indegree': CI}
    #nest.Connect(nodes_in, nodes_ex+nodes_in, conn_params_in, syn_spec =syn_dict_in)
    for edge in edges[800000:]:
        nest.Connect([edge[0]], [edge[1]], syn_spec =syn_dict_in)
    endbuild=time.time()

    print("Simulating")
    # New option to simulate in chunk and check a special condition in between
    # Thus we avoid unwanted condition
    if chunk:
        chunk_times = np.arange(0,simtime+1,chunk_size)[1::]
        if len(chunk_times) != int(simtime/chunk_size):
           raise ChunkError('check the simulation length') 
        for ch in chunk_times: 
            nest.Simulate(chunk_size)
            events_ex = nest.GetStatus(espikes,"n_events")[0]
            rate_ex   = events_ex/ch*1000.0/N_rec
            if rate_ex>60.0:
                print('Rate is too high. Finishing simulation at %s ms'%(ch))
                simtime = ch
                break
    else:
        nest.Simulate(simtime)
    endsimulate= time.time()

    events_ex = nest.GetStatus(espikes,"n_events")[0]
    events_in = nest.GetStatus(ispikes,"n_events")[0]

    rate_ex   = events_ex/simtime*1000.0/N_rec
    rate_in   = events_in/simtime*1000.0/N_rec


    num_synapses = (nest.GetDefaults("excitatory")["num_connections"] +
                   nest.GetDefaults("inhibitory")["num_connections"])

    
    build_time = endbuild-startbuild
    sim_time   = endsimulate-endbuild
    print("Brunel network simulation (Python)")
    print("Number of neurons : {0}".format(N_neurons))
    print("Number of synapses: {0}".format(num_synapses))
    print("       Exitatory  : {0}".format(int(CE * N_neurons) + N_neurons))
    print("       Inhibitory : {0}".format(int(CI * N_neurons)))
    print("Excitatory rate   : %.2f Hz" % rate_ex)
    print("Inhibitory rate   : %.2f Hz" % rate_in)
    print("Building time     : %.2f s" % build_time)
    print("Simulation time   : %.2f s" % sim_time)
    if deg:
        conn_meta = get_Noutdeg(N_rec)
        return ispikes,espikes,conn_meta
    
    return ispikes,espikes

def get_Noutdeg(N_rec):
    #out_deg = np.zeros(N_rec)
    out_deg = []#
    for indx,i in enumerate(np.arange(1,N_rec+1)):
        conn = nest.GetConnections(source=[i])
        #out_deg[indx] = len(nest.GetStatus(conn))
        out_deg.append(nest.GetStatus(conn))
    return out_deg
