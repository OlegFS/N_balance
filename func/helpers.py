import nest
import nest.raster_plot
import time
from numpy import exp
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scipy.stats as sps
from itertools import product
from pathlib import Path
import networkx as nx
import pandas
from os import listdir
from os.path import isfile, join
import os
import h5py as h5py
import pandas
from scipy.interpolate import splrep, sproot, splev
na = np.array

def unfoldTimes(ts,gids,n_units):
    ts_unit= []
    ts_unit = [ts[gids==unit] for unit in np.arange(1,n_units) if np.any(gids==unit)]
    return ts_unit

#ts_unit = unfoldTimes(ts1,gids1,8000)

def get_stat(simulation,n_units,sim_time,verbose=False):
    t1,t2 = (0,sim_time)
    ts1,gids1 =load_ts(simulation,t1,t2)
    sim_time = np.max(ts1)
    N = n_units
    rate = np.sum(ts1.astype('bool'))/((sim_time/1000)*N)
    #ts_unit = unfoldTimes(ts1,gids1,n_units)
    #d_ts_unit = [np.diff(i) for i in ts_unit  if len(i)>2]
    #var_unit = [np.var(unit) for unit in d_ts_unit]
#   ISI_var_unit = 
    #mean_unit = [np.mean(unit) for unit in d_ts_unit]
    #ISI_fano = np.mean(var_unit)/np.mean(mean_unit)   
#     Var = np.var(np.diff(ts1))
    
    #Fano = np.var(np.diff(ts1))/np.mean(np.diff(ts1))

    # Spike_count for al neurons
#     sc,_=np.histogram(ts1,np.arange(0,sim_time+1,1))
#     sc = sc/N

#     #Spike_count per neuron
#     sc_unit = []
#     calc = 0
#     sc_unit = [np.histogram(ts,np.arange(0,sim_time,1))[0] for ts in ts_unit if np.any(ts)]
#     fano_unit = [(np.var(unit)/np.mean(unit)) for unit in sc_unit]
#     var_unit = [np.var(unit) for unit in sc_unit]
    #for unit in sc_unit:
    #    fano_unit.append(np.var(unit)/np.mean(unit))
    #    var_unit.append(np.var(unit))

#     chi= np.sqrt(np.var(sc)/np.mean(var_unit))
#     if verbose:
#         print('Rate: %s Hz'%(rate))
#         #print('Var: %s Hz'%(rate))
#     #     print('Mean SC :%s'%(np.mean(sc)))
#     #     print('Variance of SC: %s '%(np.var(sc)))
#     #     print('Fano mean: %s '%(np.mean(fano_unit)))

#     #     print('Chi (synchrony measure): %s '%(chi))
#         print('ISI Variance: %s'%(np.mean(var_unit)))
#         print('ISI Fano:%s'%(ISI_fano))
    return rate#,np.mean(var_unit)
      
    
    
######### General Read and write ###########

class SizeError(Exception): pass 
class Loading(Exception): pass

def read_gdf(mypath,simulation,t):
    """Read spikes from native NEST recordings
    
    Args:
            mypath (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            t (tuple: int,int): time of the simulation to load (ms)
            
    Returns:
            arr: spikes timestamp 
            arr: corresponding unit ids
    
    """
    t1,t2 = t
    # list all files (for multithreading)
    files = [mypath +i for i in listdir(mypath) if simulation in i and 'gdf' in i]#
    #print('Simulations: %s'% simulation)
    #print('Reading from %s files'%len(files))
    if len(files)>4:
        raise SizeError()
        
    data = from_file_pandas(files,t2)
    ts,gids = data[:,1],data[:,0]
    
    ts1 = ts[(ts>t1)*(ts<t2)]
    gids1 = gids[(ts>t1)*(ts<t2)]
    
    return ts1,gids1
    
def from_file_pandas(fname,t2, **kwargs):
    """Use pandas to read gdf file.
    This function is copied from NEST"""
    data = None
    for f in fname:
        try:
            dataFrame = pandas.read_csv(
                f, sep='\s+', lineterminator='\n',
                header=None, index_col=None,
                skipinitialspace=True) #,nrows = t2
            newdata = dataFrame.values
        except:
            continue

        if data is None:
            data = newdata
        else:
            data = np.concatenate((data, newdata))
    return data

def load_ts(simulation,t,inhibitory = False):
    """ Loads either npy or hdf5 simulation
    Example: 
    t1,t2 = (0,200000)
    ts1,gids1 =load_ts(simulation,t1,t2)
    """
    t1,t2 = t
    if inhibitory:
        stamp = 'ts_i'
        uid = 'gids_i'
    else:
        stamp = 'ts'
        uid = 'gids'

    if 'npy' in simulation.name:
        try:
            espikes = np.load(simulation).item()
        except EOFError:
            espikes = {stamp:0,uid:0}
        ts = espikes[stamp]
        gids = espikes[uid]
    else:
        try:
            with h5py.File(simulation,'r') as f:
                ts = np.array(f[stamp])
                gids = np.array(f[uid])
        except EOFError:
            espikes = {stamp:0,uid:0}
    ts1 = ts[(ts>t1)*(ts<t2)]
    gids1 = gids[(ts>t1)*(ts<t2)]
    return ts1,gids1
def stack_ts(simulation,t):
    """Stack inhibitory and excitatory timestamps
    helper
    
    Args:
            simulation (str): Directory of the simulation with "/" and
                              simulation name
            t (tuple: int,int): time of the simulation to load (ms)
            
    Returns:
            arr: concatinated spikes timestamp 
            arr: concatinated unit ids
    
    """
    ts,gids = load_ts(simulation,t,inhibitory = False)
    ts_i,gids_i = load_ts(simulation,t,inhibitory = True)
    return np.hstack((ts,ts_i)),np.hstack((gids,gids_i))
    
    

def load_sim(path,simulation,t):
    """Load spike times and unit ids 
    hdf/gids/npy compatibility
    
    Comment:
    I used hdf5 and npy at earlier stages to save simulation.
    
    Args:   
            path (str): Directory of the simulation with "/".
            simulation (str): simulation name
            t (tuple: int,int): time of the simulation to load (ms)
            
    Returns:
            arr: concatinated spikes timestamp 
            arr: concatinated unit ids
    
    """
    l = 0 # ugly solution to check the right version
    try:
        ts,gids =stack_ts(Path(path+simulation+'.hdf5'),t)
        l = 1
    except:
        pass
    
    try:
        ts,gids =stack_ts(Path(path+simulation+'.npy'),t)
        l = 1
    except:
        pass
    
    if l == 0:
        try:
            ts,gids = read_gdf(path,simulation,t)
        except:    
            raise Loading('No dataset found: check the path and file name')
            
    return ts,gids

def read_voltage(mypath,
                 simulation,
                 N,
                 t,
                 nice_format = False):
    
    """Return a matrix of [unit X time, time , voltage]
    remark: Highly inefficient
    Args:
            mypath (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            N (int): Number of neurons
            t (tuple: int,int): time of the simulation to load (ms)
            TODO: fix the time (now it reads the whole simulation)
            nice_format (bool): reshape into [unit X time]
            
    Returns:
            array: matrix of [unit X time, time , voltage] or
                   [unit X time] 
                   
    """
    # list all files (for multithreading)
    t1,t2 = t
    print(mypath)
    files = [mypath +i for i in listdir(mypath) if simulation in i and 'dat' in i]#
    print('Simulations: %s'% simulation)
    print('Reading from %s files'%len(files))
    data = np.zeros([0,3])
    if len(files)>4:
        raise SizeError()
    for file in files:   
        print(file)
        data = np.concatenate((data,np.loadtxt(file)),axis=0)
    
    # reformat the data
    if nice_format:
        unit = 2
        length = data[data[:,0]==unit][:,2].shape[0]
        voltage = np.zeros([N,length]) 
        for i,unit in enumerate(np.arange(1,N+1)):
            voltage[i,:] = data[data[:,0]==unit][:,2]
        return voltage
    else:
        return data


####### Basic Stats Helpers ######### 
def lazy_rate(path,
             simulation,
             N = 1000,
             sim_time = 200
             ):
    """Computes an average firing rate of the NEST simulation
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            N (int): Number of neurons
            sim_time (int): length of the simulation in s.
            
    Returns:
            float64: N of spikes per neuron per s
    """
    ts,_ = load_sim(path,simulation,(0,sim_time*1000))
    return len(ts)/N/sim_time

def lazy_unit_rate(path,
             simulation,
             N = 1000,
             sim_time = 200
             ):
    """Computes an average firing rate of the NEST simulation
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            N (int): Number of neurons
            sim_time (int): length of the simulation in s.
            
    Returns:
            float64: N of spikes per neuron per s
    """
    ts,gid = load_sim(path,simulation,(0,sim_time*1000))
    fr = np.zeros([N])
    for ind,i in enumerate(np.arange(1,N+1)):
        fr[ind] = len(ts[np.where(gid==i+1)[0]])/(sim_time)
    
    return fr

def lazy_unit_isi(path,
             simulation,
             N = 1000,
             sim_time = 200
             ):
    """Computes an average firing rate of the NEST simulation
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            N (int): Number of neurons
            sim_time (int): length of the simulation in s.
            
    Returns:
            float64: N of spikes per neuron per s
    """
    ts,gid = load_sim(path,simulation,(0,sim_time*1000))
    isi = np.zeros([N])
    for ind,i in enumerate(np.arange(1,N+1)):
        isi[ind] = np.mean(np.diff(ts[np.where(gid==i+1)[0]]))
    
    return isi

    
def return_sc(path,
            simulation,
            t= (0,5000),
            N=1000,
            bin_size= 20,
            cells = None):
    """Computes an spike count in bins
    Args:
        path (str): Directory of the simulation with "/".
        simulation (str): Simulation name
        t (tuple): simulation time (from,until) in ms
        N (int): Number of neurons
        bin_size (int): ms in one bins
        cells (string): 'ex' - take only excitatory neuron, 'inh' 
                        take only inhibitory neurons, anything else- 
                        take the whole population;
    Returns:
        arr: spike count
    """
    t1,t2 = t 
    ts1,gids1 = load_sim(path,simulation, t)
    if cells == 'ex':
        sc,_=np.histogram(ts1[np.where(gids1<=800)[0]],np.arange(t1,t2+1,bin_size))
    elif cells == 'in':
        sc,_=np.histogram(ts1[np.where(gids1>800)[0]],np.arange(t1,t2+1,bin_size))
    else:
        sc,_=np.histogram(ts1,np.arange(t1,t2+1,bin_size))
    return sc
   
   
    
##############  Burst Identification helpers ##############

class MultiplePeaks(Exception): pass
class NoPeaksFound(Exception): pass

def fwhm(x, y, k=10):
    """
    Determine full-with-half-maximum of a peaked set of points, x and y.

    Assumes that there is only one peak present in the datasset.  The function
    uses a spline interpolation of order k.
    """

    half_max = np.max(y)/2.0
    s = splrep(x, y - half_max, k=k)
    roots = sproot(s)
    #print(roots)
    if len(roots) > 2:
#         raise MultiplePeaks("The dataset appears to have multiple peaks, and "
#                 "thus the FWHM can't be determined.")
        return 0
    elif len(roots) < 2:
        return 0
#         raise NoPeaksFound("No proper peaks were found in the data set; likely "
#                 "the dataset is flat (e.g. all zeros).")
    else:
        return abs(roots[1] - roots[0])
    
    

def burst_times(path,
                simulation,
                bin_size=20,
                t= (0,5000),
                thr = None,
                add_size = False):
    """Finds Burst times with 5 sigma threshold
    sigma  = median(|x|/0.6745)
    
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            N (int): Number of neurons
            binsize(int): hist bins in ms
            t (tuple: int,int): tuple of time start, time stop
            thr (int or None): Threshold in spikes. 
                               if None, estimates the threshold 
                               with 5 sigma
            add_size (Bool): return sizes of bursts additionally
            
    Returns:
            either list: Burst times
            or tuple: Burst times and sizes (if size is True)
            
    """
    t1,t2 = t
    ts,gids = load_sim(path,simulation,t)
    if np.max(ts)<t2/4:
        return []
    sc,_=np.histogram(ts,np.arange(t1,t2,bin_size))
    if thr == None:
        sigma =np.median(sc/0.6745)
        thr = 5*sigma
    indx = np.where(sc>(thr))[0]
    if len(indx)>0:
        indx = indx[np.hstack([np.array([True]),np.diff(indx)>2])]
        if add_size:
            size = [sc[ind] for ind in indx]
    elif add_size:
        size = [0]
    if add_size:
        return (indx*bin_size, size)
    else:
        return indx*bin_size

def burst_size(path,simulation,bin_size=20,lag= 5):
        """Finds Burst sizes
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
            N (int): Number of neurons
            binsize(int): hist bins in ms
            t (tuple: int,int): tuple of time start, time stop
            thr (int or None): Threshold in spikes. 
                               if None, estimates the threshold 
                               with 5 sigma
            
    Returns:
            list: Burst times
            """
        t1,t2 = t
        ts,gids = load_sim(path,simulation,t)
        sc,_=np.histogram(ts,np.arange(t1,t2,bin_size))
        sigma =np.median(sc/0.6745)
        indx = np.where(sc>5*sigma)[0]
        if np.max(ts)<100000:
            return []
        if len(indx)>0:
            indx = indx[np.hstack([np.array([True]),np.diff(indx)>2])]
            size = [sc[ind] for ind in indx]
        else:
            size =[0]
        return size#indx*bin_size#burst_t *bin_size



def bursts(simulation,bin_size=20,lag= 5, verbose=False):
    "Return:  1) a number of bursts \
    2)burst durations calculated as FWHM"
#     if verbose:
#         plt.figure()
    t1,t2 = (0,200000)
    ts,gids = load_ts(simulation,t1,t2)
    sc,_=np.histogram(ts,np.arange(t1,t2,bin_size))
    sigma =np.median(sc/0.6745)
    indx = np.where(sc>5*sigma)[0]
    if np.max(ts)<100000:
        return None, None
    if len(indx)<1:
        return 0, None
    
    if len(indx)>0:
        indx = indx[np.hstack([np.array([True]),np.diff(indx)>2])]
    #lag = 5
    burst = []
    for i in indx:
        if i>=lag:
            n_i = i#np.argmax(sc[i-lag:i+lag])
            if verbose:
                
                plt.plot(sc[(i+n_i)-lag*2:(i+n_i)+lag],'k',alpha= 0.5)
                #plt.plot(np.diff(sc[(i+n_i)-lag*2:(i+n_i)+lag]),'r',alpha= 0.5)
            burst.append(sc[(i+n_i)-lag*2:(i+n_i)+lag])
            
    if verbose:
        burst = align(burst)
        try:
            plt.plot(np.mean(burst,0),linewidth=5)
        except:
            plt.plot(np.mean(pad(burst),0),linewidth=5)
    BurstDuration = np.array(burst_duration(burst))*bin_size
   #BurstDuration= [fwhm(np.arange(0,len(burst[i]))*bin_size,burst[i],k=3) for i in range(len(burst))]
    #[fwhm(np.arange(0,len(burst[i]))*bin_size,burst[i],k=3) for i in range(len(burst))]
    if verbose:
        print(simulation)
        print('Number of bursts %s'%(len(burst)))
        print('Mean burst duration: %s ms'%(np.mean(BurstDuration)))
    return len(burst), np.mean(BurstDuration)

def burst_duration(bursts):
    l_burst = []
    for burst in bursts:
        thr= np.max(burst)*0.1
        indx = np.where(burst>thr)[0]
        l_burst.append(indx[-1]-indx[0])
    return l_burst

def pad(burst):
    max_b= np.max([len(i) for i in burst])
    burst =[np.pad(A,(0,max_b-len(A)), 'constant') for A in burst]
    return burst

def interBurstInterval(sc,indx,thr = 0.05):
    """Finds Inter-Burst Interval with 5% of the maximum to find the burst lengtg
    
    Args:
            indx (str): indices of burst peaks
            size (str): amplitude of burst peaks
            N (int): Number of neurons
            binsize(int): hist bins in ms
            t (tuple: int,int): tuple of time start, time stop
            thr (int or None): Threshold in spikes. 
                               if None, estimates the threshold 
                               with 5 sigma
            size (Bool): return sizes of bursts additionally
            
    Returns:
            either list: Burst times
            or tuple: Burst times and sizes (if size is True)
            
    """
    ibi = []

    for i in range(len(indx)-1):
        #start of the next burst
        size = np.max(sc[indx[i+1]-1000:indx[i+1]+2000])
        start_ = (indx[i+1]-1000)+ np.where(sc[indx[i+1]-1000:indx[i+1]+2000]>(size*thr))[0][0]
        #end of the current burst
        stop = (indx[i]-1000)+np.where(sc[indx[i]-1000:indx[i]+2000]>(size*thr))[0][-1]
        ibi.append(start_-stop)
    return ibi

###### Plotting Functions #########
def out_deg_plot(path,
                simulation,
                out_deg,
                t,
                N,
                u_id=(0,8000),
                marker=1,
                leader_marker = 5,
                inh= False):    
    
    id_1,id_n= u_id
    t1,t2 = t
    ts1,gids1 = read_gdf(path,simulation,t)
    # make the scatter
    #out_deg = (out_deg/np.max(out_deg))*1
    #out_deg = 1.5
    indi = np.argsort(out_deg)
    gid = gids1.copy()
    for ind,i in enumerate(indi):
        gid[np.where(gids1==i+1)[0]] = ind

    #for pos,i in enumerate(indi):
        
    plt.plot(ts1,gid,'.k',markersize =marker)

    #plt.plot(ts1_i,gids1_i,'.r',markersize =4)    
    #plt.figure(figsize=(15,5))
    #plt.plot(ts1,gids1,'.k',markersize =marker)
   
    #print(gids1[np.where(gids1==leader)[0]])
    plt.xlim([t1,t2])
    plt.ylim([id_1,id_n])
    plt.ylabel('units (by outdegree)')

def mark_leader(path,
                simulation,
                leader,
                t,
                N,
                u_id=(0,8000),
                marker=1,
                leader_marker = 5,
                inh= False):    
    id_1,id_n= u_id
    t1,t2 = t
    ts1,gids1 = read_gdf(path,simulation,t)
    #plt.plot(ts1_i,gids1_i,'.r',markersize =4)    
    #plt.figure(figsize=(15,5))
    plt.plot(ts1,gids1,'.k',markersize =marker)
    for n,l in enumerate(leader):
        #plt.plot(ts1[np.where(gids1==l)[0]],gids1[np.where(gids1==l)[0]],'.r',markersize =leader_marker)
        tt = ts1[np.where(gids1==l)[0]]
        idi = np.zeros(len(tt))
        idi[:]=leader[n]
        plt.plot(tt,idi,'.r',markersize =leader_marker)
    #print(gids1[np.where(gids1==leader)[0]])
    plt.xlim([t1,t2])
    plt.ylim([id_1,id_n])
    plt.ylabel('unit id')
    
def hist_leader(simulation,
                leader,
                t1,t2,
                N,
                bin_size,
                u_id=(0,8000),
                marker=1,
                inh= False):    
    id_1,id_n= u_id
    ts1,gids1 = read_gdf('sim/out_deg/',simulation,t1,t2)
    #plt.plot(ts1_i,gids1_i,'.r',markersize =4)    
    #plt.figure(figsize=(15,5))
    ts_l = np.array([])
    for l in leader:
        ts_l = np.hstack([ts_l, ts1[np.where(gids1==l)[0]]])
        
    sc,_=np.histogram(ts_l,np.arange(t1,t2+1,bin_size))
    
    plt.plot(sc/len(leader))
    plt.xticks(np.arange(0,len(sc)+1,10000/bin_size), np.arange(0,(len(sc)+1)*bin_size,10000))
    #sns.despine(trim =40)
    plt.ylabel('rate')
    #print(np.mean(isi_unit_variance(ts,gids,100)))
    plt.xlabel('time (ms)')
    plt.xlim([0,len(sc)])
    
def plot_raster(path,
                simulation,
                t,
                N,
                u_id=(0,8000),
                marker=1,
                inh= False):    
    id_1,id_n = u_id
    t1,t2 = t
    #ts1,gids1 = load_ts(simulation,t1,t2)
    ts1,gids1 = load_sim(path,simulation,t)

    #plt.figure(figsize=(15,5))
    plt.plot(ts1,gids1,'.k',markersize =marker)
    plt.xlim([t1,t2])
    plt.ylim([id_1,id_n])
    plt.ylabel('unit id')
    #sns.despine(trim =40)
    
def plot_noize(path,
                simulation,
               noise_sim,
                t,
                N,
                u_id=(0,8000),
                marker=1,
                inh= False):    
    id_1,id_n = u_id
    t1,t2 = t
    #ts1,gids1 = load_ts(simulation,t1,t2)
    ts1,gids1 = load_sim(path,simulation,t)

    #plt.figure(figsize=(15,5))
    plt.plot(ts1,gids1,'.k',markersize =marker)
    
    ts1,gids1 = load_sim(path,noise_sim,t)

    #plt.figure(figsize=(15,5))
    plt.plot(ts1,gids1,'r+',markersize =marker,alpha =1)
    plt.xlim([t1,t2])
    plt.ylim([id_1,id_n])
    plt.ylabel('unit id')
    #sns.despine(trim =40) 
    
def plot_sc(path,
            simulation,
            t= (0,5000),
            N=1000,
            bin_size= 20,
            **kwargs):
    t1,t2 = t 
    ts1,gids1 = load_sim(path,simulation,t)
    #ts1,gids1 = load_ts(simulation,t1,t2)
    sc,_=np.histogram(ts1,np.arange(t1,t2+1,bin_size))
    sc = sc#/N
    plt.plot(sc, **kwargs)
    plt.xticks(np.arange(0,len(sc)+1,10000/bin_size), np.arange(0,(len(sc)+1)*bin_size,10000))
    sigma =np.median(sc/0.6745)
    plt.plot(np.arange(0,len(sc)),[5*sigma]*len(sc),'--',linewidth =1)
    #sns.despine(trim =40)
    plt.ylabel('rate')
    #print(np.mean(isi_unit_variance(ts,gids,100)))
    plt.xlabel('time (ms)')
    plt.xlim([0,len(sc)])
    
    
def plotSpectrum(y,Fs):
    from scipy import fft, arange
    from numpy import sin, linspace, pi
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    n = len(y) # length of the signal
    k = arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    print(n/2)
    frq = frq[range(int(n/2))] # one side frequency range

    Y = fft(y)/n # fft computing and normalization
    Y = Y[range(int(n/2))]

    plt.plot(frq,abs(Y)) # plotting the spectrum
    plt.xlabel('Freq (Hz)')

    plt.ylabel('|Y(freq)|')
    
    
def norm_hist(x,**kwarg):
    weights = np.ones_like(x)/float(len(x))
    plt.hist(x,weights=weights,**kwarg)
