# -*- coding: utf-8 -*-
import time
from numpy import exp
import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import networkx as nx
import pandas
import h5py as h5py

def FixedInGraph(N,degrees,replace = True):
    """Generate connectivity with fixed in-degree
    Args:
            degree (int): target in-degree
            N (int): number of neurons
            replace (bool): double connections between two nodes 
    Returns:
            arr: connectivity array [source, target] """
        
    conn = np.zeros([N*degrees,2])
    nodes = np.arange(0,N)
    ii = 0
    for i in range(N):
    #     for deg in range(degrees):
        source = np.random.choice(nodes,size =degrees,replace=replace)
        conn[ii:ii+len(source),0] = source
        conn[ii:ii+len(source),1] = i
        ii+=len(source) 
    return conn


def get_deg(nodes = 10000,edges=None):
    """get in and out-degrees and a networkx class for a list of edges
    Args:
            nodes (int): number of nodes
            edges (arr): dim [source, target]: list of edges
    Returns:
            1) arr: in_degrees per node
            2) arr: out_degrees per node
            3) networkx class"""
    
    D= nx.MultiDiGraph()
    D.add_nodes_from(np.arange(1,nodes))
    D.add_edges_from(edges);
    return np.array(D.in_degree())[:,1],np.array(D.out_degree())[:,1],D


def fixedDiGraph(conn_,out_deg,target_degree,N):
    """Generate connectivity with fixed in and out-degree
    Args:
            conn_ (arr): dim: Nx2 [source, target]:  Set of edges (same as in networkx) 
            out_deg (arr):dim N: out_degrees per neurons
            target_degree (int): target out-degree
            N (int): number of neurons
    Returns:
            arr: connectivity array [source, target] """
        
    d_o = (out_deg[0:N]-target_degree) #outdegree difference
    # [1:NE] - ignores the first dummy element 
    pos_ind = np.where(d_o>0)[0] # positive difference (remove connections)
    neg_ind = np.where(d_o<0)[0] # negative difference (add connection)

    add_ind= []
    for ind in pos_ind:
        while len(conn_[conn_[:,0]==ind])>target_degree: # +1 might need to go away
            add_ind.append(conn_[conn_[:,0]==ind][-1,1])
            # delete a raw
            ind_d = np.where(conn_[:,0]==ind)[0][-1]
            conn_ = np.delete(conn_,ind_d,0)
    ii = 0
    for ind in neg_ind:
        while len(conn_[conn_[:,0]==ind])<target_degree:
            conn_ = np.vstack((conn_,[ind, add_ind[ii]]))
            ii+=1

    return conn_



############ Connectivity helpers################

def read_out_deg(path = 'sim/out_deg/', simulation= 'test-out_deg'):
    """Reads outdegrees saved elsewhere
    
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
    Returns:
            arr: out degrees
    """
    with h5py.File('sim/out_deg/'+'/'+simulation+'-out_deg','r') as f:                                                                                                                      
        out_deg=np.array(f['out_deg']) 
    return out_deg


def read_in_deg(simulation):
    with h5py.File('sim/out_deg/'+'/'+simulation+'-in_deg','r') as f:                                                                                                                      
        out_deg=np.array(f['in_deg']) 
    return out_deg

def read_conn(path, simulation):
    """Reads connectivity saved from NEST
    
    Args:
            path (str): Directory of the simulation with "/".
            simulation (str): Simulation name
    Returns:
            arr: connectivity array [source, target]
    """
    with h5py.File(path+simulation+'-conn','r') as f: 
        print(path+simulation+'-conn' )
        out_deg=np.array(f['conn']) 
    return out_deg



'''
Plot weight matrices example
----------------------------

This example demonstrates how to extract the connection strength
for all the synapses among two populations of neurons and gather
these values in weight matrices for further analysis and visualization.

All connection types between these populations are considered, i.e.,
four weight matrices are created and plotted.
'''

'''
First, we import all necessary modules to extract, handle and plot
the connectivity matrices
'''

'''
We now specify a function which takes as arguments lists of neuron gids
corresponding to each population
'''


def plot_weight_matrices(E_neurons, I_neurons):
    '''
    Function to extract and plot weight matrices for all connections
    among E_neurons and I_neurons
    '''

    '''
    First, we initialize all the matrices, whose dimensionality is
    determined by the number of elements in each population
    Since in this example, we have 2 populations (E/I), 2^2 possible
    synaptic connections exist (EE, EI, IE, II)
    '''

    import numpy as np
    import pylab
    import nest
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    W_EE = np.zeros([len(E_neurons), len(E_neurons)])
    W_EI = np.zeros([len(I_neurons), len(E_neurons)])
    W_IE = np.zeros([len(E_neurons), len(I_neurons)])
    W_II = np.zeros([len(I_neurons), len(I_neurons)])

    '''
    Using `GetConnections`, we extract the information about all the
    connections involving the populations of interest. `GetConnections`
    returns a list of arrays (connection objects), one per connection.
    Each array has the following elements:
    [source-gid target-gid target-thread synapse-model-id port]
    '''

    a_EE = nest.GetConnections(E_neurons, E_neurons)

    '''
    Using `GetStatus`, we can extract the value of the connection weight,
    for all the connections between these populations
    '''

    c_EE = nest.GetStatus(a_EE, keys='weight')

    '''
    Repeat the two previous steps for all other connection types
    '''

    a_EI = nest.GetConnections(I_neurons, E_neurons)
    c_EI = nest.GetStatus(a_EI, keys='weight')
    a_IE = nest.GetConnections(E_neurons, I_neurons)
    c_IE = nest.GetStatus(a_IE, keys='weight')
    a_II = nest.GetConnections(I_neurons, I_neurons)
    c_II = nest.GetStatus(a_II, keys='weight')

    '''
    We now iterate through the list of all connections of each type.
    To populate the corresponding weight matrix, we begin by identifying
    the source-gid (first element of each connection object, n[0])
    and the target-gid (second element of each connection object, n[1]).
    For each gid, we subtract the minimum gid within the corresponding
    population, to assure the matrix indices range from 0 to the size of
    the population.

    After determining the matrix indices [i, j], for each connection
    object, the corresponding weight is added to the entry W[i,j].
    The procedure is then repeated for all the different connection types.
    '''

    for idx, n in enumerate(a_EE):
        W_EE[n[0] - min(E_neurons), n[1] - min(E_neurons)] += c_EE[idx]
    for idx, n in enumerate(a_EI):
        W_EI[n[0] - min(I_neurons), n[1] - min(E_neurons)] += c_EI[idx]
    for idx, n in enumerate(a_IE):
        W_IE[n[0] - min(E_neurons), n[1] - min(I_neurons)] += c_IE[idx]
    for idx, n in enumerate(a_II):
        W_II[n[0] - min(I_neurons), n[1] - min(I_neurons)] += c_II[idx]

    '''
    We can now specify the figure and axes properties. For this specific
    example, we wish to display all the weight matrices in a single
    figure, which requires us to use ``GridSpec`` (for example)
    to specify the spatial arrangement of the axes.
    A subplot is subsequently created for each connection type.
    '''

    fig = pylab.figure()
    fig.suptitle('Weight matrices', fontsize=14)
    gs = gridspec.GridSpec(4, 4)
    ax1 = pylab.subplot(gs[:-1, :-1])
    ax2 = pylab.subplot(gs[:-1, -1])
    ax3 = pylab.subplot(gs[-1, :-1])
    ax4 = pylab.subplot(gs[-1, -1])

    '''
    Using ``imshow``, we can visualize the weight matrix in the corresponding
    axis. We can also specify the colormap for this image.
    '''

    plt1 = ax1.imshow(W_EE, cmap='jet')

    '''
    Using the ``axis_divider`` module from ``mpl_toolkits``, we can
    allocate a small extra space on the right of the current axis,
    which we reserve for a colorbar.
    '''

    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", "5%", pad="3%")
    pylab.colorbar(plt1, cax=cax)

    '''
    We now set the title of each axis and adjust the axis subplot parameters
    '''

    ax1.set_title('W_{EE}')
    pylab.tight_layout()

    '''
    Finally, the last three steps are repeated for each synapse type
    '''

    plt2 = ax2.imshow(W_IE)
    plt2.set_cmap('jet')
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", "5%", pad="3%")
    pylab.colorbar(plt2, cax=cax)
    ax2.set_title('W_{EI}')
    pylab.tight_layout()

    plt3 = ax3.imshow(W_EI)
    plt3.set_cmap('jet')
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", "5%", pad="3%")
    pylab.colorbar(plt3, cax=cax)
    ax3.set_title('W_{IE}')
    pylab.tight_layout()

    plt4 = ax4.imshow(W_II)
    plt4.set_cmap('jet')
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", "5%", pad="3%")
    pylab.colorbar(plt4, cax=cax)
    ax4.set_title('W_{II}')
    pylab.tight_layout()


