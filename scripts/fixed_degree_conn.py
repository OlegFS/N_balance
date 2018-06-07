from func.brunel import delta_brunel
from func.brunel_meta import meta_brunel,fixedOutDeg_brunel
from func.helpers import *

"""Generate and save a connectivity with fixed in and out degree
Script works only for 1000 units network. 
TODO: test"""


# Load a sample distibutuon with fixed indegreee
g= 4.0
eta =0.50
j = 0.88*1.5
# t1 = 184000
# t2 = 197000
t1,t2 = 0,400000
# simulation = Path('sim/delta_coupled/tau40/brunel_esp_ex_8.5_g_0.5_test01dt_delta_[0.85, 0.86, 0.87, 0.88, 0.89, 0.9].npy')
simulation ='brunel_esp_ex_%s_g_%s_test01dt_delta_%s-1000-'%(g,eta,j)
conn = read_conn(simulation)
#in_deg = read_in_deg(simulation)
in_deg,out_deg,D = get_deg(1000,conn)


conn_= conn.copy()
conn_ = conn_[conn_[:,0]<=1000,:]
conn_ = conn_[conn_[:,1]<=1000,:]


# Start with the excitatory neurons
d_o = (out_deg[0:800]-101) #outdegree difference
pos_ind = np.where(d_o>0)[0] # positive difference (remove connections)
neg_ind = np.where(d_o<0)[0] # Negative difference (add connection)

add_ind= []
for ind in pos_ind:
    while len(conn_[conn_[:,0]==ind+1])>100:
        add_ind.append(conn_[conn_[:,0]==ind+1][-1,1])
        # delete a raw
        ind_d = np.where(conn_[:,0]==ind+1)[0][-1]
        conn_ = np.delete(conn_,ind_d,0)
ii = 0
for ind in neg_ind:
    while len(conn_[conn_[:,0]==ind+1])<100:
        conn_ = np.vstack((conn_,[ind+1, add_ind[ii]]))
        ii+=1
        
 # Inhibitory neurons
d_o = (out_deg[800:1000]-101) #outdegree difference
pos_ind = np.where(d_o>0)[0]+800 # positive difference (remove connections)
neg_ind = np.where(d_o<0)[0]+800 # Negative difference (add connection)

add_ind= []
for ind in pos_ind:
    while len(conn_[conn_[:,0]==ind+1])>100:
        add_ind.append(conn_[conn_[:,0]==ind+1][-1,1])
        # delete a raw
        ind_d = np.where(conn_[:,0]==ind+1)[0][-1]
        conn_ = np.delete(conn_,ind_d,0)
ii = 0
for ind in neg_ind:
    while len(conn_[conn_[:,0]==ind+1])<100:
        conn_ = np.vstack((conn_,[ind+1, add_ind[ii]]))
        ii+=1
        
       

        #add a raw
        
in_deg_,out_deg_,D_ = get_deg(1000,conn_)
#add test 


# Sort, clear and save 
conn__ = conn_[conn_[:,0].argsort()]
conn__=conn__[conn__[:,0]<1001,:] 
np.save('clear',conn__)
    