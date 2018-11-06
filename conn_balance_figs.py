SC1=  []
rate = []
J =[650.]#, 650.,,650.,650.,650.]#[600.,600.,600.,600.]#,450.,470.,480.,490.,495.,500.,550.,600.,700.,900.]#[600.,600.,764.,1208.,2416.]#[1.97,2.6,2.55,3.6,16.0,60.0]# [400.,400.,400.,400.,400.]#
Taus = [90000.]
NI_ = [500,800,2000,2000,3200,3200]#,500,500,800]#,500,500,500,500,500,500,500,500,500,500]#[50,100,200,300,500,700,800,950]#[100,200,500,800,990]#100,200,500,800,
G  =[7.0,4.0,1.0,4.0,0.25,4.0]#[4.0,1.0,4.0,4.0]
U =  0.035
tau_fac = 0.0
scaling = [1.]#[0.2,0.5,0.75,1,1.25,1.5,2,4]
lines_ = np.zeros([len(NI_)])
lines_std= np.zeros([len(NI_)])
ETA = 0.1
sim_time =600000
Peaks = []
for i,NI in enumerate(NI_):
    for tau_rec in Taus:
        #for i in range(len(J)):
        for n,scale in  enumerate(scaling):
            j =J[0]
            NE = int(4000-NI)
            eta = 0.1#ETA[i]#0.2#0.5
            d =  [3.5]
            g =G[i]#(NE/NI)*scale
            path   = 'sim/tsodyks_inh_scaling_Poiss/'
            simulation='tsodyks_j=%s_g=%s_eta=%s_NI=%s_u=%s_tau_rec=%s_tau_fac=%s_Nt=%s'%(j,g,eta,NI,U,tau_rec,tau_fac,NI+NE)
            print(simulation)
            bin_size = 20
            SC1.append(return_sc(path,simulation,(100,sim_time),bin_size=bin_size))#/int((1000-NI_[i])                
            rate.append(lazy_rate(path,simulation,sim_time = (sim_time/1000), N=NI+NE))
#             thr =4000#NE# NE#5*np.std(SC1[-1][0:10000])
#             print(thr*1.5)
#             tau =2000/bin_size
#             signal  = np.convolve(SC1[i],np.exp(-np.arange(10000/bin_size)/tau),'valid')
            tau =2000/bin_size
#             signal  = np.convolve(SC1[-1],np.exp(-np.arange(10000/bin_size)/tau),'valid')
#             Sca = signal/np.mean(signal)
#             peakTimes_,peakAmp = find_peaks(Sca,height=3*np.std(Sca),width=1,distance=50)#4000-NI_[i]

            peakTimes_,peakAmp = find_peaks(SC1[-1],height=thr,width=1,distance=200)#4000-NI_[i]
            if len(peakTimes_)>1:
                lines_[i] = np.nanmean(np.diff(peakTimes_))#len(peakTimes)/300#len(peakTimes_)#
                Peaks.append(peakTimes_)
                lines_std[i] = np.std(np.diff(peakTimes_))#len(peakTimes)/300#len(peakTimes_)#
plt.figure(figsize=(10,5))
plt.plot(SC1[2]/(4000-NI_[0]))
plt.plot((SC1[3]/(4000-NI_[1]))-5)
sns.despine()
plt.title('50% inh. neurons')
plt.yticks([0,-5],['PSP-balance','Balanced input'])
plt.tight_layout()
plt.savefig('figs/conn_balance50.png')

# plt.plot((SC1[2]/(4000-NI_[2]))-10)
# plt.plot((SC1[3]/(4000-NI_[3]))-25)
# plt.plot((SC1[4]/(4000-NI_[4]))-35)
# plt.plot((SC1[5]/(4000-NI_[5]))-45)
# plt.xlim([0,5000])
# plt.ylim([-6,10])


plt.plot((lines_std[1:]/lines_[1:]),'o')
sns.despine()
plt.ylabel('Coff. of variation')
plt.xticks(np.arange(0,5),['20','50% B-PSP','50 B-Conn','80% B-PSP','80 B-Conn',],fontsize = 12)
plt.tight_layout()
plt.savefig('figs/conn_balance_var.png')


plt.plot(((lines_[1:]*20)/1000),'o')
plt.ylabel('IBI')
sns.despine()
plt.xticks(np.arange(0,5),['20','50% B-PSP','50 B-Conn','80% B-PSP','80 B-Conn',],fontsize = 12)
plt.tight_layout()
plt.savefig('figs/conn_balance_mean.png')
