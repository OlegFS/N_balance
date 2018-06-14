
import nest
import nest.raster_plot
import numpy as np
import matplotlib.pyplot as plt
from func.brunel_meta import *
from func.helpers import *
j=0.5
J= [0.1,0.5]
ind = 0
meanV = np.zeros(len(J))
stdV = np.zeros(len(J))
rate = np.zeros(len(J))
simulation = 'VoltTest_j=%s_n=1000'%(j)
s_voltage = read_voltage('sim/voltage_test/',
                     simulation,
                     1000,
                     (0,2180),
                     nice_format=True)

#np.save('sim/voltage_test/sampled_voltage',s_voltage)
print('pass1')
mm = []
print('pass2')
sstd = []
print('pass3')
for i in range(0,s_voltage.shape[1]):
    mm.append(np.mean(s_voltage[:,i]))
    sstd.append(np.std(s_voltage[:,i]))
print('pass4')
meanV[ind] = np.mean(mm[5000:])
print('pass5')
stdV[ind] = np.mean(sstd[5000:])
print('pass6')
np.save('sim/voltage_test/sampled_meanV',meanV)
print('pass7')
np.save('sim/voltage_test/sampled_stdV',meanV)

