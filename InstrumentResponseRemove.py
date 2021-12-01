# Define math defaults
#from __future__ import division #allows real devision without rounding

# Retrieve modules needed
from obspy.core import read
import numpy as np
import matplotlib.pyplot as plt

#%% Choose and import data
str = read('your file path')
print(str) #show imported data
print(str[1].stats) #show stats for the intended trace



# create dictionary of poles and zeros,come from the paz file for each station

TrillC = {'gain': 800.0,
        'poles': [complex(-3.691000e-02,3.712000e-02),
                  complex(-3.691000e-02,-3.712000e-02),
                  complex(-3.739000e+02,4.755000e+02),
                  complex(-3.739000e+02,-4.755000e+02),
                  complex(-5.884000e+02,1.508000e+03),
                  complex(-5.884000e+02,-1.508000e+03)],
        'sensitivity': 8.184000E+11,
        'zeros': [0 -4.341E+02]}


#%% Remove instrument response

str1 = str.copy() #make a copy of data
str1.simulate(paz_remove=TrillC, paz_simulate=None, water_level=60.0)
str1_m = str1.merge()

for npp in range(len(str)):
    plt.figure()
    plt.subplot(2,1, 1)
    plt.plot(str[npp].data, label=str[npp])
    plt.legend(loc=0)

    plt.subplot(2,1, 2)
    plt.plot(str1[npp].data, label=str1[npp])
    plt.legend(loc=0)

print("Instrument Response Removed")

tr1 = str1_m[1]
tr2= str1_m[2]

#tr_orig = st_orig[0]

t = np.arange(tr1.stats.npts) / tr1.stats.sampling_rate

np.savetxt('name.csv', np.transpose([t, tr1, tr2]),fmt='%.3g')

plt.show()
