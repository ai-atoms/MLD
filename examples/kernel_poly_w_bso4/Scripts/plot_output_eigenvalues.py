import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import math
import scipy
from sklearn.metrics import r2_score
from sklearn import datasets, linear_model
from scipy.stats import gaussian_kde

"""
import numpy as np
import matplotlib.pyplot as plt

# Generate fake data
x = np.random.normal(size=1000)
y = x * 3 + np.random.normal(size=1000)

# Calculate the point density
xy = np.vstack([x,y])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

fig, ax = plt.subplots()
ax.scatter(x, y, c=z, s=50, edgecolor='')
"""


f1=open('train_sigma.eigenvalues','r')
data_int=np.loadtxt(f1, usecols=[0])
train=data_int[:]
#cs_y_1=CubicSpline(x_1,y_1)
f1.close()
f2=open('full_sigma.eigenvalues','r')
data_int=np.loadtxt(f2, usecols=[0])
full=data_int[:]
f2.close()

"""
temp_min=x_1[0]
temp_max=x_1[-1]
x_all=np.linspace(temp_min, temp_max, num=1000, endpoint=True)
"""


#ydata_new=np.log(data_new[:,1])
#yerr_new=np.log(data_new[:,2])

#yerrp=np.log(data_new[:,1]+data_new[:,2]) - np.log(data_new[:,1])
#yerrm=np.log(data_new[:,1]) -np.log(data_new[:,1]-data_new[:,2])

#Method to use for optimization. See least_squares for more details.
#Default is 'lm' for unconstrained problems and
#Default is 'trf' if bounds are provided.
#The method 'lm' won't work when the number of observations is less
#than the number of variables, use 'trf' or 'dogbox' in this case.

fig=plt.figure()
ax=fig.add_subplot(111)



#ax.plot(x_1**0.5, le_1, 'o',  color='salmon', markersize=1)
ax.plot(train,  'o',  color='salmon', markersize=4, label=r"$\Sigma_{bulk}$")
ax.plot(full,   'o',  color='blue', markersize=4, label=r"$\Sigma_{full}$")

#ax.set_xlim([0,70.0])
#ax.set_ylim([0,3.0])
ax.set_yscale('log')
#ax.set_xlabel(r" reaction coordinate ($a0 \sqrt{2/3}$ units)", fontsize=20)
ax.set_xlabel(r" $i^{th}$ eigenvalue ", fontsize=20)
#ax.set_ylabel(r" $d_{MAH}$ with $\Sigma_{full}$ ", fontsize=20) #, color='slateblue')
ax.set_ylabel(r" $ \lambda_i$ ", fontsize=20) #, color='slateblue')
#ax.errorbar(xana, yana, yerr=[yerrm,yerrp], fmt='o',markersize=10, color='red', ecolor='red', capthick=2, label='anharmonic DFT')
#ax.tick_params(axis='x', labelsize='16')
ax.tick_params(axis='y', labelsize='14')#,colors='#88BBF0')
ax.tick_params(axis='x', labelsize='14')#,colors='#88BBF0')

# ad the second xaxis on the top of the graph.
ax.legend(loc='best', fontsize=20)



#majorLocator = matplotlib.ticker.MultipleLocator(1000.0)
#majorFormatter = matplotlib.ticker.FormatStrFormatter('%5.1f')
#minorLocator = matplotlib.ticker.AutoMinorLocator(5)
#ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
#ax.yaxis.set_minor_locator(minorLocator)


majorLocator = matplotlib.ticker.MultipleLocator(10.0)
majorFormatter = matplotlib.ticker.FormatStrFormatter('%0.0f')
minorLocator = matplotlib.ticker.AutoMinorLocator(5)
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.xaxis.set_minor_locator(minorLocator)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='black')





plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=7)
plt.tick_params(which='minor', length=4, color='black')





defsize=fig.get_size_inches()
defsize[0]=4.8
defsize[1]=4.8
print(defsize)
fig.set_size_inches( (defsize[0]*1.5,defsize[1]*1.5) )
fig.savefig('test.png')
fig.savefig('test.eps')

plt.show()
