import matplotlib as mpl
mpl.use('TkAgg')
import  matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import math
import scipy
from sklearn.metrics import r2_score
from sklearn import datasets, linear_model
from scipy.stats import gaussian_kde
import seaborn as sns
from sklearn.metrics import r2_score
from sklearn import datasets, linear_model
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


f1=open('full_sigma.pca','r')
data_int=np.loadtxt(f1, usecols=[0,1,2,3,4])
pca0=data_int[:,0]
pca1=data_int[:,1]
pca2=data_int[:,2]
db_class = data_int[:,4]
#cs_y_1=CubicSpline(x_1,y_1)
f1.close()



f1=open('train_sigma.pca','r')
data_int=np.loadtxt(f1, usecols=[0,1,2,3,4])
spca0=data_int[:,0]
spca1=data_int[:,1]
spca2=data_int[:,2]
db_class = data_int[:,4]
#cs_y_1=CubicSpline(x_1,y_1)
f1.close()

f1=open('kernel_full_sigma.pca','r')
data_int=np.loadtxt(f1, usecols=[0,1,2,3,4])
kfpca0=data_int[:,0]
kfpca1=data_int[:,1]
kfpca2=data_int[:,2]
#db_class = data_int[:,4]
#cs_y_1=CubicSpline(x_1,y_1)
f1.close()

f1=open('kernel_train_sigma.pca','r')
data_int=np.loadtxt(f1, usecols=[0,1,2,3,4])
ktpca0=data_int[:,0]
ktpca1=data_int[:,1]
ktpca2=data_int[:,2]
#db_class = data_int[:,4]
#cs_y_1=CubicSpline(x_1,y_1)
f1.close()




"""
temp_min=x_1[0]
temp_max=x_1[-1]
x_all=np.linspace(temp_min, temp_max, num=1000, endpoint=True)
"""

print(db_class.astype(int))

cmm=db_class.astype(int)

fig=plt.figure()
sns.reset_orig()  # get default matplotlib styles back
clrs = sns.color_palette('husl', n_colors=20)

for i in range(1,3):
    if (i==1):
        labelg=r"$\Sigma_{train}$"
        x=spca0
        y=spca1
        xk=ktpca0
        yk=ktpca1

    if (i==2):
        labelg=r"$\Sigma_{full}$"
        x=pca0
        y=pca1
        xk=kfpca0
        yk=kfpca1
    ax=fig.add_subplot(1,2,i)
    #ax.plot(x_1**0.5, le_1, 'o',  color='salmon', markersize=1)
    ax.plot(xk, yk,  'o', color="salmon")
    scatter=ax.scatter(x, y,   c=cmm,  s=2, cmap="viridis", label=labelg)
    #legend1 = ax.legend(*scatter.legend_elements(),
    #                loc="lower left", title="Classes")
    #ax.add_artist(legend1)

#ax.set_xlim([0,70.0])
#ax.set_ylim([0,3.0])
#ax.set_yscale('log')
#ax.set_xlabel(r" reaction coordinate ($a0 \sqrt{2/3}$ units)", fontsize=20)
    ax.set_xlabel(r" $\sigma_1$  ", fontsize=20)
#ax.set_ylabel(r" $d_{MAH}$ with $\Sigma_{full}$ ", fontsize=20) #, color='slateblue')
    ax.set_ylabel(r" MCD$_1$ ", fontsize=20) #, color='slateblue')
#ax.errorbar(xana, yana, yerr=[yerrm,yerrp], fmt='o',markersize=10, color='red', ecolor='red', capthick=2, label='anharmonic DFT')
#ax.tick_params(axis='x', labelsize='16')
    ax.tick_params(axis='y', labelsize='14')#,colors='#88BBF0')
    ax.tick_params(axis='x', labelsize='14')#,colors='#88BBF0')

# ad the second xaxis on the top of the graph.
    ax.legend(loc='best', fontsize=20)
sns.set(style="white", color_codes=True)
sns.jointplot(x=pca0, y=pca1,  kind='scatter');
#majorLocator = matplotlib.ticker.MultipleLocator(1000.0)
#majorFormatter = matplotlib.ticker.FormatStrFormatter('%5.1f')
#minorLocator = matplotlib.ticker.AutoMinorLocator(5)
#ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
#ax.yaxis.set_minor_locator(minorLocator)


#majorLocator = matplotlib.ticker.MultipleLocator(10.0)
#majorFormatter = matplotlib.ticker.FormatStrFormatter('%0.0f')
#minorLocator = matplotlib.ticker.AutoMinorLocator(5)
#ax.xaxis.set_major_locator(majorLocator)
#ax.xaxis.set_major_formatter(majorFormatter)
#ax.xaxis.set_minor_locator(minorLocator)
#plt.tick_params(which='both', width=2)
#plt.tick_params(which='major', length=7)
#plt.tick_params(which='minor', length=4, color='black')

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
