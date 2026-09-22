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

f1=open('train_sigma.stat','r')
data_int=np.loadtxt(f1)
dd = data_int
"""
d0=data_int[:,0]
d1=data_int[:,1]
d2=data_int[:,2]
d3=data_int[:,3]
d4=data_int[:,4]
d5=data_int[:,5]
"""
le=dd[:,6]
#d5=data_int[:,5]
#cs_y_1=CubicSpline(x_1,y_1)
f1.close()


f1=open('train_sigma.mcd','r')
data_int=np.loadtxt(f1, usecols=[0,1,2,3,4])
x_1=data_int[:,2]
le_1=data_int[:,3]


fig=plt.figure()

"""
# Calculate the point density
xy = np.vstack([xi,yi])
z = gaussian_kde(xy)(xy)

# Sort the points by density, so that the densest points are plotted last
idx = z.argsort()
x, y, z = xi[idx], yi[idx], z[idx]
"""
yref=le_1

#ax.plot(x_1**0.5, le_1, 'o',  color='salmon', markersize=1)
regr = linear_model.LinearRegression()
for i in range(0,6):
    ax=fig.add_subplot(1, 6, i+1)
    xi=dd[:,i]**0.5
    ax.plot(xi, yref, 'o',  color='salmon', markersize=1)

    xx=np.array(xi).reshape(-1,1)
    yy=np.array(yref).reshape(-1,1)
    regr.fit(xx,yy)
    score=regr.score(xx,yy)
    print("regression score", score)
    ypred=regr.predict(xx)
    #print("R2 score:", r2_score(xx, yy, multioutput='raw_values'))
    ax.plot(xx, ypred, '-',  color='lightsteelblue', markersize=1)
    ax.text(0, 2.5, r"$R^2$=%0.2f"%score)
    ax.set_xlabel(r" $d_{MCD}$-%s"%i, fontsize=10)
    ax.set_ylim([-0.4,3.0])
    if  (i==0):
        ax.set_ylabel(r" $ \epsilon_l$ ", fontsize=20) #, color='slateblue')
    if (i==0):
        ax.text(0, 2.75, r"$\Sigma \langle i_{\star} | m\rangle^2$" )

    if (i==1):
        ax.text(0, 2.75, r"$\Sigma \langle i_{\star} | m\rangle^2/\lambda_m$" )


    if (i==2):
        ax.text(0, 2.75, r"$\Sigma \langle i_{\star} | m\rangle^2/\lambda_m^2$" )

    if (i==3):
        ax.text(0, 2.75, r"$\Sigma \langle i_{\star} | m\rangle^2/\lambda_m^3$" )


    if (i==4):
        ax.text(0, 2.75, r"$\Sigma \langle i_{\star} | m\rangle^2 \lambda_m$" )

    if (i==5):
        ax.text(0, 2.75, r"$\Sigma \langle i_{\star} | m\rangle^2 \lambda_m^2$" )

    ax.tick_params(axis='y', labelsize='8')#,colors='#88BBF0')
    ax.tick_params(axis='x', labelsize='8')#,colors='#88BBF0')

#ax.plot(temp_WEAM2, a0_WEAM2/a0_WEAM2[0]-1, 'o', markersize=3, color='slateblue')


#ax.plot(xinvtempdata, ydata, 'o',markersize=10, color='#88BBF0')
#ax.plot(xinvtempdata_new, ydata_new, 'o',  markersize=1, color='#88BBF0')

#ax2=ax.twinx()
#ax2.set_ylabel(r" $\Delta V_b$ (meV)", color='seagreen', fontsize='20')
#ax2.tick_params(axis='y', colors='seagreen',labelsize='16')

#xerrm=1.0/(temp_xc-0.5*temp_xc_error)-1.0/temp_xc
#xerrp=1.0/temp_xc-1.0/(temp_xc+0.5*temp_xc_error)
#ax2.set_zorder(ax.get_zorder()-1)
#ax2.patch.set_visible(False)
#ax2.errorbar(1.0/temp_xc, 1000.0*deltav_m, xerr=[xerrm,xerrp], yerr=1000.0*deltav_m_error, fmt='o',markersize=8, alpha=0.4, color='seagreen', ecolor='seagreen', capthick=1)

#ax.errorbar(xinvtempdata_new, ydata_new, yerr=[yerrm,yerrp], fmt='o',markersize=10, color='#88BBF0', ecolor='black', capthick=2)



#ax.set_xlim([min(x_1),max(x_1)])
#ax.set_xlim([0,70.0])
#ax.set_ylim([0.0, a0_WEAM2[-1]/a0_WEAM2[0]-1])
#ax.set_ylim([min(y_1)-ene_ref, max(y_1)-ene_ref])
#ax.set_ylim([min(y_1), max(y_1)])
#ax.set_ylim([min(y_1), max(y_1)])
#ax.set_ylim([0, 4])

#xana=[2000.0, 3500.0]
#yana=[1.01, 1.023]
#yerrm=[0.0005, 0.0015]
#yerrp=[0.0005, 0.0015]


#ax.set_xlabel(r" reaction coordinate ($a0 \sqrt{2/3}$ units)", fontsize=20)
#ax.set_ylabel(r" $d_{MAH}$ with $\Sigma_{full}$ ", fontsize=20) #, color='slateblue')
#ax.errorbar(xana, yana, yerr=[yerrm,yerrp], fmt='o',markersize=10, color='red', ecolor='red', capthick=2, label='anharmonic DFT')
#ax.tick_params(axis='x', labelsize='16')

# ad the second xaxis on the top of the graph.
ax.legend(loc='best')


#majorLocator = matplotlib.ticker.MultipleLocator(1.0)
#majorFormatter = matplotlib.ticker.FormatStrFormatter('%5.1f')
minorLocator = matplotlib.ticker.AutoMinorLocator(5)
#ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)


#majorLocator = matplotlib.ticker.MultipleLocator(10.0)
#majorFormatter = matplotlib.ticker.FormatStrFormatter('%0.1f')
minorLocator = matplotlib.ticker.AutoMinorLocator(5)
#ax.xaxis.set_major_locator(majorLocator)
#ax.xaxis.set_major_formatter(majorFormatter)
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
