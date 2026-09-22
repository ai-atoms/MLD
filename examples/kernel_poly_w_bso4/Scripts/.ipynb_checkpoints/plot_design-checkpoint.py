import numpy as np
import seaborn as sns
import matplotlib.pylab as plt

uniform_data = np.loadtxt("design_matrix.dat", dtype='float')
data=uniform_data[:,1:-4]
#print(data)
#ax = sns.heatmap(data, linewidth=0.0, center=0, cmap="YlGnBu", robust=True)
#plt.show()
data=uniform_data[1:400,40:-4]
ax = sns.heatmap(data, linewidth=0.0, center=0, cmap="YlGnBu", robust=True)
plt.show()
