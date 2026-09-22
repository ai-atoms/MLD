import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit
from scipy import stats
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error


#-------------------------------------
#   Some inputs that can be modified  #
#-------------------------------------

#plot_name = 'g2-bi-SO4, 60-dim, 30xradial, jmax=3, eta_g2_max=2.0 '
#plot_name = 'milady$^2$ , B-SO4, dim_fit 1641  '

# plot depending on clases only 'class' or classes en sub_classes 'full'
type = 'class'

# choice of map color for plot
map_color = 'rainbow'

# Names of files with E, F, S to plot
efile1, efile2 = 'train_energy.out','test_energy.out'
ffile1, ffile2 = 'train_force.out', 'test_force.out'
sfile1, sfile2 = 'train_stress.out', 'test_stress.out'

# Names of the DB-classes and associated colors
#data_names =  ['bcc elastic', 'bcc distorted',    'MD 800 K ', '$V_1$ ', 'surfaces', '$\gamma$-surfaces', 'dislo']
#data_colors = ['green',        'navy', 'turquoise',       'darkorange', 'yellow',      'crimson', 'gray']


stress_calc = (os.path.exists(sfile1) and os.path.exists(sfile2))
#stress_calc = False



#-------------------------------------
#    Define some useful  functions  #
#-------------------------------------

class color_label :

    def __init__(self):
        self.dic = {}
        self.nb_color = 0

    def read_db_model_in_full(self,path_db):
        sub_classe = []
        name = ''
        r = open(path_db,'r').readlines()
        for l in r :
            first=False
            if not l.startswith('##'):
                print(l)
                if l.startswith('#'):
                    if l.startswith('# @'):
                        sub_classe = []
                        data_line = l.split()
                        print(len(data_line), data_line)
                        for k in range(3,6):
                            if data_line[k] != '0':
                                sub_classe.append(data_line[k])
                            else :
                                sub_classe.append('/')
                        name = l.split('#')[-1][1:-1]
                        first=True
                    else:
                        sub_classe = []
                        data_line = l.split()
                        name = l.split('#')[-1][1:-1]


                else :
                    if not l.startswith('#'):
                        tmp = l.split()
                        if not '%s_%s'%(tmp[0],tmp[1]) in self.dic.keys():
                            if tmp[1] == '000':
                                label = r''+name
                            else :
                                label = r''+name+' ('
                                for compt, id in enumerate(tmp[1]):
                                    if (first) :
                                        if sub_classe[compt] != '/' :
                                            label += '%s %s, '%(id,sub_classe[compt])
                                label = label[:-2]+')'
                            new_dic = {'%s_%s'%(tmp[0],tmp[1]):[label,None]}
                            self.dic.update(new_dic)

        return


    def read_db_model_in_class(self,path_db):
        r = open(path_db,'r').readlines()
        for l in r :
            if not l.startswith('##'):
                print(l)
                if l.startswith('# @'):
                    name = l.split('#')[-1][1:-1]

                else :
                    if not l.startswith('#'):
                        tmp = l.split()
                        if not '%s'%(tmp[0]) in self.dic.keys():
                            label = r''+name
                            new_dic = {'%s'%(tmp[0]):[label,None]}
                            self.dic.update(new_dic)
        return

    def update_color(self,map_color):
        self.nb_color = len(self.dic.keys())
        Cmap = plt.get_cmap(map_color)
        color_key = [Cmap(v) for v in np.linspace(0,1,self.nb_color) ]
        for compt, key in enumerate(self.dic.keys()):
            self.dic[key][1] = color_key[compt]

        return

    def extract_list_color_label(self):
        list_color = []
        list_label = []
        for key in self.dic.keys():
            list_label.append(self.dic[key][0])
            list_color.append(self.dic[key][1])

        return list_label, list_color
################################################################################

#def read_x_y_class(filename):
#    if (os.path.exists(filename)):
#      x = np.loadtxt(filename, usecols = (0))
#      y = np.loadtxt(filename, usecols = (1))
#      cl = np.loadtxt(filename, usecols = (2))
#    else:
#      x=[0]
#      y=[0]
#      cl="01"
#    return x, y, cl

def read_x_y_class_new(filename):
    if (os.path.exists(filename)):
      x = np.loadtxt(filename, usecols = (0))
      y = np.loadtxt(filename, usecols = (1))
      cl = np.genfromtxt(filename, usecols = (3) ,dtype='str')
      if type == 'full':
          cl = [c.split('/')[1][:6] for c in cl ]
      if type == 'class':
          cl = [c.split('/')[1][:2] for c in cl ]

    else:
      print('!file %s not found!'%(filename))
      exit(0)
      x=[0]
      y=[0]
      cl="01_000"
    return x, y, cl


#def point_colors(cl_array, data_colors_array):
#    colors = []
#    for n in range (len(cl_array)):
#        indx=int(cl_array[n])-1
#        colors.append(data_colors_array[indx])
#    return colors

def point_colors_new(cl_array,dic):
    colors = []
    for label in cl_array :
        #print(dic.dic[label][1])
        colors.append(dic.dic[label][1])
    return colors

def fit_func(x, b):
    return x + b

#---------------------------------------------------
#   Read information from the files and set colors #
#---------------------------------------------------

color_dic = color_label()
if type == 'full' :
    color_dic.read_db_model_in_full('db_model.in')
if type == 'class' :
    color_dic.read_db_model_in_class('db_model.in')
color_dic.update_color(map_color)

data_names, data_colors = color_dic.extract_list_color_label()

# Energy
Xe_train, Ye_train, class_e_train = read_x_y_class_new(efile1)
Xe_test, Ye_test, class_e_test = read_x_y_class_new(efile2)
# Force
Xf_train, Yf_train, class_f_train = read_x_y_class_new(ffile1)
Xf_test, Yf_test, class_f_test = read_x_y_class_new(ffile2)
# Stress
if stress_calc == True:
    Xs_train, Ys_train, class_s_train = read_x_y_class_new(sfile1)
    Xs_test, Ys_test, class_s_test = read_x_y_class_new(sfile2)



# --  Set differenet colors for points from different classes --

# Energy
colors_e_train = point_colors_new(class_e_train, color_dic)
print(class_e_test)
print(color_dic)
colors_e_test = point_colors_new(class_e_test, color_dic)
# Force
colors_f_train = point_colors_new(class_f_train, color_dic)
colors_f_test = point_colors_new(class_f_test, color_dic)
# Stress
if stress_calc == True:
    colors_s_train = point_colors_new(class_s_train, color_dic)
    colors_s_test = point_colors_new(class_s_test, color_dic)


#---------------------------------------------------
#                   Plot the data                  #
#---------------------------------------------------

point_size = 25
transparency = 0.6
line_width = 1.0
line_style = 'k--'

if stress_calc == True:
    plt.figure(figsize=(8, 25))
else:
    plt.figure(figsize=(10, 10))

plt.subplots_adjust(top=0.92, bottom=0.05, left=0.1, right=0.98, hspace=0.32,wspace=0.30)
#plt.suptitle('%s' %(plot_name), fontsize=16, fontweight='bold')


# Create x for the linear function: y=x that will appear on the plot
xdata_e_train=(np.linspace(np.amin(Xe_train)-1, np.amax(Xe_train)+1, 10))
xdata_e_test=(np.linspace(np.amin(Xe_test)-1, np.amax(Xe_test)+1, 10))
xdata_f_train=(np.linspace(np.amin(Xf_train)-1, np.amax(Xf_train)+1, 10))
xdata_f_test=(np.linspace(np.amin(Xf_test)-1, np.amax(Xf_test)+1, 10))

if stress_calc == True:
    xdata_s_train=(np.linspace(np.amin(Xs_train)-1, np.amax(Xs_train)+1, 10))
    xdata_s_test=(np.linspace(np.amin(Xs_test)-1, np.amax(Xs_test)+1, 10))

#slope, intercept, r_value, p_value, std_err = stats.linregress(X_train, Y_train)
#print"r-squared:", r_value**2
#print slope, intercept
#plt.plot(xdata, intercept + slope*xdata, 'k', lw=1.0)


# Fit the functions: Y_train = X_train & Y_test = X_test
popt_e_train, pcov_e_train = curve_fit(fit_func, Xe_train, Ye_train, bounds=(0, 0.0001))
popt_e_test, pcov_e_test = curve_fit(fit_func, Xe_test, Ye_test, bounds=(0, 0.0001))
popt_f_train, pcov_f_train = curve_fit(fit_func, Xf_train, Yf_train, bounds=(0, 0.0001))
popt_f_test, pcov_f_test = curve_fit(fit_func, Xf_test, Yf_test, bounds=(0, 0.0001))

if stress_calc == True:
    popt_s_train, pcov_s_train = curve_fit(fit_func, Xs_train, Ys_train, bounds=(0, 0.0001))
    popt_s_test, pcov_s_test = curve_fit(fit_func, Xs_test, Ys_test, bounds=(0, 0.0001))


# Plot subplots

if stress_calc == True:
    pattern = '32'  # 2 plots along x , 3 plots along y
    i1 = 3
    i2 = 2
else:
    pattern = '22'  # 2 plots along x , 2 plots along y
    i1 = 2
    i2 = 2


# ------- E - train ------
#ax = plt.subplot('%s1' %pattern)
ax = plt.subplot(i1, i2, 1)
plt.scatter(Xe_train, Ye_train, s=point_size, c=colors_e_train, alpha=transparency, lw=line_width)
plt.plot(xdata_e_train, fit_func(xdata_e_train, *popt_e_train), line_style, lw=line_width)
#plt.text(0.28*(np.amax(Xe_train)-np.amin(Xe_train)), 0.05*(np.amax(Ye_train)-np.amin(Ye_train)),
#                     'Standard deviation: %.3f' %(np.sqrt(np.diag(pcov_e_train))), fontsize=12)
rmse_e_train = np.sqrt(mean_squared_error(Xe_train, Ye_train))
mae_e_train = mean_absolute_error(Xe_train, Ye_train)
plt.text(0.7,0.1,r'RMSE = %1.3f eV'%(rmse_e_train), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.text(0.7,0.05,r'MAE = %1.3f eV'%(mae_e_train), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

plt.title('Energy train', fontsize=11, fontweight='bold')
plt.xlabel('Computed energy, eV', fontsize=11)
plt.ylabel('DFT-reference, eV', fontsize=11)

# Legend
if len(data_names) > 8 :
    font = 7
else :
    font = 10
handles = [mpatches.Patch(color=colour, label=label) for label, colour in zip(data_names, data_colors)]
plt.legend(handles=handles, loc='best', fontsize=font, shadow=False, scatterpoints=1)


# ------- E - test ------
#ax = plt.subplot('%s2' %pattern)
ax = plt.subplot(i1, i2, 2)
plt.plot(xdata_e_test, fit_func(xdata_e_test, *popt_e_test), line_style, lw=line_width)
plt.scatter(Xe_test, Ye_test, s=point_size, c=colors_e_test, alpha=transparency, lw=line_width)
#plt.text(0.28*(np.amax(Xe_test)-np.amin(Xe_test)), 0.05*(np.amax(Ye_test)-np.amin(Ye_test)),
#                      'Standard deviation: %.3f' %(np.sqrt(np.diag(pcov_e_test))), fontsize=12)
rmse_e_test = np.sqrt(mean_squared_error(Xe_test, Ye_test))
mae_e_test = mean_absolute_error(Xe_test, Ye_test)
plt.text(0.7,0.1,r'RMSE = %1.3f eV'%(rmse_e_test), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.text(0.7,0.05,r'MAE = %1.3f eV'%(mae_e_test), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

plt.title('Energy test', fontsize=11, fontweight='bold')
plt.xlabel('Computed energy, eV', fontsize=11)
plt.ylabel('DFT-reference, eV', fontsize=11)


# ------- F - train ------
#ax = plt.subplot('%s3' %pattern)
ax = plt.subplot(i1, i2, 3)
plt.plot(xdata_f_train, fit_func(xdata_f_train, *popt_f_train), line_style, lw=line_width)
plt.scatter(Xf_train, Yf_train, s=point_size, c=colors_f_train, alpha=transparency, lw=line_width)

rmse_f_train = np.sqrt(mean_squared_error(Xf_train, Yf_train))
mae_f_train = mean_absolute_error(Xf_train, Yf_train)
plt.text(0.7,0.1,r'RMSE = %1.3f eV/$\AA$'%(rmse_f_train), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.text(0.7,0.05,r'MAE = %1.3f eV/$\AA$'%(mae_f_train), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

plt.title('Force train', fontsize=11, fontweight='bold')
plt.xlabel('Computed force, eV/$\AA$', fontsize=11)
plt.ylabel('DFT-reference, eV/$\AA$', fontsize=11)
#plt.text(0.28*(np.amax(Xf_train)-np.amin(Xf_train))+np.amin(Xf_train), 0.02*(np.amax(Yf_train)-np.amin(Yf_train))+np.amin(Yf_train),
#                      'Standard deviation: %.3f' %(np.sqrt(np.diag(pcov_f_train))), fontsize=12)


# ------- F - test ------
#ax = plt.subplot('%s4' %pattern)
ax = plt.subplot(i1, i2, 4)
plt.plot(xdata_f_test, fit_func(xdata_f_test, *popt_f_test), line_style, lw=line_width)
plt.scatter(Xf_test, Yf_test, s=point_size, c=colors_f_test, alpha=transparency, lw=line_width)

rmse_f_test = np.sqrt(mean_squared_error(Xf_test, Yf_test))
mae_f_test = mean_absolute_error(Xf_test, Yf_test)
plt.text(0.7,0.1,r'RMSE = %1.3f eV/$\AA$'%(rmse_f_test), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
plt.text(0.7,0.05,r'MAE = %1.3f eV/$\AA$'%(mae_f_test), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)


plt.title('Force test', fontsize=11, fontweight='bold')
plt.xlabel('Computed force, eV/$\AA$', fontsize=11)
plt.ylabel('DFT-reference, eV/$\AA$', fontsize=11)
#plt.text(0.28*(np.amax(Xf_test)-np.amin(Xf_test))+np.amin(Xf_test), 0.02*(np.amax(Yf_test)-np.amin(Yf_test))+np.amin(Yf_test),
#                  'Standard deviation: %.3f' %(np.sqrt(np.diag(pcov_f_test))), fontsize=12)




if stress_calc == True:
    # ------- S - train ------
    #ax = plt.subplot('%s5' %pattern)
    ax = plt.subplot(i1, i2, 5)
    plt.plot(xdata_s_train, fit_func(xdata_s_train, *popt_s_train), line_style, lw=line_width)
    plt.scatter(Xs_train, Ys_train, s=point_size, c=colors_s_train, alpha=transparency, lw=line_width)

    rmse_s_train = np.sqrt(mean_squared_error(Xs_train, Ys_train))
    mae_s_train = mean_absolute_error(Xs_train, Ys_train)
    plt.text(0.7,0.1,r'RMSE = %1.3f eV/$\AA^3$'%(rmse_s_train), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    plt.text(0.7,0.05,r'MAE = %1.3f eV/$\AA^3$'%(mae_s_train), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    plt.title('Stress train', fontsize=11, fontweight='bold')
    plt.xlabel('Computed stress, eV/$\AA^3$', fontsize=11)
    plt.ylabel('DFT-reference, eV/$\AA^3$', fontsize=11)
#    plt.text(0.28*(np.amax(Xs_train)-np.amin(Xs_train))+np.amin(Xs_train), 0.02*(np.amax(Ys_train)-np.amin(Ys_train))+np.amin(Ys_train),
#                      'Standard deviation: %.3f' %(np.sqrt(np.diag(pcov_s_test))), fontsize=12)


    # ------- S - test ------
    #ax = plt.subplot('%s6' %pattern)
    ax = plt.subplot(i1, i2, 6)
    plt.plot(xdata_s_test, fit_func(xdata_s_test, *popt_s_test), line_style, lw=line_width)
    plt.scatter(Xs_test, Ys_test, s=point_size, c=colors_s_test, alpha=transparency, lw=line_width)

    rmse_s_test = np.sqrt(mean_squared_error(Xs_test, Ys_test))
    mae_s_test = mean_absolute_error(Xs_test, Ys_test)
    plt.text(0.7,0.1,r'RMSE = %1.3f eV/$\AA$'%(rmse_s_test), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
    plt.text(0.7,0.05,r'MAE = %1.3f eV/$\AA$'%(mae_s_test), fontsize = 10, horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)

    plt.title('Stress test', fontsize=11, fontweight='bold')
    plt.xlabel('Computed stress, eV/$\AA^3$', fontsize=11)
    plt.ylabel('DFT-reference, eV/$\AA^3$', fontsize=11)
#    plt.text(0.28*(np.amax(Xs_test)-np.amin(Xs_test))+np.amin(Xs_test), 0.02*(np.amax(Ys_test)-np.amin(Ys_test))+np.amin(Ys_test),
#                  'Standard deviation: %.3f' %(np.sqrt(np.diag(pcov_s_test))), fontsize=12)

plt.tight_layout()
plt.savefig('train_test.png')
plt.show()
