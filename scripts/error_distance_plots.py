#!/usr/bin/env python
#VERSION CHANGES
# Changed to keep all st dev within the max interpolation distance
# Changed so that it derives 1 equation for all data.
# Took out half and 1 percent files as backup

import numpy as np
from numpy import copy
import sys
from scipy import optimize
import random


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText


#Get data (.csv) from source_int.sh
data=sys.argv[1]
name=sys.argv[2]
ss_samp_den_tmp=float(sys.argv[3])
max_int_dist_tmp=(sys.argv[4])

short_name="All Terrain"
ss_samp_den=("%.3f" % ss_samp_den_tmp)
max_int_dist=int(max_int_dist_tmp)
print("Max int dist is", max_int_dist)

try:
	print("starting data import into array")
	my_data = np.loadtxt(data, delimiter=' ')
	print("created array from data with no memory issues!")
except MemoryError:
	print("Memory Error")

error=my_data[:,0]
distance=my_data[:,1]
x=distance
y=error

print("max distance of data is", np.max(distance))
#print "maxout,cov,infodict,mesg,ier=optimize.leastsq interpolation distance is", int(max_int_dist)

nbins = 10
#nbins=range(0,int(max_int_dist),10)
#nbins='auto'

n, _ = np.histogram(x, bins=nbins)
sy, _ = np.histogram(x, bins=nbins, weights=y)
sy2, _ = np.histogram(x, bins=nbins, weights=y*y)

with np.errstate(invalid='ignore'):
	mean = sy / n
	std = np.sqrt(sy2/n - mean*mean)
#print "x bin distance is", (_[1:] + _[:-1])/2
#print "st. dev is", std
#print "sample size is", n

fig = plt.figure()
ax = plt.subplot(111)

plt_data=ax.scatter(x,y, zorder=1, label='Measurements', marker=".", color="black", s=20)
plt.show
plt_data_uncert=ax.errorbar((_[1:] + _[:-1])/2, mean, yerr=std, fmt='r-', linewidth=3)


box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

#plt.axes.get_xaxis().set_ticks([])
plt.tick_params(
axis='x',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='on',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelbottom='on') # labels along the bottom edge are off
#handles,labels = ax.get_legend_handles_labels()
#handles = [handles[0], handles[2], handles[1]]
#labels = [labels[0], labels[2], labels[1]]

anchored_text = AnchoredText(short_name, loc=2)
anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(anchored_text)

str_ss_samp_den="Sampling Density = " + str(ss_samp_den) + " %"
anchored_text3 = AnchoredText(str_ss_samp_den, loc=4)
anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(anchored_text3)

#ax.text(right, top, 'Deep Bathy', horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, backgroundcolor='grey', color='black')

plt.legend([plt_data, plt_data_uncert], ['Interpolation Error', 'Mean +/- St. Deviation'], loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True, ncol=2, fontsize=14)
#plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=2, mode="expand", borderaxespad=0.)
#plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], loc=8)
#plt.legend([plt_data, plt_data_uncert, plt_st_dev, plt_mean, plt_st_error], ['Measurement', 'Measurement Uncertainty', 'Standard Deviation', 'Mean Elevation', 'Standard Error'], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.legend(loc='upper right')
plt.xlabel('Distance from Measurement (cells)', fontsize=14)
plt.ylabel('Interpolation Error (m)', fontsize=14)
plt.xlim(xmin=0)
plt.xlim(xmax=int(max_int_dist)+1)


#fig,ax = plt.subplots()
# Plot the data
#data_line = ax.plot(x,y, label='Data', marker='o')
# Plot the average line
#mean_line = ax.plot(x,y_mean, label='Mean', linestyle='--'
#fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
#ax.plot(x_val, elevs)
png_name = str(name)+"_scatter.png"
plt.savefig(png_name)   # save the figure to file
plt.close()



print("center of x bin distance is", (_[1:] + _[:-1])/2)
print("st. dev is", std)
print("sample size is", n)

#always keep the first three data points after 0,0
std_final=[0]
std_final.append(std[0])
std_final.append(std[1])
std_final.append(std[2])



#Get first 3 measurements
#then only add if greater


#Here is where I originally only kept st devs if they were bigger than previous distance
#changed to keep 1 smaller std and then break
# i=3
# while i < len(std):
# 	if std[i] >= std[i-1]:
# 		#print "appending", std[i]
# 		std_final.append(std[i])
# 	else:
# 		std_final.append(std[i])
# 		print "exiting"
# 		break
# 	i = i+1

# print "std_final used for deriving equations is", std_final


#Changed to keep all st devs
i=3
while i < len(std):
	std_final.append(std[i])
	i = i+1

bins_orig=(_[1:] + _[:-1])/2
bins_final=[0]

for j in bins_orig[0:len(std_final)-1]:
	#print "appending", j
	bins_final.append(j)
	

print("LENGTH BINS", len(bins_final))
print("LENGTH STD", len(std_final))

fig = plt.figure()
ax = plt.subplot(111)

plt_data=ax.scatter(bins_final,std_final, zorder=1, label='Error St. Dev.', marker="o", color="black", s=30)
plt.show
#bins_final=bins_orig[0:len(std_final)]

print("plotting data")
print("distance is", bins_final)
print("st dev is", std_final)

xdata=bins_final
ydata=std_final


fitfunc = lambda p, x: p[0] + p[1] * x ** (p[2])
errfunc = lambda p, x, y: (y - fitfunc(p, x))

#out,success = optimize.leastsq(errfunc, [0,0.5,0.5],args=(xdata, ydata),maxfev=3000)




#need 3 points to fit polynomial?
if len(xdata)==2:
	print("ONLY 2 POINTS -- Adding 3rd of 0,0 to make polyfit work")
	xdata=[0]+xdata
	ydata=[0]+ydata


coeff_guess=[0,0.5,0.5]

print("xdata is", xdata)
print("ydata is", ydata)
print("errfunc is", errfunc)
print("coeff_guess is", coeff_guess)

out,cov,infodict,mesg,ier=optimize.leastsq(errfunc,coeff_guess,args=(xdata,ydata),full_output=True)

ssErr = (infodict['fvec']**2).sum()
ssTot = ((y-y.mean())**2).sum()
rsquared = 1-(ssErr/ssTot )


#change power coefficint to be non-zero (which results in constant uncertainty)
if out[2]<0.001:
	out[2]=0.001

print("%g + %g*x^%g"%(out[0],out[1],out[2]))
print("coefficients are", out)

coeff1=float("{0:.3f}".format(out[1]))
coeff2=float("{0:.3f}".format(out[2]))
print("coeff1 is", coeff1)
print("coeff2 is", coeff2)

st_dev_fit=[]
bins_fit=range(0,int(max_int_dist)+2)
for i in bins_fit:
	#st_dev_fit.append(out[0] + out[1]*i**out[2])
	st_dev_fit.append(out[1]*i**out[2])

#print "log_coeff is", log_coeff
#print "intercept is", intercept
print("bins_fit are", bins_fit)
print("st_dev_fit are", st_dev_fit)

plt_best_fit,=ax.plot(bins_fit,st_dev_fit, zorder=1, linewidth=2.0)
plt.show

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])

#plt.axes.get_xaxis().set_ticks([])
plt.tick_params(
axis='x',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='on',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelbottom='on') # labels along the bottom edge are off
#handles,labels = ax.get_legend_handles_labels()
#handles = [handles[0], handles[2], handles[1]]
#labels = [labels[0], labels[2], labels[1]]

anchored_text = AnchoredText(short_name, loc=2)
anchored_text.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(anchored_text)

#no r2 value
anchored_text2 = AnchoredText(" $y = {%gx}^{%g}$ "%(coeff1,coeff2), loc=1)
#add r2 value using below
#anchored_text2 = AnchoredText(" $y = {%gx}^{%g}$      $r^2=%g$ "%(coeff1,coeff2,rsquared), loc=1)
anchored_text2.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(anchored_text2)

anchored_text3 = AnchoredText(str_ss_samp_den, loc=4)
anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
ax.add_artist(anchored_text3)


#anchored_text3 = AnchoredText("r2 = %g"%(rsquared), loc=1)
#anchored_text3.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
#ax.add_artist(anchored_text3)

#ax.text(right, top, 'Deep Bathy', horizontalalignment='right', verticalalignment='top', transform=ax.transAxes, backgroundcolor='grey', color='black')

plt.legend([plt_data, plt_best_fit], ['Interpolation Error St Dev', 'Best-Fit Line'], loc='upper center', bbox_to_anchor=(0.5, 1.15), fancybox=True, shadow=True, ncol=2, fontsize=14)



#plt.axes.get_xaxis().set_ticks([])
plt.tick_params(
axis='x',          # changes apply to the x-axis
which='both',      # both major and minor ticks are affected
bottom='on',      # ticks along the bottom edge are off
top='off',         # ticks along the top edge are off
labelbottom='on') # labels along the bottom edge are off
plt.xlabel('Distance from Measurement (cells)', fontsize=14)
plt.ylabel('Interpolation Error St Dev (m)', fontsize=14)
plt.xlim(xmin=0)
plt.xlim(xmax=int(max_int_dist)+1)
plt.ylim(ymin=0)
y_max=max(std_final)+(0.25*max(std_final))
print("y_max is", y_max)
plt.ylim(ymax=y_max)



png_name = str(name)+"_best_fit.png"
plt.savefig(png_name)   # save the figure to file
plt.close()


#print "rsquared value is", rsquared
#print (0, coeff1, coeff2)

sys.exit([0, coeff1, coeff2])






