from pylab import * 
import pandas as pd 
from pylab import * 
import numpy as np

first_number = input('What numerical extension does your first file/trigger have?\n')
last_number =  input('What numberical extension does your last file/trigger have?\n')
bin_number = input('What bin number would you like to use for your histogram?\n')
log_scale = input('Would you like a log scale for your histogram? Type True or False\n')
normalized = input('Would you like normalized histogram? Type True or False\n')
minimum = input('Would you like to plot the minimum in every trigger? Type True or False\n')

for f in range(first_number,last_number,1):
	
	globals()['ts%d' % f] = np.array([])
	globals()['ss%d' % f] = np.array([])
	 
for i in range(first_number,last_number,1):
	
	if(0<i<10):
		#print("did you make it")
		globals()['fns%d' % i] = pd.read_csv("C1--Table--0000"+str(i)+".csv",  usecols=[3,4],names=['colA', 'colB'],header=None)
		globals()['col1s%d' % i] = globals()['fns%d' % i]["colA"].values
		globals()['col2s%d' % i] = globals()['fns%d' % i]["colB"].values
		globals()['ts%d' % i] = np.append(globals()['col1s%d' % i],globals()['ts%d' %i])
		globals()['ss%d' % i] = np.append(globals()['col2s%d' % i],globals()['ss%d' %i])
	
	if(10<=i<=99):
		print('did you come here')
		globals()['fns%d' % i] = pd.read_csv("C1--Table--000"+str(i)+".csv",  usecols=[3,4],names=['colA', 'colB'],header=None)
		globals()['col1s%d' % i] = globals()['fns%d' % i]["colA"].values
		globals()['col2s%d' % i] = globals()['fns%d' % i]["colB"].values
		globals()['ts%d' % i] = np.append(globals()['col1s%d' % i],globals()['ts%d' %i])
		globals()['ss%d' % i] = np.append(globals()['col2s%d' % i],globals()['ss%d' %i])
	
arrays = np.array([])
for j in range(first_number,last_number,1):
	if minimum == True:
		if(0<j<10):
			arrays = np.append(globals()['ss%d' %j].min(),arrays)
		if(10<=j<=99):
			arrays = np.append(globals()['ss%d' %j].min(),arrays)
	if minimum == False:
		if(0<j<10):
			arrays = np.append(globals()['ss%d' %j],arrays)
		if(10<=j<=99):
			arrays = np.append(globals()['ss%d' %j],arrays)		
	print(globals()['ss%d' % j],j)
print(len(arrays))

figure(0)
hist(arrays*1000,bins=bin_number,normed=normalized,color='r',log=log_scale)
xlabel('min amplitude [mV]')
grid('on')
show()
	

