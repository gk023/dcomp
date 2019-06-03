#!/user/bin/env python
#
# file: chisq_test.py 
#
# Creating a function that runs dcomp, gets the data, and checks the 
# chi squared 
#
# Usage: python run_dcomp.py <dcomp card file> <cross_section data> <datatype>
#
# Data type can be 1 for elastic, 2 for reaction, and 3 for total cross sections
#
# Data is expected to be: Energy <space> Data <space> Error, but modifications
# could be added to handle data delimited in other ways. For now, error is 
# assumed to be relative, but more options could be added 
#
###################################################################

from subprocess import Popen, PIPE, STDOUT
import numpy as np
from scipy.interpolate import CubicSpline
import sys
import matplotlib.pyplot as plt

##################### function definitions ########################
#Takes in a parameter list from a card file and runs dcomp
#Returns numpy array: [energy, elastic xs, rxn xs, and total xs]
def dcomp(parameters):
        
	my_command='../dcomp_main.x'
	print(my_command)
	p=Popen([my_command,'-u'], stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1)
	print p.stdout.readline()

	for i in range(len(parameters)):
		print >>p.stdin, parameters[i]
		p.stdin.flush()
		print p.stdout.readline()

	with p.stdout:
		for line in iter(p.stdout.readline, b''):
			print line,
	p.wait()

	cross_sections=[]

	f=open('cross_sections.out','r')
	for line in f.readlines()[1:]:
		line=line.strip().split()
		line=[float(x) for x in line]
		cross_sections.append(line)
        f.close()
	cross_sections=np.asarray(cross_sections)

	return(cross_sections)

#################################################################


#################### Main Code ###################################

parameter_list = []
cardfile = sys.argv[1]
datafile = sys.argv[2] 
datatype = int(sys.argv[3])

#Read in parameters from properly formatted card file
f=open(cardfile,'r')

for line in f.readlines():
	param=line.strip().split()[-1] #Grabs the last thing in the line
	parameter_list.append(param)

f.close()


th_data = dcomp(parameter_list)

th_data = th_data[:,[0, datatype]]  #slice to get the right data type

th_interp = CubicSpline(th_data[:,0],th_data[:,1])

exp_data=[]

f=open(datafile,'r')

for line in f.readlines():
	line=line.strip().split()
	
	
	line=[float(x) for x in line]

	line[2]=line[2]*line[1]

	exp_data.append(line)

f.close()

exp_data=np.asarray(exp_data)

th_prediction = np.asarray( [th_interp(x) for x in exp_data[:,0]] )

#print(th_prediction)

#compute the chi square
chisq = np.sum( ((th_prediction - exp_data[:,1])/exp_data[:,2])**2 )

print("\n Chi Square/d.o.f:\n"+str(chisq/len(exp_data)))

fig=plt.figure(figsize=(12,8))

plt.rc('text',usetex=True)

plt.plot(exp_data[:,0],th_prediction,lw=2,color='darkblue')
plt.errorbar(exp_data[:,0],exp_data[:,1],yerr=exp_data[:,2],color='k',fmt='o',markersize=20,elinewidth=2)

plt.xlabel(r"$E$ (MeV)",fontsize=50)
plt.ylabel(r"$\sigma$ (mb)",fontsize=50)

plt.tick_params(axis='both',labelsize=50)

plt.tight_layout()
plt.show()

