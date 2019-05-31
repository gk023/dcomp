##!/user/bin/env python
#
# file: mcmc_run_dcomp.py 
#
# Creating a function that runs dcomp using MCMC to get parameter
# distributions and stores the cross sections for the runs
#
# Usage: python run_dcomp.py <search card> <datatype>
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
from scipy.stats import multivariate_normal
import sys

##################### function definitions ########################

#Takes in a parameter list from a card file and runs dcomp
#Requires numpy array of enegies
#Returns numpy array: [energy, elastic xs, rxn xs, and total xs]
def dcomp(parameters,energies):
        

	my_command='/home/gk023/dcomp/dcomp_main.x'
	energies = energies.tolist()
	energies = [str(x) for x in energies]
	number_of_energies = str(len(energies))
	parameters.append(number_of_energies)
        inputs = parameters + energies
        #print("Running DCOMP with the following inputs:\n")
        #print(inputs,'\n')
	p=Popen([my_command,'-u'], stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1)
	print p.stdout.readline()

	for i in range(len(inputs)):
		print >>p.stdin, inputs[i]
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

#Takes in np arrays of the dcomp predictions and experimental data to get a chi square
#per degree of freedom. Theory values are what is returned by the dcomp function. Experimental
#data are read in from a file in the main code below. data_type picks the right kind of data
#in the array of dcomp predictions

def chisq(theory,exp,data_type):

	#trim the stuff you don't need from dcomp
	theory = theory[:,[0,data_type]]

	#Evaluate the chi square
	chisq_per_dof = np.sum( ((theory[:,1] - exp[:,1])/exp[:,2])**2 )/len(exp_data)

	print("\n Chi Square/d.o.f:\n"+str(chisq_per_dof))
	
	return(theory,chisq_per_dof)


#The likelihood function

def likelihood(theory,exp,data_type):

	chi_square = chisq(theory,exp,data_type)[1]

	return( np.exp(-0.5*chi_square), chi_square )

def prior(center,sigma,sample):

	return( multivariate_normal.pdf(x=sample,mean=center,cov=np.diag(sigma)) )

def take_step(parameters,step_size):

	parameters = np.asarray(parameters)
	step = np.random.random_sample((len(parameters),))
	step_scaled = 2.*step_size*parameters*(step-1./2.)

	return(step_scaled)

def update_parameters(parameters,step):
	parameters = [float(x) for x in parameters]
	parameters = np.asarray(parameters)
        parameters[4:6] = parameters[4:6] + step[0:2]
	parameters[9:14] = parameters[9:14] + step[2:]
	parameters = parameters.tolist()
	parameters = [str(x) for x in parameters]
	return(parameters)

#################################################################
numerics=[]

parameter_list = [] #list of dcomp parameters
exp_data=[] #list for the experimental data

saved=[]
xs_predictions=[]
naccepted=0
nrejected=0

cardfile = sys.argv[1]
datatype = int(sys.argv[2])


#Read your parameters in from properly formatted card file
f=open(cardfile,'r')

isParameters = True
isNumerics = False

for line in f.readlines():
	
	line=line.strip().split()
	
	if(isParameters == True and len(line) > 0 and line[0] != '&Data'):
		parameter_list.append(line[-1])

        if(len(line) > 0 and line[0] == '&Data'):
		isParameters = False

	if(isParameters == False and isNumerics == False and len(line) > 0 and line[0] != '&Data' and line[0] != '&Numerics' ):
		
		line=[float(x) for x in line]
		line[2]=line[2]*line[1]

		exp_data.append(line)

	if(len(line) > 0 and line[0] == '&Numerics'):
		isNumerics = True

	if(len(line) >0 and isNumerics == True and line[0] != '&Numerics'):
		numerics.append(line[-1])

f.close()



exp_data=np.asarray(exp_data)

#print(exp_data[:,0])

#Set up the numerics for the run
max_accept = int(numerics[0])
skip = float(numerics[1])
burn_in = int(numerics[2])
step_scale = float(numerics[3])

print("Pulling "+str(max_accept)+" parameters with:\nThinning: "+str(skip)+"\nBurn-in: "+str(burn_in)+"\nStep Scaling: "+str(step_scale))

#Grab the parameters we're varying in the run
start=parameter_list[4:6]+parameter_list[9:14]
start=[float(x) for x in start] #make them floats
centers = np.asarray(start) #The means for the priors

#Get data, likelihood, prior*likelihood combo for inital data set
start_th = dcomp(parameter_list,exp_data[:,0])

start_like = likelihood(start_th,exp_data,datatype)[0]

start_priorlike = start_like*prior(centers,centers,centers)


count=0


while(len(saved) < max_accept):

	count +=1 

	#Take a random step
	step = take_step(start,step_scale)
	new_list = update_parameters(parameter_list,step)[:20]

	#Repeat process to get prior*likelihood for new set
	new = new_list[4:6]+new_list[9:14]
	new = [float(x) for x in new]
	new = np.asarray(new)
	
        #print(new_list)		

	new_th = dcomp(new_list,exp_data[:,0])
	new_like, new_chisquare = likelihood(new_th,exp_data,datatype)
	new_priorlike = new_like*prior(centers,centers,new)

	ratio = new_priorlike/start_priorlike #ratio for Metroplois Hastings
        
    	#Use the Metropolis-Hastings Algorithm
	if (ratio >= np.random.uniform()):
		parameter_list = new_list
		start = new
		start_priorlike = new_priorlike
		
		#Only take sets after the burn in one every 'skip' acceptances
		if((count > burn_in)):
			naccepted += 1
		
			if(naccepted % skip == 0):	
				saved.append(new.tolist()+[new_chisquare,new_like])
                        	xs_predictions.append(new_th[:,[0,datatype]].tolist())
		

	elif(count > burn_in):
		nrejected += 1

	if(count == burn_in):
		print("BURN IN PERIOD ENDED")
	if(count >= burn_in and len(saved) != max_accept):
		print(str(len(saved))+" parameters accepted out of "+str(max_accept))
		print("Run number: " + str(count+1))

print("\n\n\n% Accepted: "+str(100.*(float(naccepted)/float(count-burn_in))))
print("% Rejected: "+str(100.*(float(nrejected)/float(count-burn_in))))

#print(saved)
	

g=open("parameters.txt","w")

for line in saved:
	for i in range(len(line)):
		g.write(str(line[i])+"   ")
	g.write("\n")
g.close()	


g=open("xs_bands.txt","w")

for line in xs_predictions:
	for i in range(len(line)):
		for j in range(len(line[i])):
			g.write(str(line[i][j])+"   ")
	g.write("\n&\n")
g.close()	
	



