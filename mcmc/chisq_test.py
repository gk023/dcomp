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
		print p.stdout.readlines()

	with p.stdout:
		for line in iter(p.stdout.readline,b''):
			print line.rstrip()
	p.wait()

	cross_sections=[]

	f=open('cross_sections.out','r')
	for line in f.readlines()[:1]:
		line=line.strip().split()
		cross_sections.append(line)
        f.close()
	cross_sections=np.asarray(cross_sections)

	return(cross_sections)

#################################################################


#################### Main Code ###################################

parameter_list = []
cardfile = sys.argv[1]
datafile = sys.argv[2] 
datatype = sys.argv[3]

#Read in parameters from properly formatted card file
f=open(cardfile,'r')

for line in f.readlines():
	param=line.strip().split()[-1] #Grabs the last thing in the line
	parameter_list.append(param)

f.close()

print(parameter_list)

data = dcomp(parameter_list)

print(data)
