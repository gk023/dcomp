#!/usr/bin/env python
#
# file: run_dcomp.py
#
# This python script runs the DCOMP code. 
# Reads in parameters from card file and passes them
# onto dcomp_main
#
# Usage: python run_dcomp.py <card file>
#
########################################################

from subprocess import Popen, PIPE, STDOUT
import sys
import matplotlib.pyplot as plt

parameter_list = []

cardfile = sys.argv[1]

# Read in parameters from properly formatted cardfile
f=open(cardfile,'r')

for line in f.readlines():
	param=line.strip().split()[-1]
	parameter_list.append(param)

f.close()

my_command = './dcomp_main.x'
p = Popen([my_command, '-u'], stdout=PIPE, stdin=PIPE, stderr=STDOUT, bufsize=1)
print p.stdout.readline()

for i in range(len(parameter_list)):
	print >>p.stdin, parameter_list[i]
        p.stdin.flush()
        print p.stdout.readline()

with p.stdout:
	for line in iter(p.stdout.readline, b''):
		print line,
p.wait()

