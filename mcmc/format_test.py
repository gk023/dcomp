#!/user/bin/env python
#
# file: format_test.py
#
# Makes a toy data file with a dcomp run
#
# Usage: format_test.py <data> <output> <type>
#
# Data type: 1. Elastic, 2. Reaction, 3. Total
#
############################################################

import sys
import numpy as np

############################################################

infile = sys.argv[1]
outfile = sys.argv[2]
data_type=int(sys.argv[3])

noise_level=0.1

file_data=[]

f=open(infile,'r')
for line in f.readlines():
	line=line.strip().split()
	file_data.append(line)

f.close()

file_data = np.asarray(file_data)

f=open(outfile,'w')

for i in range(1,len(file_data)):
	mean=float(file_data[i][data_type])
        noisy=mean + np.random.normal(0.,noise_level*mean)
	f.write(file_data[i][0]+'   '+str(noisy)+'   '+str(noise_level)+'\n')

f.close()
