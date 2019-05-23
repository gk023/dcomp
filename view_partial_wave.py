#Code to plot wave functions
#Run: python view_partial_wave.py <path/to/wavefunction/output>

import matplotlib.pyplot as plt
import sys

r=[]
Re=[]
Im=[]

f=open(sys.argv[1],'r')

for line in f.readlines()[1:]:
	line=line.strip().split()
	r.append(float(line[0]))
	Re.append(float(line[1]))
	Im.append(float(line[2]))

f.close()

plt.rc('text',usetex=True)

fig=plt.figure(figsize=(12,8))

plt.plot(r,Re,color="darkblue",label="Re[$\chi_0(R)$]",ls='-',lw=3)
plt.plot(r,Im,color="orangered",label="Im[$\chi_0(R)$]",ls='--',lw=3)
plt.xlabel('$R$ (arb units)',fontsize=50)
plt.xlabel('$\chi_O(R)$ (arb units)',fontsize=50)
plt.tick_params(axis='both',which='major',labelsize=45)
plt.legend(loc='best',fontsize=45)

plt.show()
#plt.savefig('wave.png') #Uncomment if you want to save the wavefucntion
