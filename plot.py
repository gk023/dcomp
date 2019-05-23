# a plotting code for DCOMP cross sections
#
#############################################
import matplotlib.pyplot as plt
import sys

print("Reaction (1) or total (2) cross section?")
response = int(raw_input(">>"))

infile=sys.argv[1]

E=[]
xs1=[]
xs2=[]

f=open(infile,'r')


for line in f.readlines()[1:]:
	line = line.strip().split()
        E.append(float(line[0]))
        xs1.append(float(line[response+1]))

f.close()

#plt.rc('text',usetex=True)

fig=plt.figure(figsize=(12,8))


plt.plot(E,xs1,lw=3,ls='-',marker='o',ms=15,color='darkblue',label='Cross Sections')
plt.xlabel('$E$ (MeV)',fontsize=50)
plt.ylabel('$\sigma$  (mb)', fontsize=50)
plt.tick_params(axis='both',labelsize=50)
plt.xlim([0,200])
plt.legend(loc='best',fontsize='30')

plt.tight_layout()
plt.savefig("cross_sections.png")

plt.show()
