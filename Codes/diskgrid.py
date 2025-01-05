#python code to generate disks 

#first we will import some basic libraries
import math
import os
import time
import sys
import numpy as np

path=os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)

epsl0=8.85*(10**(-18)) #in Farad micrometer-inverse units
#
N=5 #defined length ==> number of disks
#
X=np.random.uniform(0,100,N) #chosen from a uniform distribution between 0 and 100 (X-coordinate)
Y=np.random.uniform(0,100,N) #chosen from a uniform distribution between 0 and 100 (Y-coordinate)
R=abs(np.random.normal(3,1,N)) #chosen from a Gaussian distribution with mean=3,std=1 (radius)
V0=np.random.normal(0,10**(-3),N) #choose from gaussian distribution centred at 0
#
Q=[16*epsl0*R[i]*V0[i] for i in range(N)] #defining Q's for the disks
#

#Creating File to store information about disks
nmdisk=path+'\\Results\\'+f'DiskGrid(N={N}).txt' #information of the disks
fd1=open(nmdisk,"w")
fd1.close()
fd2=open(nmdisk,"a")
#

#storing information about the disk (x,y,radi,Q,V0)
for i in range(N): fd2.write(f'{X[i]} {Y[i]} {0.0} {R[i]} {V0[i]} {Q[i]}\n')