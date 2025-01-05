# energy correction due to stark effect

import math
import cmath
import sys
import time
import os
import numpy as np
import scipy
from scipy.special import spherical_jn
from scipy.special import spherical_yn
from scipy.special import spherical_in
from scipy.special import spherical_kn

# all external functions needed for the code

# READ FILE
def read_matrix(x): #more than one column #parameter: x ==> name of file
    f=open(x,'r') #'r' ==> read only
    X=[[float(num) for num in line.split(' ')] for line in f]
    f.close()
    return(X)
#

# APPEND FILE
def append_file(name, str): #arguments: name ==> name of file, str ==> string to append
    f=open(name, 'a') #'a' ==> append file
    f.write(str)
    f.close()
    #
#

#MOD SQUARE OF A COMPLEX NUMBER
def modulus(x): #x ==> complex number
   real=x.real**2
   imag=x.imag**2
   modsq=real+imag 
   return(modsq)
#

#MOD SQUARE OF A COMPLEX NUMBER
def modulus(x): #x ==> complex number
   real=x.real**2
   imag=x.imag**2
   modsq=real+imag 
   return(modsq)
#

# Integration
def simpson(f, lower, upper, h, name):
    N=int(abs(upper-lower)/h)
    I=0
    X=[0 for i in range(N+1)] #for storing N+1 values, including upper and lower
    X[0]=lower
    for i in range(1,N+1):
        X[i]=lower+h*i #finding the next point
        mid=(X[i-1]+X[i])/2 #finding mid value
        h1=(X[i]-X[i-1])/2  #diving every slice in two equal slices, and finding h relevant to that slice
        #calculating the integration
        I=I+(h1/3)*(f(X[i-1])+4*f(mid)+f(X[i]))
        #
    #storing N vs I value in file
    append_file(name, f'{N} {I}\n')
    return(I)
#

# sorting array with respect to certain column ==> bubble algortihm
def arraysort(X,col):
    n=len(X)
    ncol=len(X[0])
    for i in range(n):
        for j in range(i+1,n):
            if (X[i][col]>X[j][col]): 
                for k in range(ncol) :
                    temp=X[i][k]
                    X[i][k]=X[j][k]
                    X[j][k]=temp
                #
            #
        #
    #
    return(X)
#

# SHIFT VALUES OF MATRIX BY MULTIPLCATION ==> multiply some constant/variable to terms of a matrix
def shift_matrix_multiply(M,val,col): #arguments: M ==> original matrix which needs shifting, val ==> value to add to the matrix (all terms), col ==> which column of the matrix to multiply (if 1D matrix, then col=0)
    n=len(M) #number of rows
    if col==0: 
        CopyM=[0 for i in range(n)]
        for i in range(n): 
            CopyM[i]=M[i]*val
    #
    else: 
        m=len(M[0])
        CopyM=[[0 for j in range(m)] for i in range(n)]
        for i in range(n): 
            for j in range(m):
                if j ==col : CopyM[i][col]=M[i][col]*val
                else : CopyM[i][j]=M[i][j]
            #
        #
    #
    return(CopyM)
#            
    
#----------------------------------------------------- Actual Code Starts -------------------------------------------------#

#clock beginning time to measure time taken by code to run
begin=time.time()

#setting path of file
path=os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)
# -------------------------------------------------------------------- #

# defining constants in natural units
confac=5.066*(10**(-3))
epsilon0=1
epsilon=5.367
k=(1/(4*math.pi*epsilon0*epsilon))
el0=1.602*(10**(-19))
elconfac=1.0/(5.289*(10**(-19)))
el=el0*elconfac
evJoule=1.602*(10**(-19))

EF2=1014223.237


# reading energy from files
nmexcitonEN=path+"\\Wavefunctions\\"+'Energy of Exciton.txt'
fecxitonEN=open(nmexcitonEN,'r')
header_line=next(fecxitonEN) # ignoring the first line of the file
EnergyVals=[[float(num) for num in line.split(' ')] for line in fecxitonEN] #read entire file
fecxitonEN.close()

# as we only careabout the first two energies for the first order energy correction, we will store them
E0=EnergyVals[0][4]
E1=EnergyVals[1][4]
Ecorr0=-0.03438 #first order correction to ground state

# first importing the ground state
# ground state
nmExc=path+"\\Wavefunctions\\"+f'FirstCorrectedExciton.txt'
fileExc=open(nmExc,'r')
# storing data in array
Exc=[[float(num) for num in line.split(' ')] for line in fileExc] #0 : position ; 3 : |\psi|^2 # l=0 for CE
fileExc.close()

# performing the two integrations
IntH1=0.0
IntH2=0.0
itr=int(1) # counter
delr=Exc[1][0]-Exc[0][0]
while itr<(len(Exc)-1): 
    r1=Exc[itr][0]
    psisq1=Exc[itr][3]
    HC2r1=r1**2
    # first point
    Int2r1=4*math.pi*psisq1*HC2r1*(r1**2)*(1/3)*delr
    #
    r2=Exc[int(itr+1)][0]
    psisq2=Exc[int(itr+1)][3]
    HC2r2=r2**2
    # second point
    Int2r2=4*math.pi*psisq2*HC2r2*(r2**2)*(1/3)*delr
    #
    IntavgH2=(Int2r1+Int2r2)/2.0
    IntH2=IntH2+IntavgH2
    #
    itr=int(itr+2)
#

EcorrStark=(-(el0**2))*2*(IntH2*(10**(-18)))/((E0+Ecorr0-E1)*evJoule)
print(f'Correction={EcorrStark}')