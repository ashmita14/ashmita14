#Python file to do statistical analysis on data collected for electrostatic energy of quantum dot due to N number of disks

import math
import sys
import numpy as np
import os
import matplotlib.pyplot as plt



#all functions being used in this code

#X represents any one-dimensional array passed to it
def Avg(X):
    sumX=0.0
    n=len(X)
    for i in range(n):
        sumX=sumX+X[i]
    #
    avgX=sumX/n
    return avgX
#

#defining a function which calculates standard deviation of dataset (as an array) passed to it
#X represents any one-dimensional array passed to it
def StDeviation(X,avg):
    n=len(X)
    sumXstd=0.0
    for i in range(n):
        sumXstd+=((X[i]-avg)**2)
    #
    StdX=math.sqrt(sumXstd/n)
    return StdX
#

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

# SHIFT VALUES OF MATRIX BY ADDITION ==> add some constant/variable to all terms of a matrix
def shift_matrix_add(M,val,col): #arguments: M ==> original matrix which needs shifting, val ==> value to add to the matrix (all terms), col ==> which column of the matrix to shift (if 1D matrix, then col=0)
    n=len(M) #number of rows
    if col==0: 
        CopyM=[0 for i in range(n)]
        for i in range(n): 
            CopyM[i]=M[i]+val
    #
    else: 
        m=len(M[0])
        CopyM=[[0 for j in range(m)] for i in range(n)]
        for i in range(n): CopyM[i][col]=M[i][col]+val
    #
    return(CopyM)
#

# RAISE DATA TO HIGHER TERMS TO SOME POWER ==> raise a column of a matrix to higher power
def power_matrixcol(M,pow,col): #arguments: M ==> original matrix, pow ==> power to which a column of that matrix is to be raised, col ==> the column that needs to be changed
    n=len(M) #number of rows
    if col==0: 
        CopyM=[0 for i in range(n)]
        for i in range(n): 
            term=M[i]
            CopyM[i]=term**pow
        #
    #
    else:
        m=len(M[0])
        CopyM=[[0 for j in range(m)] for i in range(n)] 
        for i in range(n): 
            term=M[i][col]
            CopyM[i][col]=term**pow
        #
    #
    return(CopyM)
#
    


#################################### ACTUAL CODE STARTS ##############################################

path=os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)

#We need to generate arrays by importing the file of interest and reading the data into appropriate arrays. Gor now, we are interested in the energy and "t"=>curve parameter
N=3000
nmenergy=path+'\\Results\\'+f'N={N},AllTrajData.txt' #information of the disks
energyarray=read_matrix(nmenergy)
numrows=len(energyarray) #number of rows
numcol=len(energyarray[0]) #number of columns
#defining empty arrays
Energy=[0 for i in range(numrows)]
#filling the arrays
for i in range(numrows):
    Energy[i]=float(energyarray[i][2]) # averages data in "j+1"th column, where energyarray[i][j]
#
#finding average and standard deviation
avgU=Avg(Energy)
StdU=StDeviation(Energy,avgU)
print(f'Average Ez for all trajectories is = {avgU}')
print(f'Standard Deviation of Ez for all trajectories is = {StdU}')
#
#shifiting the data by n*Standard Deviation
Eshift=shift_matrix_add(Energy,0,0) #col=0, as U is a 1D array of E_z; Ushift is E_z+n*standard deviation
#
#then squaring that column of the shifted matrix to get |E_z|^2
Eshiftsquare=power_matrixcol(Eshift,2,0) #this is |E_z|^2
#
#finding average and standard deviation of |E_z|^2
avgUSS=Avg(Eshiftsquare)
StdUSS=StDeviation(Eshiftsquare,avgUSS)
print(f'Obtained Average |Ez^2| for all trajectories is = {avgUSS}')
print(f'Obtained Standard Deviation of |Ez^2| for all trajectories is = {StdUSS}')
#
# If E_z is gaussian, then we have an expected avg and std for E_z^2
k=1
lmd=(avgU/StdU)**2
avgUSSexp=(StdU**2)*(k+lmd)
StdUSSexp=(StdU**2)*math.sqrt(2*(k+2*lmd))
print(f'Expected Average |Ez^2| for all trajectories is = {avgUSSexp}')
print(f'Expected Standard Deviation of |Ez^2| for all trajectories is = {StdUSSexp}')
#
#plotting histogram
#plt.hist(Energy,bins=300)
#plt.show()
