# python file for generating exciton wavefunction

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


#----------------------------------------------------- Actual Code Starts -------------------------------------------------#

#clock beginning time to measure time taken by code to run
begin=time.time()

#setting path of file
path=os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)
# -------------------------------------------------------------------- #

# creating big file to srtore value of exciton energies E = Eg + Ec + Eh
Eg=1.65
nmexcitonEN=path+"\\Wavefunctions\\"+'Energy of Exciton.txt'
fecxitonEN=open(nmexcitonEN,'w')
fecxitonEN.close()
append_file(nmexcitonEN,"l_{CE} l_{VH} E_{exciton} E_{CE} E_{VH}"+f' (Eg={Eg})\n')

# importing energy for conduction electrons
nmrootsCE=path+"\\Wavefunctions\\"+f'CE_Roots(B3)_final.txt'
frootsCE=open(nmrootsCE,'r')
RootsCE=[[float(num) for num in line.split(' ')] for line in frootsCE] #read entire file

nmrootsVH=path+"\\Wavefunctions\\"+f'VH_Roots(B3)_final.txt'
frootsVH=open(nmrootsVH,'r')
RootsVH=[[float(num) for num in line.split(' ')] for line in frootsVH] #read entire file


# for loop to generate all values at once
for itrCE in range(4):
    for itrVH in range(4):
        # defining l values for wavefunction (integer)
        lCE=itrCE
        lVH=itrVH
        #
        #finding and storing exciton energy
        ENexciton=(RootsCE[lCE][1]+0.1)+(RootsVH[lVH][1]+0.65)
        append_file(nmexcitonEN,f'{lCE} {lVH} {RootsCE[lCE][1]+0.1} {RootsVH[lVH][1]+0.65} {ENexciton}\n')
        #
        # reading data from appropriate wavefunction files
        nmfileCE=path+"\\Wavefunctions\\"+f'CE_wavefunc_l={lCE}.txt'
        nmfileVH=path+"\\Wavefunctions\\"+f'VH_wavefunc_l={lVH}.txt'
        fileCE=open(nmfileCE,'r')
        fileVH=open(nmfileVH,'r')
        #
        # storing in big arrays (both files should have same values of r)
        CEWave=[[float(num) for num in line.split(' ')] for line in fileCE] #read entire file
        VHWave=[[float(num) for num in line.split(' ')] for line in fileVH] #read entire file
        #
        # creating new file to store exciton wavefunction
        nmexciton=path+"\\Wavefunctions\\"+f'Exciton(lCE={lCE},lVH={lVH}).txt'
        fileexciton=open(nmexciton,'w')
        fileexciton.close()
        #
        # creating new array of appropriate length
        nCE=len(CEWave)
        nVH=len(VHWave)
        if nCE==nVH :
            n=nCE
            Exciton=[[0 for j in range(4)] for i in range(n)] # creating an empty n cross 4 matrix
            for i in range(n):
                Exciton[i][0]=CEWave[i][0] # radial value (should be same for both files)
                compexciton=complex(CEWave[i][1],CEWave[i][2])*complex(VHWave[i][1],VHWave[i][2]) # gives total complex wavefunction
                # storing the real and the imaginary terms
                Exciton[i][1]=compexciton.real
                Exciton[i][2]=compexciton.imag 
                # storing modulus value
                Exciton[i][3]=modulus(compexciton)
                # storing in file
                append_file(nmexciton, f'{Exciton[i][0]} {Exciton[i][1]} {Exciton[i][2]} {Exciton[i][3]}\n')
            #
        #
        else : print(f'The number of rows should be the same in both files.')
    #
#

# ------------------------------------------------ #
end=time.time()
print(f'Time taken for code to run = {end-begin}\n')
