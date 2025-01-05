# python file to generate corrections due to energy of quantum dot due to perturbations
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
k0=9*(10**9)
el0=1.602*(10**(-19))
elconfac=1.0/(5.289*(10**(-19)))
el=el0*elconfac
wavefuncsqnmcon=((1.974*100)**(3))
Jouletoev=1/(1.602*(10**(-19)))

# reading energy from files
nmexcitonEN=path+"\\Wavefunctions\\"+'Energy of Exciton.txt'
fecxitonEN=open(nmexcitonEN,'r')
header_line=next(fecxitonEN) # ignoring the first line of the file
EnergyVals=[[float(num) for num in line.split(' ')] for line in fecxitonEN] #read entire file
fecxitonEN.close()

# as we only careabout the first two energies for the first order energy correction, we will store them
E0=EnergyVals[0][4]
E1=EnergyVals[1][4]

# we only care about the ground state, so we will import only that
lCE=0
lVH=0
# ground state
nmCE=path+"\\Wavefunctions\\"+f'CE_wavefunc_l={lCE}.txt'
fileCE=open(nmCE,'r')
# storing data in array
CE=[[float(num) for num in line.split(' ')] for line in fileCE] #0 : position ; 3 : |\psi|^2 # l=0 for CE
fileCE.close()
nmVH=path+"\\Wavefunctions\\"+f'VH_wavefunc_l={lVH}.txt'
fileVH=open(nmVH,'r')
# storing data in array
VH=[[float(num) for num in line.split(' ')] for line in fileVH] #0 : position ; 3 : |\psi|^2 # l=0 for VH
fileVH.close()
#
# first excited state
nmCE1=path+"\\Wavefunctions\\"+f'CE_wavefunc_l={lCE}.txt'
fileCE1=open(nmCE1,'r')
# storing data in array
CE1=[[float(num) for num in line.split(' ')] for line in fileCE1] #0 : position ; 3 : |\psi|^2 # l=0 for CE
fileCE1.close()
nmVH1=path+"\\Wavefunctions\\"+f'VH_wavefunc_l={lVH+1}.txt'
fileVH1=open(nmVH1,'r')
# storing data in array
VH1=[[float(num) for num in line.split(' ')] for line in fileVH1] #0 : position ; 3 : |\psi|^2 # l=0 for VH
fileVH1.close()

# --------------------------------------------------------------------------------------------- COULOMB CORRECTION --------------------------------------------------------------------------- #
# need to perform two different integrals : |psi_g|^2 H^2 and (|psi_g|^2 H)
# if I want to do a proper integration using simpson, I need to include all the wavefunction generators. 
# Instead we can do a modified simpson, where I take the sum of integrad at two r positions, then take their average.

# first we need to normalise the wavefunction => ground state
#itr=int(1) # counter
#psisqintCE=0.0
#psisqintVH=0.0
#while itr<(len(CE)-1): 
#    r1=CE[itr][0]*confac
    # CE
#    psisq1CE=CE[itr][3]
#    Int1psiCE=4*math.pi*psisq1CE*(r1**2)
    # VH
#    psisq1VH=VH[itr][3]
#    Int1psiVH=4*math.pi*psisq1VH*(r1**2)
    #
#    r2=CE[int(itr+1)][0]*confac
    # CE
#    psisq2CE=CE[int(itr+1)][3]
#    Int2psiCE=4*math.pi*psisq2CE*(r2**2)
    # VH
#    psisq2VH=VH[int(itr+1)][3]
#    Int2psiVH=4*math.pi*psisq2VH*(r2**2)
    #
#    IntpsiavgCE=(Int1psiCE+Int2psiCE)/2.0
#    psisqintCE=psisqintCE+IntpsiavgCE
#    IntpsiavgVH=(Int1psiVH+Int2psiVH)/2.0
#    psisqintVH=psisqintVH+IntpsiavgVH
    #
#    itr=int(itr+2)
#
#NormCE=math.sqrt(1.0/psisqintCE)
#NormVH=math.sqrt(1.0/psisqintVH)
#print(NormCE, NormVH)

# first we need to normalise the wavefunction => first excited state
#itr=int(1) # counter
#psisqintCE1=0.0
#psisqintVH1=0.0
#while itr<(len(CE1)-1): 
#    r1=CE1[itr][0]*confac
    # CE
#    psisq1CE=CE1[itr][3]
#    Int1psiCE1=4*math.pi*psisq1CE*(r1**2)
    # VH
#    psisq1VH=VH1[itr][3]
#    Int1psiVH1=4*math.pi*psisq1VH*(r1**2)
    #
#    r2=CE1[int(itr+1)][0]*confac
    # CE
#    psisq2CE=CE1[int(itr+1)][3]
#    Int2psiCE1=4*math.pi*psisq2CE*(r2**2)
    # VH
#    psisq2VH=VH1[int(itr+1)][3]
#    Int2psiVH1=4*math.pi*psisq2VH*(r2**2)
    #
#    IntpsiavgCE1=(Int1psiCE1+Int2psiCE1)/2.0
#    psisqintCE1=psisqintCE1+IntpsiavgCE1
#    IntpsiavgVH1=(Int1psiVH1+Int2psiVH1)/2.0
#    psisqintVH1=psisqintVH1+IntpsiavgVH1
    #
#    itr=int(itr+2)
#
#NormCE1=math.sqrt(1.0/psisqintCE)
#NormVH1=math.sqrt(1.0/psisqintVH)
#print(NormCE, NormVH)


# ----------------------------------------------------------------------FIRST ORDER ENERGY CORRECTION -------------------------------------------------------#
# performing integration over CE by keeping VH constant at r=0 => ground state
psiVH=VH[0][3]
#print(psiVH)
Int=0.0 # total integration value for CE
delr=CE[1][0]-CE[0][0]
itr=int(1) # counter
while itr<(len(CE)-1): 
    r1=CE[itr][0]
    psisq1=CE[itr][3]
    HC1=-(k0/epsilon)*(el0**2)*(1/r1)
    Int1=4*math.pi*psisq1*HC1*(r1**2)*delr # first point
    #
    r2=CE[int(itr+1)][0]
    psisq2=CE[int(itr+1)][3]
    HC2=-(k0/epsilon)*(el0**2)*(1/r2)
    Int2=4*math.pi*psisq2*HC2*(r2**2)*delr # second point
    #
    Intavg=(Int1+Int2)/2.0
    Int=Int+Intavg
    #
    itr=int(itr+2)
#
print(Int, psiVH)
print(f'\nGround State Energy = {E0} eV')
print(f'First Order Energy Correction Due to Coulomb Interaction = {Int*psiVH*Jouletoev*(10**9)} eV\n')

# performing integration over CE by keeping VH constant at r=0 => first excited state
psiVH1=VH1[0][3]
#print(psiVH)
Int=0.0 # total integration value for CE
itr=int(1) # counter
while itr<(len(CE1)-1): 
    r1=CE1[itr][0]*(10**(-9))
    psisq1=CE1[itr][3]
    HC1=-(k0/epsilon)*(el0**2)*(1/r1)
    Int1=4*math.pi*psisq1*HC1*(r1**2)*delr # first point
    #
    r2=CE1[int(itr+1)][0]*(10**(-9))
    psisq2=CE1[int(itr+1)][3]
    HC2=-(k0/epsilon)*(el0**2)*(1/r2)
    Int2=4*math.pi*psisq2*HC2*(r2**2)*delr # second point
    #
    Intavg=(Int1+Int2)/2.0
    Int=Int+Intavg
    #
    itr=int(itr+2)
#
print(f'\nFirst Excited State Energy = {E1} eV')
print(f'First Order Energy Correction Due to Coulomb Interaction = {Int*psiVH1*Jouletoev*(10**9)} eV\n')

# -------------------------------------------------------------------- FIRST ORDER CORRECTION TO WAVEFUNCTIONS -----------------------------------------------------------------------#

# to find the first order correction to the wavefunction, we need to (technically) know all the states, but we can make an estimate to that.
# we first need to sort the energy vals array with respect to the exciton energy
EnergyValsSort=arraysort(EnergyVals,4)


# -------------------------------------------- we need to find contributiong factors for each energy level ==> for wavefunction correction 
EcorrWave=[0 for i in range(len(EnergyVals))] # storing correction factor
for i in range(1,len(EnergyValsSort)):
    IntState=0.0 # total integration value
    lCE=int(EnergyValsSort[i][0])
    lVH=int(EnergyValsSort[i][1])
    # importing CE wavefunction for that state
    nmce=path+"\\Wavefunctions\\"+f'CE_wavefunc_l={lCE}.txt'
    filece=open(nmce,'r')
    # storing data in array
    CEext=[[float(num) for num in line.split(' ')] for line in filece] #0 : position ; 3 : |\psi|^2
    filece.close()
    # importing VH wavefunction for that state
    nmvh=path+"\\Wavefunctions\\"+f'VH_wavefunc_l={lVH}.txt'
    filevh=open(nmvh,'r')
    # storing data in array
    VHext=[[float(num) for num in line.split(' ')] for line in filevh] #0 : position ; 3 : |\psi|^2
    filevh.close()
    # first we need to normalise the wavefunction
    #itr=int(1) # counter
    #psisqintCEext=0.0
    #psisqintVHext=0.0
    #while itr<(len(CEext)-1): 
    #    r1=CEext[itr][0]*confac
        # CE
    #    psisq1CEext=CEext[itr][3]
    #    Int1psiCEext=4*math.pi*psisq1CEext*(r1**2)
        # VH
    #    psisq1VHext=VHext[itr][3]
    #    Int1psiVHext=4*math.pi*psisq1VHext*(r1**2)
        #
    #    r2=CEext[int(itr+1)][0]*confac
        # CE
    #    psisq2CEext=CEext[int(itr+1)][3]
    #    Int2psiCEext=4*math.pi*psisq2CEext*(r2**2)
        # VH
    #    psisq2VHext=VHext[int(itr+1)][3]
    #    Int2psiVHext=4*math.pi*psisq2VHext*(r2**2)
        #
    #    IntpsiavgCEext=(Int1psiCEext+Int2psiCEext)/2.0
    #    psisqintCEext=psisqintCEext+IntpsiavgCEext
    #    IntpsiavgVHext=(Int1psiVHext+Int2psiVHext)/2.0
    #    psisqintVHext=psisqintVHext+IntpsiavgVHext
        #
    #    itr=int(itr+2)
    #
    #NormCEext=math.sqrt(1.0/psisqintCEext)
    #NormVHext=math.sqrt(1.0/psisqintVHext)
    # storing the r=0 val for psi_VH
    psivh1=VHext[0][1]
    psivhg1=VH[0][1]
    # performing the integration
    itr=int(1) # counter
    while itr<(len(CEext)-1): 
        r1=CEext[itr][0]
        psi1=CEext[itr][1]
        psig1=CE[itr][1]
        HC1=-(k0*(10**27)/epsilon)*(el0**2)*(1/r1) # in nm units
        IntS1=4*math.pi*psi1*psig1*HC1*(r1**2)*delr # first point
        #
        r2=CEext[int(itr+1)][0]
        psi2=CEext[int(itr+1)][1]
        psig2=CE[int(itr+1)][1]
        HC2=-(k0*(10**27)/epsilon)*(el0**2)*(1/r2)
        IntS2=4*math.pi*psi2*psig2*HC2*(r2**2)*delr # second point
        #
        IntStateavg=(IntS1+IntS2)/2.0
        IntState=IntState+IntStateavg
        #
        itr=int(itr+2)
    #
    Ei=EnergyValsSort[i][4]
    IntState=IntState/((E0-Ei)*1.602*(10**(-1))) # converting to nmJoules
    EcorrWave[i]=IntState*psivh1*psivhg1
    #print(IntState)
    #print(f'Correction factor for m={i} (lCE={lCE}, lVH={lVH}) : {IntState*psivh1*psivhg1*(10**9)} eV')
    #print(EcorrWave)
#

# ------------------------------------ generating corrected wavefunctions
nmcorr=path+"\\Wavefunctions\\"+f'FirstCorrectedExciton.txt'
fcorr=open(nmcorr,"w")
fcorr.close()
nmoriginal=path+"\\Wavefunctions\\"+f'UnperturbedExciton.txt'
forig=open(nmoriginal,"w")
forig.close()
# restoring ground state wavefunction
CorrWave=shift_matrix_multiply(CE,VH[0][1],1)
Correction=[0 for i in range(len(CE))]
# running loop to generate corrections
for i in range(1,len(EcorrWave)):
    lCE=int(EnergyValsSort[i][0])
    lVH=int(EnergyValsSort[i][1])
    # importing CE wavefunction for that state
    nmce=path+"\\Wavefunctions\\"+f'CE_wavefunc_l={lCE}.txt'
    filece=open(nmce,'r')
    # storing data in array
    CEext=[[float(num) for num in line.split(' ')] for line in filece] #0 : position ; 3 : |\psi|^2
    filece.close()
    # importing VH wavefunction for that state
    nmvh=path+"\\Wavefunctions\\"+f'VH_wavefunc_l={lVH}.txt'
    filevh=open(nmvh,'r')
    # storing data in array
    VHext=[[float(num) for num in line.split(' ')] for line in filevh] #0 : position ; 3 : |\psi|^2
    filevh.close()
    # first we need to normalise the wavefunction
    #itr=int(1) # counter
    #psisqintCEext=0.0
    #psisqintVHext=0.0
    #while itr<(len(CEext)-1): 
    #    r1=CEext[itr][0]*confac
    #    # CE
    #    psisq1CEext=CEext[itr][3]
    #    Int1psiCEext=4*math.pi*psisq1CEext*(r1**2)
        # VH
    #    psisq1VHext=VHext[itr][3]
    #    Int1psiVHext=4*math.pi*psisq1VHext*(r1**2)
        #
    #    r2=CEext[int(itr+1)][0]*confac
        # CE
    #    psisq2CEext=CEext[int(itr+1)][3]
    #    Int2psiCEext=4*math.pi*psisq2CEext*(r2**2)
        # VH
    #    psisq2VHext=VHext[int(itr+1)][3]
    #    Int2psiVHext=4*math.pi*psisq2VHext*(r2**2)
        #
    #    IntpsiavgCEext=(Int1psiCEext+Int2psiCEext)/2.0
    #    psisqintCEext=psisqintCEext+IntpsiavgCEext
    #    IntpsiavgVHext=(Int1psiVHext+Int2psiVHext)/2.0
    #    psisqintVHext=psisqintVHext+IntpsiavgVHext
        #
    #    itr=int(itr+2)
    #
    #NormCEext=math.sqrt(1.0/psisqintCEext)
    #NormVHext=math.sqrt(1.0/psisqintVHext)
    # generating corrected wavefunctions
    if EcorrWave[i]!=0.0:
        for itr in range(len(CEext)):
            Correction[itr]=Correction[itr]+EcorrWave[i]*(CEext[itr][1])*(VHext[0][1])
        #
    #
#
# corrected wavefunctions
for i in range(len(CorrWave)):
    append_file(nmoriginal,f'{CorrWave[i][0]} {CorrWave[i][1]} {0.0} {modulus(CorrWave[i][1])}\n')
    CorrWave[i][1]=CorrWave[i][1]+Correction[i]
    append_file(nmcorr,f'{CorrWave[i][0]} {CorrWave[i][1]} {0.0} {modulus(CorrWave[i][1])}\n')
# ------------------------------------------------ #
end=time.time()
print(f'Time taken for code to run = {end-begin}\n')
