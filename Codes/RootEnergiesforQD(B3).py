# Generating difference of wavefunctions of Quantum Dot and finding what is the root energy for CE or VH for various l values (Boundary 3 - infinite well)

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
import matplotlib
from matplotlib import pyplot
 

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

# DERIVATIVE
def derivative(f, x, h): # arguments = equation, point at which derivative is to be found, h
    y=(f(x+h)-f(x-h))/(2*h)
    return(y)
#

# NEWTON RAPHSON
def newton_raphson(f, a, max_itr, nm): #arguments ==> formula, trial root, maximum iterations, name of file to store errors
    err=pow(10, -15)
    n=0
    X=[0 for i in range(max_itr)]
    while n<max_itr: #so that it does not cross max number of allowed iterations
        if n==0:
            X[n]=a
        else:
            h=X[n-1]*0.0001
            X[n]=X[n-1]-(f(X[n-1])/derivative(f, X[n-1], h))
            #error = x_{n+1}-x_{n} o equivalently, x_{n}-x_{n-1}
            append_file(nm, f'{n} {abs(X[n] - X[n - 1])}\n') #appends absolute error value with iteration number to file
            if abs(X[n]-X[n-1])<err: return(X[n], True)
            #
        n+=1
        #
    return(X[n-1], False) #if root not obtained even after max_iterations
#

def simpson(f, lower, upper, N, h, name):
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


############################################# ACTUAL CODE STARTS ######################################################

#clock beginning time to measure time taken by code to run
begin=time.time()

#setting path of file
path=os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)
# -------------------------------------------------------------------- #
#defining all constants => masses (natural units, electron volts)
hbar=1.0
m0=0.511*(10**6) #rest mass of electon in natural units
mecore=0.11*m0
meshell=0.14*m0
mhcore=0.8*m0
mhshell=0.51*m0
# ----------------------------------------------------------------------- #


#defining all constants => radii (in nanometers)
rcorenm=1.5
shellthicknessnm=1.4
rshellnm=rcorenm+shellthicknessnm
# converting to natural units
confac=5.066*(10**(-3))
rcore=confac*rcorenm
rshell=confac*rshellnm
# ---------------------------------------------------------------------- #

#defining "k" as a function of energy, potential and mass (for e/h in c/v bands respectively)
def k(E,V,M):
    hbar=1.0
    if (E-V)>=0 : kval=math.sqrt((2*(E-V)*M)/(hbar**2))
    else : kval=cmath.sqrt((2*(E-V)*M)/(hbar**2))
    #print(f'E-V={E-V}\n k={kval}')
    return(kval)
#
# ---------------------------------------------------------------------- #

#defining potential for conduction (for e) band as a function of r (in eV)
def Vc(r,core,shell):
    #defining all constants => potentials (natural units, electron volts)
    Vcshell=0.0 #wrt Vcshell=0
    Vccore=-0.1
    #note, Vccore and Vvcore are not zero wrt each other
    #
    if r<=core: V=Vccore
    elif r>core and r<=shell: V=Vcshell
    else: V=10**15 #tends to infinity
    return(V)
#

#defining potential for valence (for h) band as a function of r (in eV)
def Vv(r,core,shell):
    #defining all constants => potentials (natural units, electron volts)
    Vvcore=-0.1-1.65 #wrt Vcshell=0
    Vvshell=Vvcore-0.65 
    #note, Vccore and Vvcore are not zero wrt each other
    #
    if r<=core: V=Vvcore
    elif r>core and r<=shell: V=Vvshell
    else: V=-10**15 #tends to infinity
    return(V)
#
# --------------------------------- functions for wavefunction (CE and VH)------------------------------------------- #

# defining wavefunction of core for electrons (in conduction band)
def psi_CE_core(Acecore,l,kcecore,r):
    psicecore=Acecore*spherical_jn(l, kcecore*r)
    return(psicecore)
#

# defining wavefunction of shell for electrons (in conduction band)
def psi_CE_shell(Aceshell,Bceshell,l,kceshell,r):
    if isinstance(kceshell,complex)==True : psiceshell=Aceshell*spherical_in(l,(kceshell*r).imag)+Bceshell*spherical_kn(l,(kceshell*r).imag)
    else : psiceshell=Aceshell*spherical_jn(l,kceshell*r)+Bceshell*spherical_yn(l,kceshell*r)
    return(psiceshell)
#

#defining wavefunction of core for holes (in valence band)
def psi_VH_core(Avhcore,l,kvhcore,r):
    psivhcore=Avhcore*spherical_jn(l,kvhcore*r)
    return(psivhcore)
#

# defining wavefunction for shell for holes (in valence band)
def psi_VH_shell(Avhshell,Bvhshell,l,kvhshell,r):
    psivhshell=Avhshell*spherical_jn(l,kvhshell*r)+Bvhshell*spherical_yn(l,kvhshell*r)
    return(psivhshell)
#
# ----------------------------------------- A,B functions (CE and VH) -------------------------------------------------- #

# defining A_CE_shell
def A_CE_shell(Acecore,l,kcecore,kceshell,rcore,Mcecore,Mceshell):
    if isinstance(kceshell, complex)==True:
        numceshell=(spherical_jn(l,kcecore*rcore).item()*spherical_kn(l,(kceshell*rcore).imag,derivative=True).item()-((Mceshell*kcecore)/(Mcecore*(kceshell.imag)))*spherical_jn(l,kcecore*rcore,derivative=True).item()*spherical_kn(l,(kceshell*rcore).imag).item())
        denomceshell=spherical_in(l,(kceshell*rcore).imag).item()*spherical_kn(l,(kceshell*rcore).imag,derivative=True).item()-spherical_in(l,(kceshell*rcore).imag,derivative=True).item()*spherical_kn(l,(kceshell*rcore).imag).item()
    #
    else :
        numceshell=(spherical_jn(l,kcecore*rcore).item()*spherical_yn(l,kceshell*rcore,derivative=True).item()-((Mceshell*kcecore)/(Mcecore*kceshell))*spherical_jn(l,kcecore*rcore,derivative=True).item()*spherical_yn(l,kceshell*rcore).item())
        denomceshell=spherical_jn(l,kceshell*rcore).item()*spherical_yn(l,kceshell*rcore,derivative=True).item()-spherical_jn(l,kceshell*rcore,derivative=True).item()*spherical_yn(l,kceshell*rcore).item()
    #
    #print((kceshell*rcore).imag,spherical_in(l,(kceshell*rcore).imag).item(),spherical_kn(l,(kceshell*rcore).imag),spherical_in(l,kceshell*rcore,derivative=True),spherical_kn(l,kceshell*rcore,derivative=True).tolist(),denomceshell)
    Aceshell=Acecore*(numceshell/denomceshell)
    return(Aceshell)
#

# defining B_CE_shell
def B_CE_shell(Acecore,l,kcecore,kceshell,rcore,Mcecore,Mceshell):
    if isinstance(kceshell,complex)==True:
        Kceshellrcore=(kceshell*rcore).imag
        numceshell=(spherical_jn(l,kcecore*rcore).item()*spherical_in(l,Kceshellrcore,derivative=True).item()-((Mceshell*kcecore)/(Mcecore*(kceshell.imag)))*spherical_jn(l,kcecore*rcore,derivative=True).item()*spherical_in(l,Kceshellrcore).item())
        denomceshell=spherical_in(l,Kceshellrcore,derivative=True).item()*spherical_kn(l,Kceshellrcore).item()-spherical_in(l,Kceshellrcore).item()*spherical_kn(l,Kceshellrcore,derivative=True).item()
    #
    else :
        numceshell=(spherical_jn(l,kcecore*rcore)*spherical_jn(l,kceshell*rcore,derivative=True)-((Mceshell*kcecore)/(Mcecore*kceshell))*spherical_jn(l,kcecore*rcore,derivative=True)*spherical_jn(l,kceshell*rcore))
        denomceshell=spherical_jn(l,kceshell*rcore,derivative=True)*spherical_yn(l,kceshell*rcore)-spherical_jn(l,kceshell*rcore)*spherical_yn(l,kceshell*rcore,derivative=True)
    #
    Bceshell=Acecore*(numceshell/denomceshell)
    return(Bceshell)
#

# defining A_VH_shell
def A_VH_shell(Avhcore,l,kvhcore,kvhshell,rcore,Mvhcore,Mvhshell):
    numvhshell=(spherical_jn(l,kvhcore*rcore)*spherical_yn(l,kvhshell*rcore,derivative=True)-((Mvhshell*kvhcore)/(Mvhcore*kvhshell))*spherical_jn(l,kvhcore*rcore,derivative=True)*spherical_yn(l,kvhshell*rcore))
    denomvhshell=spherical_jn(l,kvhshell*rcore)*spherical_yn(l,kvhshell*rcore,derivative=True)-spherical_jn(l,kvhshell*rcore,derivative=True)*spherical_yn(l,kvhshell*rcore)
    Avhshell=Avhcore*(numvhshell/denomvhshell)
    return(Avhshell)
#

# defining B_CE_shell
def B_VH_shell(Avhcore,l,kvhcore,kvhshell,rcore,Mvhcore,Mvhshell):
    numvhshell=(spherical_jn(l,kvhcore*rcore)*spherical_jn(l,kvhshell*rcore,derivative=True)-((Mvhshell*kvhcore)/(Mvhcore*kvhshell))*spherical_jn(l,kvhcore*rcore,derivative=True)*spherical_jn(l,kvhshell*rcore))
    denomvhshell=spherical_jn(l,kvhshell*rcore,derivative=True)*spherical_yn(l,kvhshell*rcore)-spherical_jn(l,kvhshell*rcore)*spherical_yn(l,kvhshell*rcore,derivative=True)
    Bvhshell=Avhcore*(numvhshell/denomvhshell)
    return(Bvhshell)
#

#After this, all functions will be in a for loops. Outermost loops will be "l" and "m" which define the orders. 
#Then, the loops will be for "r" and "E". This is the final generation of data section.


# -------------------------------------------------- CONDUCTION ELECTRONS -------------------------------------------------- #
# defining big file to store roots of CE wavefunctions
nmrootsCE=path+"\\Wavefunctions\\"+f'CE_Roots(B3).txt'
frootsCE=open(nmrootsCE, "w")
frootsCE.close()
# defining function to call for generating roots
def CE_rootsfunc(lwave):
    # defining constants and logistics about terms to generate (step sizes, initial points, final points, number of terms)
    hE=0.001 #step size in energy (in eV)
    E0=-0.1
    NE=100
    # first need to define a simple function of psidiff as E to easily call Newton Rhapson later
    def psiCEShellasE(E):
        Vceshell=Vc(rshell,rcore,rshell) #potential for the iteration
        Vcecore=Vc(rcore,rcore,rshell)
        kcecoreO=k(E,Vcecore,mecore)
        kceshellO=k(E,Vceshell,meshell)
        AceshellO=A_CE_shell(1,lwave,kcecoreO,kceshellO,rcore,mecore,meshell)
        BceshellO=B_CE_shell(1,lwave,kcecoreO,kceshellO,rcore,mecore,meshell)
        psiceshell0=psi_CE_shell(AceshellO,BceshellO,lwave,kceshellO,rshell)
        return(psiceshell0)
    #
    # naming files ==> may not use always in the code
    nmdiffCE=path+"\\Wavefunctions\\"+"\\EnergyVals\\"+f'CE_l(B3)={lwave}.txt' 
    fdiffCE=open(nmdiffCE,"w")
    fdiffCE.close()
    #now calling the above function => we will only do this to generate plots (rest of the time, this will be turned off)
    for iCE in range(NE):
        E=E0+(iCE**(3/2))*hE #defines E for that iteration
        psival=psiCEShellasE(E)
        append_file(nmdiffCE,f'{E} {psival.real} {psival.imag}\n')
    #
    # from plots, we define a guess root, then call Newton Rhapson
    maxiterCE=1000
    trialrootCE=float(input(f'Please enter value of trial root for l={lwave} (CE).\n'))
    nmCEerr=path+"\\Wavefunctions\\"+"\\EnergyVals\\"+f'CE_err_l={lwave}.txt'
    fCEerr=open(nmCEerr,"w")
    fCEerr.close()
    #as the function outputs a value and a bool function, the first term is the root (if found), the second is the bool
    ErootCE,ErootCEbool=newton_raphson(psiCEShellasE,trialrootCE,maxiterCE,nmCEerr)
    if ErootCEbool==True: 
        #print(f'\nRoot found for CE, and the value is given as = {ErootCE}\n')
        append_file(nmrootsCE, f'{lwave} {ErootCE}\n')
    #
    else: print(f'\nRoot not found.\n')    
    return(0)
#

# ------------------------------------------------ VALENCE HOLES ----------------------------------------------- #
# defining big file to store roots of CE wavefunctions
nmrootsVH=path+"\\Wavefunctions\\"+f'VH_Roots(B3).txt'
frootsVH=open(nmrootsVH, "w")
frootsVH.close()
def VH_rootsfunc(lwave):
    # defining constants and logistics about terms to generate (step sizes, initial points, final points, number of terms)
    hE=10**(-5) #step size in energy (in eV)
    E0=-0.1
    NE=100
    # first need to define a simple function of psidiff as E to easily call Newton Rhapson later
    def psiVHShellasE(E):
        Vvhshell=Vv(rshell,rcore,rshell) #potential for the iteration
        Vvhcore=Vv(rcore,rcore,rshell)
        kvhcoreO=k(E,Vvhcore,mhcore)
        kvhshellO=k(E,Vvhshell,mhshell)
        AvhshellO=A_VH_shell(1,lwave,kvhcoreO,kvhshellO,rcore,mhcore,mhshell)
        BvhshellO=B_VH_shell(1,lwave,kvhcoreO,kvhshellO,rcore,mhcore,mhshell)
        psivhshell0=psi_VH_shell(AvhshellO,BvhshellO,lwave,kvhshellO,rshell)
        return(psivhshell0)
    #
    # naming files ==> may not use always in the code
    nmdiffVH=path+"\\Wavefunctions\\"+"\\EnergyVals\\"+f'VH_l(B3)={lwave}.txt' 
    fdiffVH=open(nmdiffVH,"w")
    fdiffVH.close()
    #now calling the above function => we will only do this to generate plots (rest of the time, this will be turned off)
    for iVH in range(NE):
        E=E0+(iVH**(3/2))*hE #defines E for that iteration
        append_file(nmdiffVH,f'{E} {psiVHShellasE(E)}\n')
    #
    # from plots, we define a guess root, then call Newton Rhapson
    maxiterVH=1000
    trialrootVH=float(input(f'Please enter value of trial root for l={lwave} (VH).\n'))
    nmVHerr=path+"\\Wavefunctions\\"+"\\EnergyVals\\"+f'VH_err_l={lwave}.txt'
    fVHerr=open(nmVHerr,"w")
    fVHerr.close()
    #as the function outputs a value and a bool function, the first term is the root (if found), the second is the bool
    ErootVH,ErootVHbool=newton_raphson(psiVHShellasE,trialrootVH,maxiterVH,nmVHerr)
    if ErootVHbool==True: 
        #print(f'\nRoot found for CE, and the value is given as = {ErootCE}\n')
        append_file(nmrootsVH, f'{lwave} {ErootVH}\n')
    #
    else: print(f'\nRoot not found.\n')
    return(0)
#

# now, we can generate roots for any of the two ==> CE of VH
lmax=4
for lw in range(lmax):
    ceroots=CE_rootsfunc(lw)
    #vhroots=VH_rootsfunc(lw)
#
    
# ------------------------------------------------ #
end=time.time()
print(f'Time taken for code to run = {end-begin}\n')
