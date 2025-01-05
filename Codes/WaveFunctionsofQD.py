# Generating wavefunctions of Quantum Dot

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

# DERIVATIVE
def derivative(f, x, h): # arguments = equation, point at which derivative is to be found, h
    y=(f(x+h)-f(x-h))/(2*h)
    return(y)
#

# NEWTON RAPHSON
def newton_raphson(f, a, max_itr, nm): #arguments ==> formula, trial root, maximum iterations, name of file to store errors
    err=pow(10, -9)
    n=0
    X=[0 for i in range(max_itr)]
    while n<max_itr: #so that it does not cross max number of allowed iterations
        if n==0:
            X[n]=a
        else:
            h=X[n-1]*0.001
            X[n]=X[n-1]-(f(X[n-1])/derivative(f, X[n-1], h))
            #error = x_{n+1}-x_{n} o equivalently, x_{n}-x_{n-1}
            append_file(nm, f'{n} {abs(X[n] - X[n - 1])}\n') #appends absolute error value with iteration number to file
            if abs(X[n]-X[n-1])<err: return(X[n], True)
            #
        n+=1
        #
    return(X[n-1], False) #if root not obtained even after max_iterations
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
wavefuncsqnmcon=1/((1.974*100)**(3))
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
    Vcshell=0.0 
    Vccore=-0.1 #wrt Vcshell=0
    #note, Vccore and Vvcore are not zero wrt each other
    #
    if abs(r)<=core: V=Vccore
    elif abs(r)>core and abs(r)<=shell: V=Vcshell
    else: V=10**15 #tends to infinity
    return(V)
#

#defining potential for valence (for h) band as a function of r (in eV)
def Vv(r,core,shell):
    #defining all constants => potentials (natural units, electron volts)
    Vvcore=-0.65 #wrt Vcshell=0
    Vvshell=0.0 
    #note, Vccore and Vvcore are not zero wrt each other
    #
    if abs(r)<=core: V=Vvcore
    elif abs(r)>core and abs(r)<=shell: V=Vvshell
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
    psiceshell=Aceshell*spherical_jn(l,kceshell*r)+Bceshell*spherical_yn(l,kceshell*r)
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
    numceshell=(spherical_jn(l,kcecore*rcore).item()*spherical_yn(l,kceshell*rcore,derivative=True).item()-((Mceshell*kcecore)/(Mcecore*kceshell))*spherical_jn(l,kcecore*rcore,derivative=True).item()*spherical_yn(l,kceshell*rcore).item())
    denomceshell=spherical_jn(l,kceshell*rcore).item()*spherical_yn(l,kceshell*rcore,derivative=True).item()-spherical_jn(l,kceshell*rcore,derivative=True).item()*spherical_yn(l,kceshell*rcore).item()
    #
    #print((kceshell*rcore).imag,spherical_in(l,(kceshell*rcore).imag).item(),spherical_kn(l,(kceshell*rcore).imag),spherical_in(l,kceshell*rcore,derivative=True),spherical_kn(l,kceshell*rcore,derivative=True).tolist(),denomceshell)
    Aceshell=Acecore*(numceshell/denomceshell)
    return(Aceshell)
#

# defining B_CE_shell
def B_CE_shell(Acecore,l,kcecore,kceshell,rcore,Mcecore,Mceshell):
    numceshell=(spherical_jn(l,kcecore*rcore).item()*spherical_jn(l,kceshell*rcore,derivative=True).item()-((Mceshell*kcecore)/(Mcecore*kceshell))*spherical_jn(l,kcecore*rcore,derivative=True).item()*spherical_jn(l,kceshell*rcore).item())
    denomceshell=spherical_jn(l,kceshell*rcore,derivative=True).item()*spherical_yn(l,kceshell*rcore).item()-spherical_jn(l,kceshell*rcore).item()*spherical_yn(l,kceshell*rcore,derivative=True).item()
    Bceshell=Acecore*(numceshell/denomceshell)
    return(Bceshell)
#

# defining A_VH_shell
def A_VH_shell(Avhcore,l,kvhcore,kvhshell,rcore,Mvhcore,Mvhshell):
    numvhshell=(spherical_jn(l,kvhcore*rcore).item()*spherical_yn(l,kvhshell*rcore,derivative=True).item()-((Mvhshell*kvhcore)/(Mvhcore*kvhshell))*spherical_jn(l,kvhcore*rcore,derivative=True).item()*spherical_yn(l,kvhshell*rcore).item())
    denomvhshell=spherical_jn(l,kvhshell*rcore).item()*spherical_yn(l,kvhshell*rcore,derivative=True).item()-spherical_jn(l,kvhshell*rcore,derivative=True).item()*spherical_yn(l,kvhshell*rcore).item()
    Avhshell=Avhcore*(numvhshell/denomvhshell)
    return(Avhshell)
#

# defining B_CE_shell
def B_VH_shell(Avhcore,l,kvhcore,kvhshell,rcore,Mvhcore,Mvhshell):
    numvhshell=(spherical_jn(l,kvhcore*rcore).item()*spherical_jn(l,kvhshell*rcore,derivative=True).item()-((Mvhshell*kvhcore)/(Mvhcore*kvhshell))*spherical_jn(l,kvhcore*rcore,derivative=True)*spherical_jn(l,kvhshell*rcore).item())
    denomvhshell=spherical_jn(l,kvhshell*rcore,derivative=True).item()*spherical_yn(l,kvhshell*rcore).item()-spherical_jn(l,kvhshell*rcore).item()*spherical_yn(l,kvhshell*rcore,derivative=True).item()
    Bvhshell=Avhcore*(numvhshell/denomvhshell)
    return(Bvhshell)
#
# ---------------------------------------- defining constants and building up wavefunctions ---------------------------------------- #

#defining logistics about terms to generate (step sizes, initial points, final points, number of terms)
hr=10**(-3)*confac #step size in r (in converted units)
hE=10**(-7) #step size in energy (in eV)
r0=0.0
E0=10**(-10)
rf=confac*3.0
Ef=0.0002
Nr=int(math.ceil((rf-r0)/hr))
NE=int(math.ceil((Ef-E0)/hE))
#------------------------------------------------------------------------------#

# We need to import the value of l and correponding roots from the file.

# ----------------------------------------- CONDUCTION ELECTRONS --------------------------------------------- #
nmrootsCE=path+"\\Wavefunctions\\"+f'CE_Roots(B3)_final.txt'
frootsCE=open(nmrootsCE,'r')
RootsCE=[[float(num) for num in line.split(' ')] for line in frootsCE] #read entire file

#After this, all functions will be in a for loops. Outermost loops will be "l" and "m" which define the orders. 
#Then, the loops will be for "r" and "E". This is the final generation of data section.

for itr in range(len(RootsCE)):
    lwave=RootsCE[itr][0]
    ErootCE=round(RootsCE[itr][1],4)
    print(lwave, ErootCE)
    #
    #Next, we need to normalise the complete wavefunction to find Acore, which will allow us to obtain all Ace and Bce, thus psice (first for conduction electrons)
    # for this, we need a function which squares psicoreCE*r^2 for given root energy value
    #def psicoreCEmodsquare(rce):
    #    if abs(rce)>rcore : return(0)
    #    else:
    #        # defining potential and K for iteration
    #        Vce=Vc(rce,rcore,rshell) #potential for the iteration
    #        kcecoreMod=k(ErootCE,Vce,mecore) # k for iteration
    #        psiCErsquare=(rce**2)*modulus(psi_CE_core(1,lwave,kcecoreMod,rce))*4*math.pi
    #        #print(Vce,kcecoreMod,psiCErsquare)
    #        return(psiCErsquare)
        #
    #
    def psiCEmodsquare(rce):
        if abs(rce)<=rcore:
            # defining potential and K for iteration
            Vce=Vc(rce,rcore,rshell) #potential for the iteration
            kcecoreMod=k(ErootCE,Vce,mecore) # k for iteration
            psiCErcoresquare=(rce**2)*modulus(psi_CE_core(1,lwave,kcecoreMod,rce))*4*math.pi
            #print(Vce,kcecoreMod,psiCErsquare)
            return(psiCErcoresquare)
        #
        elif abs(rce)>rcore and abs(rce)<=rshell:
            VceshellOO=Vc(rce,rcore,rshell) #potential for the iteration
            VcecoreOO=Vc(rcore,rcore,rshell) #constant potential of core
            kcecoreOO=k(ErootCE,VcecoreOO,mecore)
            kceshellOO=k(ErootCE,VceshellOO,meshell)
            AceshellOO=A_CE_shell(1,lwave,kcecoreOO,kceshellOO,rcore,mecore,meshell)
            BceshellOO=B_CE_shell(1,lwave,kcecoreOO,kceshellOO,rcore,mecore,meshell)
            psiCErshellsquare=(rce**2)*modulus(psi_CE_shell(AceshellOO,BceshellOO,lwave,kceshellOO,rce))*4*math.pi
            return(psiCErshellsquare)
        #
        else: return(0)
    #
    #now calling the numerical integration function
    nmpsiCEnorm=path+"\\Wavefunctions\\"+f'Normalisation of psiCEcore.txt'
    fpsiCEnorm=open(nmpsiCEnorm,"w")
    fpsiCEnorm.close()
    IntCEcore=simpson(psiCEmodsquare,r0,rshell,hr,nmpsiCEnorm)  
    #print(IntCEcore)  
    AcoreCE=math.sqrt(1/IntCEcore)
    print(f'Normalization constant AcoreCE = {AcoreCE}\n')
    #
    #now generating the complete wavefunction for given lwave, CE for a given energy value =>ErootCE
    # we will also generate |psi|^2
    nmCEwave=path+"\\Wavefunctions\\"+f'CE_wavefunc_l={int(lwave)}.txt'
    fCEwave=open(nmCEwave,"w")
    fCEwave.close()
    for rvals in range(Nr):
        r=r0+rvals*hr
        if abs(r)<=rcore:
            VcecoreOO=Vc(r,rcore,rshell) #potential for the iteration
            kcecoreOO=k(ErootCE,VcecoreOO,mecore)
            psiCErcore=psi_CE_core(AcoreCE,lwave,kcecoreOO,r)
            append_file(nmCEwave,f'{r/confac} {psiCErcore.real*(wavefuncsqnmcon**(1/2))} {psiCErcore.imag*(wavefuncsqnmcon**(1/2))} {modulus(psiCErcore)*(wavefuncsqnmcon)}\n')
        #
        elif abs(r)>rcore and abs(r)<=rshell:
            VceshellOO=Vc(r,rcore,rshell) #potential for the iteration
            VcecoreOO=Vc(rcore,rcore,rshell) #constant potential of core
            kcecoreOO=k(ErootCE,VcecoreOO,mecore)
            kceshellOO=k(ErootCE,VceshellOO,meshell)
            AceshellOO=A_CE_shell(AcoreCE,lwave,kcecoreOO,kceshellOO,rcore,mecore,meshell)
            BceshellOO=B_CE_shell(AcoreCE,lwave,kcecoreOO,kceshellOO,rcore,mecore,meshell)
            psiCErshell=psi_CE_shell(AceshellOO,BceshellOO,lwave,kceshellOO,r)
            append_file(nmCEwave,f'{r/confac} {psiCErshell.real*(wavefuncsqnmcon**(1/2))} {psiCErshell.imag*(wavefuncsqnmcon**(1/2))} {modulus(psiCErshell)*(wavefuncsqnmcon)}\n')
        #
        else: append_file(nmCEwave,f'{r/confac} {0} {0} {0}\n')
    #
#

# -------------------------------------------------------------------------- VALENCE HOLES --------------------------------------------------------------- #

nmrootsVH=path+"\\Wavefunctions\\"+f'VH_Roots(B3)_final.txt'
frootsVH=open(nmrootsVH,'r')
RootsVH=[[float(num) for num in line.split(' ')] for line in frootsVH] #read entire file

#After this, all functions will be in a for loops. Outermost loops will be "l" and "m" which define the orders. 
#Then, the loops will be for "r" and "E". This is the final generation of data section.

for itr in range(len(RootsVH)):
    lwave=RootsVH[itr][0]
    ErootVH=round(RootsVH[itr][1],4)
    print(lwave, ErootVH)
    #
    #Next, we need to normalise the core wavefunction to find Acore, which will allow us to obtain all Ace and Bce, thus psice (first for conduction electrons)
    # for this, we need a function which squares psicoreCE*r^2 for given root energy value
    def psicoreVHmodsquare(rvh):
        if abs(rvh)<=rcore:
            # defining potential and K for iteration
            Vvh=Vv(rvh,rcore,rshell) #potential for the iteration
            kvhcoreMod=k(ErootVH,Vvh,mhcore) # k for iteration
            psiVHrcoresquare=(rvh**2)*modulus(psi_VH_core(1,lwave,kvhcoreMod,rvh))*4*math.pi
            return(psiVHrcoresquare)
        #
        elif abs(rvh)>rcore and abs(rvh)<=rshell:
            VvhshellOO=Vv(rvh,rcore,rshell) #potential for the iteration
            VvhcoreOO=Vv(rcore,rcore,rshell) #constant potential of core
            kvhcoreOO=k(ErootVH,VvhcoreOO,mhcore)
            kvhshellOO=k(ErootVH,VvhshellOO,mhshell)
            AvhshellOO=A_VH_shell(1,lwave,kvhcoreOO,kvhshellOO,rcore,mhcore,mhshell)
            BvhshellOO=B_VH_shell(1,lwave,kvhcoreOO,kvhshellOO,rcore,mhcore,mhshell)
            psiVHrshellsquare=(rvh**2)*modulus(psi_VH_shell(AvhshellOO,BvhshellOO,lwave,kvhshellOO,rvh))*4*math.pi
            return(psiVHrshellsquare)
        #
        else: return(0)
    #
    #now calling the numerical integration function
    nmpsiVHnorm=path+"\\Wavefunctions\\"+f'Normalisation of psiVHcore.txt'
    fpsiVHnorm=open(nmpsiVHnorm,"w")
    fpsiVHnorm.close()
    IntVHcore=simpson(psicoreVHmodsquare,r0,rshell,hr,nmpsiVHnorm)  
    AcoreVH=math.sqrt(1/IntVHcore)
    print(f'Normalization constant AcoreVH = {AcoreVH}\n')
    #
    #now generating the complete wavefunction for given lwave, CE for a given energy value =>ErootCE
    # we will also generate |psi|^2
    nmVHwave=path+"\\Wavefunctions\\"+f'VH_wavefunc_l={int(lwave)}.txt'
    fVHwave=open(nmVHwave,"w")
    fVHwave.close()
    for rvals in range(Nr):
        r=r0+rvals*hr
        if abs(r)<=rcore:
            VvhcoreOO=Vv(r,rcore,rshell) #potential for the iteration
            kvhcoreOO=k(ErootVH,VvhcoreOO,mhcore)
            psiVHrcore=psi_VH_core(AcoreVH,lwave,kvhcoreOO,r)
            append_file(nmVHwave,f'{r/confac} {psiVHrcore.real*(wavefuncsqnmcon**(1/2))} {psiVHrcore.imag*(wavefuncsqnmcon**(1/2))} {modulus(psiVHrcore)*(wavefuncsqnmcon)}\n')
        #
        elif abs(r)>rcore and abs(r)<=rshell:
            VvhshellOO=Vv(r,rcore,rshell) #potential for the iteration
            VvhcoreOO=Vv(rcore,rcore,rshell) #constant potential of core
            kvhcoreOO=k(ErootVH,VvhcoreOO,mhcore)
            kvhshellOO=k(ErootVH,VvhshellOO,mhshell)
            AvhshellOO=A_VH_shell(AcoreVH,lwave,kvhcoreOO,kvhshellOO,rcore,mhcore,mhshell)
            BvhshellOO=B_VH_shell(AcoreVH,lwave,kvhcoreOO,kvhshellOO,rcore,mhcore,mhshell)
            psiVHrshell=psi_VH_shell(AvhshellOO,BvhshellOO,lwave,kvhshellOO,r)
            append_file(nmVHwave,f'{r/confac} {psiVHrshell.real*(wavefuncsqnmcon**(1/2))} {psiVHrshell.imag*(wavefuncsqnmcon**(1/2))} {modulus(psiVHrshell)*(wavefuncsqnmcon)}\n')
        #
        else: append_file(nmVHwave,f'{r/confac} {0} {0} {0}\n')
    #
#

# ------------------------------------------------ #
end=time.time()
print(f'Time taken for code to run = {end-begin}\n')




# ------------------------- unsed part in main code ----------------------------- #

#first, we are purely testing if the spherical bessel and spherical harmonics generated are correct or not. (VERIFIED)
#lmax=4
#to generate spherical bessel
#for l in range(lmax):
#    #defining names and locations of files to generate
#    nmspbessel1=path+"\\Wavefunctions\\" + "\\SpBessel1\\" +f'SphericalBessel1_l={l}.txt'
#    nmspbessel2=path+"\\Wavefunctions\\"+ "\\Spbessel2\\" +f'SphericalBessel2_l={l}.txt'
#    fsphb1=open(nmspbessel1,"w")
#    fsphb2=open(nmspbessel2,"w")
#    fsphb1.close()
#    fsphb2.close()
#    for rpoint in range(Nr):
#        r=r0+rpoint*hr
#        append_file(nmspbessel1, f'{r} {spherical_jn(l,r)}\n')
#        append_file(nmspbessel2, f'{r} {spherical_yn(l,r)}\n')
#    #
#
# we will now test the spherical harmonic function 
#lmax=3
#theta0=0.0
#phi0=0.0
#N=1000
#htheta=math.pi/N
#hphi=2*math.pi/N
#for l in range(lmax):
#    for m in range(-l,l+1):
#        nmsphar=path+"\\Wavefunctions\\"+"\\SpHarmonic\\"+f'SpHar_l={l},m={m}.txt'
#        fsphar=open(nmsphar,"w")
#        fsphar.close()
#        for inc in range(N):
#            theta=theta0+inc*htheta
#            phi=phi0
#            append_file(nmsphar,f'{theta} {phi} {modulus(scipy.special.sph_harm(m,l,phi,theta))}\n')
#        #
#    #
#
# Then, we find exact smallest root using Newton-Rhapson method
