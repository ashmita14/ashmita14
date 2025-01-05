#python code to obtain electrostatic energies, generating polynomials through recurrence relation and changes Al's and Bl's to prevent blowups.
#Also, this code attempts to parallelize the functioning

#first we will import some basic libraries
import math
import os
import time
import sys
import numpy as np
import multiprocessing

path=os.path.dirname(os.path.realpath(__file__))
sys.path.append(path)

#measuring beginning time to find time taken by code to run
begin=time.time()

#defining additional functions which might be needed
#
#Legendre Polynomial of given order
def P(n,x,Plnm1,Plnm2):
    if n>=2 : Pln=(1/n)*((2*n-1)*x*Plnm1-(n-1)*Plnm2)
    if n==1 : Pln=Plnm1
    if n==0 : Pln=Plnm2
    return Pln
#

#declare fixed constants, epsilon0 and position of the dot on z-axis
epsl0=8.85*(10**(-18))
zdot=0.2

#defining polarizability of complete CdSe-CdS dot ==> all in SI units
epsl0SI=epsl0*(10**6) #epsilon naught in SI units
b=2.9*(10**(-9)) #radius of quantum dot in meters
a=1.5*(10**(-9)) # radius of core of quantum dot in meters
nCDSEref=2.3847 #real part of refractive index of CdSe for 1 micron wavelength
kCDSEref=0.065147 #imaginary part of refractive index of CdSe for 1 micron wavelength
nCDSref=2.307 #real part of refractive index of CdS for 1 micron wavelength
kCDSref=0.0 #imaginary part of refractive index of CdS for 1 micron wavelength
epslCdSeR=((nCDSEref**2)-(kCDSEref**2))*epsl0SI #real part of dielectric constant of CdSe dot
epslCdSR=((nCDSref**2)-(kCDSref**2))*epsl0SI
numeff=(b**3)*(epslCdSeR+2*epslCdSR)+ 2*(a**3)*(epslCdSeR-epslCdSR)
denomeff=(b**3)*(epslCdSeR+2*epslCdSR)-(a**3)*(epslCdSeR-epslCdSR)
epslEff=epslCdSR*(numeff/denomeff)
alpha=(4*math.pi*epsl0SI*(b**3))*((epslEff-epsl0SI)/(epslEff+2*epsl0SI)) #polarizability for CdSe-CdS dot 

#b=10*(10**(-9)) #radius of quantum dot in meters
#nref=2.3847 #real part of refractive index of CdSe for 1 micron wavelength
#kref=0.065147 #imaginary part of refractive index of CdSe for 1 micron wavelength
#epslCdSeR=((nref**2)-(kref**2))*epsl0SI #real part of dielectric constant of CdSe dot
#alpha=(4*math.pi*epsl0SI*(b**3))*((epslCdSeR-epsl0SI)/(epslCdSeR+2*epsl0SI)) #polarizability for CdSe dot (not shell structure yet)

#defining A's, B's, Erho's and Ez's as functions
def A(Q,l): return (Q/(8*math.pi*epsl0))*(((-1)**(l/2))/(l+1)) #only for even numbers
def B(Q,l): return (Q/(8*math.pi*epsl0))*(((-1)**((l+1)/2))/(l)) #only for odd numbers
#
def Erhoin(Q,a,l,rho,z,Pln,Plnm1):
    #cscth=abs(z/(math.sqrt(z**2 + rho**2))) #cos theta term ==> coefficient of the Legendre polynomials
    return B(Q,l)*l*(((rho**2 + z**2)/(a**2))**(l/2))*(1/((rho**2 +z**2)*a))*((z/rho)*math.sqrt(z**2 + rho**2)*Plnm1-((rho**2 + z**2)/rho)*Pln)
#
def Ezin(Q,a,l,rho,z,Pln,Plnm1):
    #cscth=abs(z/(math.sqrt(z**2 + rho**2))) #cos theta term ==> coefficient of the Legendre polynomials
    return -B(Q,l)*l*(((rho**2 + z**2)/(a**2))**(l/2))*(1/(a*math.sqrt(rho**2 + z**2)))*Plnm1
#
def Erhoout(Q,a,l,rho,z,Pln,Plnm1):
    cscth=abs(z/(math.sqrt(z**2 + rho**2))) #cos theta term ==> coefficient of the Legendre polynomials
    return A(Q,l)*(((a**2)/(rho**2 + z**2))**(l/2))*((rho**2 +z**2)**(-3/2))*(((l*(rho**2 - z**2)+rho**2)/rho)*Pln+l*(z/rho)*math.sqrt(z**2 + rho**2)*Plnm1)
#
def Ezout(Q,a,l,rho,z,Pln,Plnm1):
    cscth=abs(z/(math.sqrt(z**2 + rho**2))) #cos theta term ==> coefficient of the Legendre polynomials
    return A(Q,l)*(((a**2)/(rho**2 + z**2))**(l/2))*((rho**2 +z**2)**(-3/2))*(z*(2*l+1)*Pln-l*math.sqrt(z**2 + rho**2)*Plnm1)
#

#next, we need to import grid data from grid text file
N=3000 #number of disks
nmdisk=path+'\\Results\\'+f'DiskGrid(N={N}).txt' #information of the disks
fd=open(nmdisk,'r')
gridarray=[[float(num) for num in line.split(' ')] for line in fd] #read entire file
numrows=len(gridarray) #number of rows
numcol=len(gridarray[0]) #number of columns
#defining empty arrays
X=[0 for i in range(numrows)]
Y=[0 for i in range(numrows)]
R=[0 for i in range(numrows)]
V0=[0 for i in range(numrows)]
Q=[0 for i in range(numrows)]
#storing values appropriately
for i in range(numrows):
    X[i]=gridarray[i][0]
    Y[i]=gridarray[i][1]
    R[i]=gridarray[i][3]
    V0[i]=gridarray[i][4]
    Q[i]=gridarray[i][5]
#
#next, we need to define the trajectory of the dot moving with zdot constant
#
#If we have a straight and parallel line, this section of the code will be used. We need to choose the type of line. specify choice by converting desired line variable to "True"
parallelx=True
parallely=False
slopeline=False
#
# Now, we will define a function to generate results for one line at a time. We will then call this function using a multiprocessing tool. 
def MultiLine(ln):
    # If lines parallel to x-axis (constant y lines) ==> parallelx==True
    if parallelx==True:
        c0=0.0; c1=1.0 #coefficients for x (c0 gives value of x for constant x line)
        d0=lowlim+prec*ln; d1=0.0 #coefficients of y (d0 gives value of y for constant y line) 
        #defining file name
        nmener=path+'\\Results\\'+'\\N=3000(old)\\'+f'Energy(N={N},y0={round(d0,1)},Constant Y line).txt' #energy data file
    #
    # If lines parallel to y-axis ==> parallely==True
    if parallely==True:
        c0=lowlim+prec*ln; c1=0.0 #coefficients for x (c0 gives value of x for constant x line)
        d0=0.0; d1=1.0 #coefficients of y (d0 gives value of y for constant y line)
        #defining file name
        nmener=path+'\\Results\\'+'\\N=3000(old)\\'+f'Energy(N={N},x0={round(c0,1)},Constant X line).txt' #energy data file
    #
    # If line has a non-zero slope (and non-infinity) ==> slopeline==True
    if slopeline==True:
        m=1.0 #slope of line
        c=0.0 #y-intercept
        c0=0.0; c1=1.0 #coefficients for x (c0 also determines x-intercept)
        d0=m*c0+c; d1=m*c1 #coefficients for y
        #defining file name
        nmener=path+'\\Results\\'+f'Energy(N={N},m={m},x0={c0},y0={d0},Sloped line).txt' #energy data file
    #
    t0=0.0; tmax=100.0 #limits on t
    inc=0.05
    t=t0
    #
    #defining file to store values of energy in
    fn1=open(nmener,"w")
    fn1.close()
    fn2=open(nmener,"a")
    #fn2.write('Xdot Ydot Zdot tDot Udot\n')
    #
    #I am choosing to define it as a straight line. ==> This is the first (and outermost loop).
    while t<=tmax:
        #calculate x and y from t
        x=round(c0+c1*t,3)
        y=round(d0+d1*t,3)
        #calculate rho for dot
        rhodot=math.sqrt((x**2)+(y**2))
        #defining total electric field 
        Etotx=0.0
        Etoty=0.0
        Etotz=0.0
        Utot=0.0
        #now we need to run a loop through all disks ==> This is the second loop
        for i in range(N):
            #rholim=(R[i]**2) - (zdot**2)
            rhodotdisk=((X[i]-x)**2)+((Y[i]-y)**2) #rho squared of point with respect to disk
            disdotdisk=rhodotdisk+(zdot**2) #actual distance of dot from center of disk
            rhodotdisksqrt=math.sqrt(rhodotdisk)
            Etotrhodisk=0.0
            Etotzdisk=0.0
            if disdotdisk<=(R[i]**2): #in this condition we have B's
                Pl=[]
                Pl.append(1)
                Pl.append(zdot/(math.sqrt(disdotdisk)))
                #
                l=1 #beginning count on number of terms
                if disdotdisk<=(0.000001*(R[i]**2)): #doing an additional cutoff check ==> if the point is very close to center, only one term is enough
                    Ezinterm=Ezin(Q[i],R[i],l,rhodotdisksqrt,zdot,Pl[1],Pl[0])
                    Erhointerm=Erhoin(Q[i],R[i],l,rhodotdisksqrt,zdot,Pl[1],Pl[0])
                    Etotrhodisk+=Erhointerm
                    Etotzdisk+=Ezinterm
                #
                else:
                    Ein=False #variable to decide whether to end while loop to count number of terms
                    while Ein==False: #loop runs until sum is under a tolerable level
                        if l==1: Pln=Pl[1]; Plnm1=Pl[0]
                        else : 
                            Plnm1=P(l-1,zdot/math.sqrt(disdotdisk),Pl[l-2],Pl[l-3])
                            Pl.append(Plnm1)
                            Pln=P(l,zdot/math.sqrt(disdotdisk),Pl[l-1],Pl[l-2])
                            Pl.append(Pln)
                        #
                        Ezinterm=Ezin(Q[i],R[i],l,rhodotdisksqrt,zdot,Pln,Plnm1)
                        Erhointerm=Erhoin(Q[i],R[i],l,rhodotdisksqrt,zdot,Pln,Plnm1)
                        #print(x,y,l,Plnm1, f'Erhointerm={Erhointerm}', f'Erhoinsum={Etotrhodisk}', f'Ezterm={Ezinterm}', f'Ezsum={Etotzdisk}')
                        if abs(Erhointerm)<=(0.0001*abs(Etotrhodisk)) and abs(Ezinterm)<=(0.0001*abs(Etotzdisk)): Ein=True #if wihtin tolerance, end loop
                        Etotrhodisk+=Erhointerm
                        Etotzdisk+=Ezinterm
                        l+=2
                    #
                #   
            #
            else: #in this condition we have A's
                Pl=[]
                Pl.append(1)
                Pl.append(zdot/(math.sqrt(disdotdisk)))
                #
                l=0 #beginning count on number of terms
                if disdotdisk>=(10000*(R[i]**2)): #doing additiobnal cutoff check ==> if point is very far away from disk, use only one term
                    Erhooutterm=Erhoout(Q[i],R[i],l,rhodotdisksqrt,zdot,Pl[0],0)
                    Ezoutterm=Ezout(Q[i],R[i],l,rhodotdisksqrt,zdot,Pl[0],0)
                    Etotrhodisk+=Erhooutterm
                    Etotzdisk+=Ezoutterm
                #   
                else:
                    Ein=False #variable to decide whether to end while loop to count number of terms
                    while Ein==False: #loop runs until sum is under a tolerable level
                        #print(l,P(l,zdot/math.sqrt(disdotdisk)))
                        if l==0: Pln=Pl[0]; Plnm1=0
                        elif l==2: 
                            Plnm1=Pl[1]
                            Pln=P(l,zdot/math.sqrt(disdotdisk),Pl[l-1],Pl[l-2])
                            Pl.append(Pln)
                        #
                        else : 
                            Plnm1=P(l-1,zdot/math.sqrt(disdotdisk),Pl[l-2],Pl[l-3])
                            Pl.append(Plnm1)
                            Pln=P(l,zdot/math.sqrt(disdotdisk),Pl[l-1],Pl[l-2])
                            Pl.append(Pln)
                        #
                        Erhooutterm=Erhoout(Q[i],R[i],l,rhodotdisksqrt,zdot,Pln,Plnm1)
                        Ezoutterm=Ezout(Q[i],R[i],l,rhodotdisksqrt,zdot,Pln,Plnm1)
                        if abs(Erhooutterm)<=(0.0001*abs(Etotrhodisk)) and abs(Ezoutterm)<=(0.0001*abs(Etotzdisk)): Ein=True #if wihtin tolerance, end loop
                        Etotrhodisk+=Erhooutterm
                        Etotzdisk+=Ezoutterm
                        l+=2
                    #
                #
            #
            #Erho needs to be added vectorially, so we need to convert it to Ex and Ey
            Etotxdisk=((x-X[i])/(rhodotdisksqrt))*Etotrhodisk
            Etotydisk=((y-Y[i])/(rhodotdisksqrt))*Etotrhodisk
            Etotx+=Etotxdisk
            Etoty+=Etotydisk
            Etotz+=Etotzdisk
            #this gets repeated for all disks
        #
        Utot=alpha*((Etotx**2)+(Etoty**2)+(Etotz**2))*(10**(12)) #total energy at the point(x,y,zdot) due to all disks
        #storing the values in a file
        fn2.write(f'{round(x,3)} {round(y,3)} {zdot} {round(t,3)} {Utot} {Etotx*(10**6)} {Etoty*(10**6)} {Etotz*(10**6)}\n')
        #increments go in last
        t=t+inc
        #if blowup==True: fcon2.write(f'\n\nPoint Change\n')
    #
    fn2.close()
#
#
# Now that we have our function completeley defined, we will call the function from a multiprocessing unit.
# range of lines to generate, with precision
prec=0.2
uplim=0.0
lowlim=0.0
noflines=int((uplim-lowlim)/prec) #number of lines to generate
if __name__ == '__main__':
    pool_obj=multiprocessing.Pool(processes=5)
    pool_obj.map(MultiLine,range(0,noflines+1))
    pool_obj.close()
#
#measuring end time
end=time.time()
#printing time taken for code to run
print(f'Time taken by code to run = {end-begin}')
