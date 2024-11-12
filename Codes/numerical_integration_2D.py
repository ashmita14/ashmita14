# Library containing codes required to solve numerical integration in two dimensional

import random
import math
import handling_files
import derivatives

#Trapeziodal Method in 2D
def trapezoidal_2D(f, lowx, upx, lowy, upy, Nx, Ny):
    h=(lowx-upx)/Nx
    k=(lowy-upy)/Ny
    X=[0 for i in range(Nx+1)]
    X[0]=lowx
    X[Nx]=upx
    Y=[0 for i in range(Ny+1)]
    Y[0]=lowy
    Y[Ny]=upy
    t1=0
    t2=0
    t3=0
    T1=((h*k)/4)*(f(lowx,lowy)+f(lowx,upy)+f(upx,lowy)+f(upx,upy))
    for i in range(1,Nx):
        X[i]=X[i-1]+h
        t1=t1+f(X[i],lowy)+f(X[i],upy)
        #
    for j in range(1,Ny):
        Y[j]=Y[j-1]+k
        t2=t2+f(lowx,Y[j])+f(upx,Y[j])
        #
    for i in range(1,Nx):
        for j in range(1,Ny):
            t3=t3+f(X[i],Y[j])
            #
        #
    T2=((h*k)/2)*(t1+t2)
    T3=(h*k)*t3
    I=T1+T2+T3
    return(I)
