#Library containing codes required to solve numerical integration

import random
import math
import handling_files
import derivatives

#Midpoint Method
def midpoint(f, lower, upper, N, name):
    h=(upper-lower)/N
    I=0
    for i in range(1,N): #if divided into N equal parts, then it has N-1 midpoints
        m=(lower+h*(i-1)+lower+h*(i))/2
        I=I+h*f(m)
        #
    # storing N vs I value in file
    handling_files.append_file(name, f'{N} {I}\n')
    return(I, h)

########################################################################

#Trapezoidal Method
def trapezoidal(f, lower, upper, N, name):
    h=(upper-lower)/N
    I=0
    X=[0 for i in range(N+1)] #for storing N+1 values, including upper and lower
    X[0]=lower
    for i in range(1,N+1):
        X[i]=lower+h*i #finding next point
        #finding integration value
        I=I+(h/2)*(f(X[i-1])+f(X[i]))
        #
    # storing N vs I value in file
    handling_files.append_file(name, f'{N} {I}\n')
    return(I, h)

#########################################################################

#Simpson Method
def simpson(f, lower, upper, N, name):
    h=(upper-lower)/N
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
    handling_files.append_file(name, f'{N} {I}\n')
    return(I, h)

##############################################################################

#Midpoint method with error given
def midpoint_error(f, lower, upper, err, name):
    f_dd=float(input('Please input the maximum value of second derivative of the function you want to integrate.\n'))
    N_cal=math.sqrt((math.pow(upper-lower,3)/(24*err))*abs(f_dd))
    N=math.ceil(N_cal)
    return midpoint(f, lower, upper, N, name)
    #

###############################################################################

#Trapezoidal method with error given
def trapezoidal_error(f, lower, upper, err, name):
    f_dd=float(input('Please input the maximum value of second derivative of the function you want to integrate.\n'))
    N_cal=math.sqrt((math.pow(upper-lower,3)/(12*err))*abs(f_dd))
    N=math.ceil(N_cal)
    return trapezoidal(f, lower, upper, N, name)
    #

###############################################################################

#Simpson method with error given
def simpson_error(f, lower, upper, err, name):
    f_d4=float(input('Please input the maximum value of fourth derivative of the function you want to integrate.\n'))
    N_cal=math.pow((math.pow(upper-lower,5)/(180*err))*abs(f_d4),1/4)
    N=math.ceil(N_cal)
    return simpson(f, lower, upper, N, name)
    #


################################################################################

#Monte Carlo Method
def montecarlo_uniform(f, lower, upper, N, name_I, name_err):
    F=0 #for storing the sum of f(X)
    F2=0 #for storing sum of [f(X)]^2
    for i in range(1,N+1): #loops from 1 to N
        R = random.random()  # random number between 0 and 1
        X = lower + (upper - lower) * R  # converting to random number between lower and upper
        F=F+f(X)
        F2=F2+math.pow(f(X), 2)
        #
    I=((upper-lower)/N)*F
    err=math.sqrt((1/N)*F2-math.pow((1/N)*F, 2))
    handling_files.append_file(name_I, f'{N} {I}\n')
    handling_files.append_file(name_err, f'{N} {err}\n')
    return(I, err)

#########################################################################################

#Trapezoidal Method in 1D for multivariable function
def trapezoidal_multi_z(f, lower, upper, N, name):
    h=(upper-lower)/N
    I=0
    X=[0 for i in range(N+1)]
    Y=[0 for i in range(N+1)]
    Z=[0 for i in range(N+1)] #for storing N+1 values, including upper and lower
    Int=[[0 for j in range(N+1)] for i in range(N+1)]
    X[0]=lower
    Y[0]=lower
    Z[0]=lower
    #finding the points
    for i in range(1,N+1):
        X[i]=lower+h*i
        Y[i]=lower+h*i
        Z[i]=lower+h*i
        #
    for i in range(N+1):
        for j in range(N+1):
            I=0
            for k in range(1,N+1):
                I = I + (h / 2) * (f(X[i], Y[j], Z[k - 1]) + f(X[i], Y[j], Z[k]))
                #
            Int[i][j]=I
            #
        #

    # storing N vs I value in file
    for i in range(N+1):
        for j in range(N+1):
            handling_files.append_file(name, f'{X[i]} {Y[j]} {Int[i][j]}\n')
            #
        #
    return(X,Y,Int)
