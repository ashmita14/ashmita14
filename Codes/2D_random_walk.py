#2 Dimensional Random Walk

import random
import math
import handling_files

#Generating the Random Walk, and calculating final radial distance
def walk(N, name):
    #starting values
    x=0
    y=0
    handling_files.append_file(name, f'0 0\n')
    for i in range(N): #generating 250 steps, with fixed step size of 1
        #chosing random direction
        theta=2*math.pi*random.random()
        x=x+math.cos(theta)
        y=y+math.sin(theta)
        handling_files.append_file(name, f'{x} {y}\n')
        #
    r=math.sqrt(x*x+y*y)
    return(r, x, y) #returning radial distance and displacements along x and y (final values of x and y)
    #

###################################################################################

#Calculating Rrms, avg displacement along x and y
def averages(n, R, X, Y):
    sum_R=0.0
    sum_x=0.0
    sum_y=0.0
    for i in range(n):
        sum_R=sum_R+math.pow(R[i], 2)
        sum_x=sum_x+X[i]
        sum_y=sum_y+Y[i]
        #
    Rrms=math.sqrt(sum_R/n)
    avgx=sum_x/n
    avgy=sum_y/n
    return(Rrms, avgx, avgy)


########################################################################################
