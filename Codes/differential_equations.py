#Codes for solving differential equations

import math
import handling_files

#FORWARD EULER (EXPLICIT EULER)
def forward_euler(f,x0,y0,h,N,name):
    x=x0
    y=y0
    for i in range(N):
        x=x+h
        K1=f(y,x)
        y=y+h*K1
        handling_files.append_file(name,f'{x} {y}\n')
        #
    #


#################################################################################

#RUNGE KUTTA 4 (RK4) for second order differential equation (double derivatives)
def RK4_double(F1,F2,x0,y0,dy0,h,N,name):
    x=x0
    y=y0
    dy=dy0
    handling_files.append_file(name,f'{x} {y}\n')
    for i in range(N):
        K1Y=F1(dy,y,x)
        K1DY=F2(dy,y,x)
        #
        K2Y=F1(dy+(K1DY*h)/2, y+(K1Y*h)/2, x+h/2)
        K2DY=F2(dy+(K1DY*h)/2, y+(K1Y*h)/2, x+h/2)
        #
        K3Y=F1(dy+(K2DY*h)/2, y+(K2Y*h)/2, x+h/2)
        K3DY=F2(dy+(K2DY*h)/2, y+(K2Y*h)/2, x+h/2)
        #
        K4Y=F1(dy+(K3DY*h), y+(K3Y*h), x+h)
        K4DY=F2(dy+(K3DY*h), y+(K3Y*h), x+h)
        #
        dy=dy+((K1DY+2*K2DY+2*K3DY+K4DY)*h)/6
        y=y+((K1Y+2*K2Y+2*K3Y+K4Y)*h)/6
        x=x+h
        handling_files.append_file(name,f'{x} {y}\n')
        #
    return(y)


####################################################################

#DIRICHELET BOUNDARY
def boundary_dirichlet(F1,F2,x0,y0,yn,h,N,name):
    guess = float(input(f'Enter guess value of slope of curve at x={x0 + h} when y({x0})={y0}\n'))
    y=RK4_double(F1,F2,x0,y0,guess,h,N,name)
    err=math.pow(10,-4)
    if y>yn: #if first guess is greater than the actual value
        guess_high=guess
        guess=float(input('Guess value of slope again, lesser than the previous guess.\n'))
        file = open(name, 'w')
        y=RK4_double(F1,F2,x0,y0,guess,h,N,name)
        while y>yn: #looping until a lower guess is found
            guess_high = guess
            guess = float(input('Guess value of slope again, lesser than the previous guess.\n'))
            file = open(name, 'w')
            y = RK4_double(F1, F2, x0, y0, guess, h, N, name)
            #
        if y==yn: return(guess,True) #if guess gives perfect value of yn at xn, then solution is obtained
        else: guess_low=guess
        #
    elif y<yn: #if first guess is lesser than the actual value
        guess_low=guess
        guess=float(input('Guess value of slope again, greater than the previous guess.\n'))
        file = open(name, 'w')
        y = RK4_double(F1, F2, x0, y0, guess, h, N, name)
        while y < yn: #looping unntil a high guess is found
            guess_high = guess
            guess = float(input('Guess value of slope again, greater than the previous guess.\n'))
            file = open(name, 'w')
            y = RK4_double(F1, F2, x0, y0, guess, h, N, name)
            #
        if y==yn: return (True) #if guess gives perfect value of yn at xn, then solution is obtained
        else:
            guess_high = guess
        #
    else: return(guess,True) #if guess gives perfect value of yn at xn, then solution is obtained
    #
    #finding y values for high and low value of guess slope
    file = open(name, 'w')
    y_low=RK4_double(F1, F2, x0, y0, guess_low, h, N, name)
    file = open(name, 'w')
    y_high=RK4_double(F1, F2, x0, y0, guess_high, h, N, name)
    #setting upper limit on iteration
    limit=int(input('Enter the upper limit on the number of iteration for finding the correct value of the slope.\n'))
    i=0 #counter
    while round(abs(y-yn),4)>err and i<limit:
        guess=guess_low+((guess_high-guess_low)/(y_high-y_low))*(yn-y_low)
        file=open(name,'w')
        y=RK4_double(F1, F2, x0, y0, guess, h, N, name)
        i=i+1 #incrementing counter by 1
        #
    if round(abs(y-yn),4)>err: return (False) #solution not found
    else: return(guess,True) #solution found




