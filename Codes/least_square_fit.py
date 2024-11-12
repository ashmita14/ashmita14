# Least Square Fitting Code
import math
import lineareqns
import handling_files

def least_square(M,nm):
    n=len(M) #no. of elements
    #counters
    x1=0
    x2=0
    y1=0
    x1y1=0
    #summing over xi, yi, xi^2 and xiyi
    for i in range(n):
        x1=x1+M[i][0]
        x2=x2+math.pow(M[i][0],2)
        y1=y1+M[i][1]
        x1y1=x1y1+(M[i][0]*M[i][1])
        #
    handling_files.append_file(nm, f'{n} {x1}\n')
    handling_files.append_file(nm,f'{x1} {x2}\n')
    Y=[0,0]
    #defining column vector [y1 x1y1]
    Y[0]=y1
    Y[1]=x1y1
    #reading the A1 matrix from file, A1*Coefficient Matrix = Y
    A1=handling_files.read_matrix(nm)
    A=lineareqns.LU_decomposition(A1) #LU decomposing A1 matrix
    #obtaining Coefficient matrix from forward-backward matrix operation
    Coef=lineareqns.forward_backward(A,Y)
    return(Coef)

def pearson(M):
    n=len(M) #no. of elements
    #counters
    x1=0
    y1=0
    #finding sums
    for i in range(n):
        x1 = x1 + M[i][0]
        y1 = y1 + M[i][1]
        #
    #finding averages
    xavg=x1/n
    yavg=y1/n
    sxx=0
    syy=0
    sxy=0
    #finding sxx, sxy, syy
    for i in range(n):
        sxx=sxx+math.pow((M[i][0]-xavg),2)
        syy=syy+math.pow((M[i][1]-yavg),2)
        sxy=sxy+(M[i][0]-xavg)*(M[i][1]-yavg)
        #
    #finding r^2 from given formula
    r2=sxy*sxy/(sxx*syy)
    return(r2)
