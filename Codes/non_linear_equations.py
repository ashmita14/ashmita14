#Library containing functions required to solve various systems of linear equations

import math
import handling_files
import derivatives
#def polynomial(X, )

# BRACKETING FUNCTION
def bracketing(f, a, b, tot_itr): #arguments = function, guesses a and b, and total number of allowed iterations
    par = 1.5  # factor to be multiplied with (b-a)
    diff = b - a
    if f(a) * f(b) < 0:  # bracketing complete
        return (a, b, True)  # returns True so that main function can analyse bracketing is complete
    if f(a) * f(b) > 0:
        if (tot_itr > 0):  # max iterations not reached yet
            if abs(f(a)) < abs(f(b)): return(bracketing(f, a-par*diff, b, tot_itr-1))  # a on right of root, so shifted to left
            if abs(f(a)) > abs(f(b)): return(bracketing(f, a, b+par*diff, tot_itr-1)) # b on left of root, so shifted to right
        else:
            return(a, b, False)
        #
    #

#######################################################################################

# BISECTION FUNCTION
def bisection_1(f, a, b, max_itr, nm): #arguments ==> formula, trial bracketing values, maximum iterations, name of file to store errors
    err=pow(10, -6)
    mid=[0 for i in range(max_itr)]
    n=0
    while n<max_itr:
        mid[n]=(a+b)/2
        #error = x_{n+1}-x_{n}, or equivalently, x_{n}-x_{n-1}
        if n>0: handling_files.append_file(nm, f'{n} {abs(mid[n] - mid[n - 1])}\n') #to avoid garbage value of mid[-1]
        if abs(b-a)<err: return(mid[n], True) #root found
        else:
            if f(a) * f(mid[n]) < 0: b = mid[n]  # root lies left of mid
            if f(a) * f(mid[n]) > 0: a = mid[n]  # root lies right of mid
            #
        n+=1
        #
    return(mid[n-1], False) #incase root not reached, returns last found value


###############################################################################

# REGULAR-FALSI
def regula_falsi(f, a, b, max_itr, x): #arguments ==> formula, trial bracketing values, maximum iterations, name of file to store errors
    err = pow(10, -6)
    n=0 #counter
    c=[0 for i in range(max_itr)]
    while n<max_itr: #so that it does not exceed max allowed iterations
        c[n]=b-(((b-a)*f(b))/(f(b)-f(a)))
        # error = x_{n+1}-x_{n}, equivalently x_{n}-x_{n-1}
        if n>0: handling_files.append_file(x, f'{n} {abs(c[n] - c[n - 1])}\n') #to avoid garbage value of c[-1]
        if abs(c[n]-c[n-1])>err:
            if f(a)*f(c[n])<0: b=c[n]
            if f(a)*f(c[n])>0: a=c[n]
            #
        else:
            return(c[n], True) #root found
        n+=1
        #
    return(c[n-1], False) #incase root not reached, returns last found value


######################################################################################

# NEWTON RAPHSON
def newton_raphson(f, a, max_itr, nm): #arguments ==> formula, trial root, maximum iterations, name of file to store errors
    err=pow(10, -6)
    n=0
    X=[0 for i in range(max_itr)]
    while n<max_itr: #so that it does not cross max number of allowed iterations
        if n==0:
            X[n]=a
        else:
            h=X[n-1]*0.001
            X[n]=X[n-1]-(f(X[n-1])/derivatives.derivative(f, X[n-1], h))
            #error = x_{n+1}-x_{n} o equivalently, x_{n}-x_{n-1}
            handling_files.append_file(nm, f'{n} {abs(X[n] - X[n - 1])}\n') #appends absolute error value with iteration number to file
            if abs(X[n]-X[n-1])<err: return(X[n], True)
            #
        n+=1
        #
    return(X[n-1], False) #if root not obtained even after max_iterations

######################################################################################

#  LAGUERRE METHOD
def laguerre(P, trial, max_itr, n): #arguments ==> formula, trial root, maximum iterations, name of file to store errors
    err=pow(10, -6)
    if round(P(trial), 6) !=0: #if trial is already the root
        G = (derivatives.derivative(P, trial, (0.001))) / P(trial)
        H = G * G - derivatives.double_derivative(P, trial, 0.001) / P(trial)
        if abs(G + math.sqrt((n - 1) * (n * H - G*G))) > abs(G - math.sqrt((n - 1) * (n * H - G*G))):
            a = n / (G + math.sqrt((n - 1) * (n * H - G*G)))
            #
        else:
            a = n / (G - math.sqrt((n - 1) * (n * H - G*G)))
            #
        if abs(a) < err: #a = trial_{n}-trial_{n+1} ; as trial_{n+1} = trial_{n}-a
            return(trial, True)
        else:
            if max_itr>0: return (laguerre(P, trial - a, max_itr-1, n))
            else: return(trial, False) #root not found in max_iterations
    else:
        return(trial, True)
    #

################################################################################

# DEFLATION BY SYNTHETIC DIVISION
def deflation(X, root, n): #arguments = coeffiecient matrix, root, order of polynomial
    f1=''
    while n>0: #first coefficient remains unchanged
        X[n-1]=X[n-1]+root*X[n] #coefficient of nth root in nth position and so on
        if X[n]<0: f1=f1+f'{X[n]}*pow(x, {n-1})'
        if X[n]>0: f1=f1+f'+{X[n]}*pow(x, {n-1})'   #building the new polynomial is string
        n-=1
        #
    r=X[0] #remainder
    X.pop(0) #removing remainder from array, and returning it separately
    if r<0: f1=f1+f'{r}'
    if r>0: f1=f1+f'+{r}'
    f=lambda x: eval(f1) #converting string of new polynomial to formula
    return(f, X, r)
