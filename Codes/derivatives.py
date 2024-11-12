#finding derivatives

# DERIVATIVE
def derivative(f, x, h): # arguments = equation, point at which derivative is to be found, h
    y=(f(x+h)-f(x-h))/(2*h)
    return(y)

##########################################################################

# DOUBLE DERIVATIVE
def double_derivative(f, x, h): # arguments = equation, point at which double derivative is to be found, h
    y=(f(x+h)+f(x-h)-2*f(x))/(2*h*h)
    return(y)

###################################################################################
