#Library for handling files and their contents

# READ FILE
def read_matrix(x): #more than one column #parameter = name of file
    f=open(x,'r') #'r' ==> read only
    X=[[float(num) for num in line.split(' ')] for line in f]
    f.close()
    return(X)

##############################################################################

# APPEND FILE
def append_file(x, str): #arguments = name of file, string to append
    f=open(x, 'a') #'a' ==> append file
    f.write(str)
    f.close()
    #

###############################################################################

# PRINT CONTENTS OF A TEXT FILE
def print_file(x): #argument = name of file
    f=open(x, 'r')
    contents=f.read()
    print(contents)

#################################################################################
