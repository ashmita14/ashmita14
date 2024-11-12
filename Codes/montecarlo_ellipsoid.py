# Volume of Ellipsoidal by Monte Carlo Method

import random
import math
import handling_files


def ellipsoidal(N, a, b, c, name):
    N_hit=0 #counting number of points in the ellipsoid
    for i in range(N):
        # generating random variables in x-axis between -a to a
        x = -a + 2 * a * random.random()
        # generating random variables in y-axis between -b to b
        y = -b + 2 * b * random.random()
        # generating random variables in z-axis between -c to c
        z = -c + 2 * c * random.random()
        if (x**2/a**2)+(y**2/b**2)+(z**2/c**2) <= 1:
            N_hit=N_hit+1
            handling_files.append_file(name, f'{x} {y} {z}\n')
            #
        #
    vol=(N_hit/N)*(2*a)*(2*b)*(2*c)
    return (vol)


