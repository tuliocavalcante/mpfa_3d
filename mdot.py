# This is the function to multiply an array of matrices and an array of arrays, but entrance to entrance.
import numpy as np

def mdot(a,B,c):

    t = []

    for i in range(len(a)):
        t.append(sum(a[i]*(B[i]@c[i])))

    t = np.array(t)

    return t
