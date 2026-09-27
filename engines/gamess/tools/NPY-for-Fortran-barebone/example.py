#!/usr/bin/env python3
import numpy as np

a_soll =  np.zeros((10,20), dtype=np.double)
for ind in range(a_soll.shape[0]):
    for jnd in range(a_soll.shape[1]):
        i = ind + 1
        j = jnd + 1
        a_soll[ind,jnd] = i * j

b_soll =  np.zeros(10, dtype=np.double)
for ind in range(b_soll.shape[0]):
    i =  ind + 1
    b_soll[ind] = 2 * i

a = np.load("a.npy")
b = np.load("b.npy")

print("A: ")
print( np.max(np.abs(a - a_soll) / a_soll) )

print("B: ")
print( np.max(np.abs(b - b_soll) / b_soll) )
