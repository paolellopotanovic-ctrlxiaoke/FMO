#!/usr/bin/env python

import numpy as np

H2kcal_mol=627.509

f=open("temp1",'r')
lines1=f.readlines()
f.close()

n=len(lines1)
e1c=np.zeros((n,))
err1=np.zeros((n,))

f=open("temp2",'r')
lines2=f.readlines()
f.close()

e1_0c_sum=0.
err_sum=0.
for line1 in lines1:
   fname=line1[:-1]+"_en"
   f=open(fname,'r')
   lines0=f.readlines()
   f.close()
   t1,i,t2,t3=line1.split("_")
   i=int(i)
   i=i-1
   e1=float(lines0[0].split()[1])
   e1_0=float(lines0[0].split()[0])
   err1[i]=float(lines0[0].split()[2])
   e1c[i]=e1-e1_0
   e1_0c_sum+=e1c[i]
   err_sum+=err1[i]

e2_0c_sum=0.
for line2 in lines2:
   fname=line2[:-1]+"_en"
   f=open(fname,'r')
   lines0=f.readlines()
   f.close()
   t1,i,j,t2=line2.split("_")
   i,j=map(int,(i,j))
   i,j=i-1,j-1
   e2=float(lines0[0].split()[1])
   e2_0=float(lines0[0].split()[0])
   err2=float(lines0[0].split()[2])
   e2c=(e2-e2_0)-e1c[i]-e1c[j]
   e2_0c_sum+=e2c
   err_sum+=err2+err1[i]+err1[j]


e,err=e1_0c_sum+e2_0c_sum,err_sum
print "%.2f +/- %.2f"%(round(H2kcal_mol*e,2),round(H2kcal_mol*err,2))
