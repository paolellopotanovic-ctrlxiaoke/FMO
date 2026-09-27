#!/usr/bin/env python

import sys, re

fname=sys.argv[1]
fname0=fname.split('.')[0]

def step1(fname0):

   key=re.compile(r'QMC CURRENT.*=\s*(\d+).*=\s*(\d+).*=\s*(\d+)')

   fname=fname0+".dat"
   f=open(fname)
   lines=f.readlines()
   f.close()

   f_t1=open("temp1",'w')
   f_t2=open("temp2",'w')

   mons=[]
   n=len(lines)
   for ind in range(n):
      line=lines[ind]
      match=key.search(line)
      if match:
         msg=match.groups()
         i,j,k=map(int,[msg[0],msg[1],msg[2]])
         nat=int(lines[ind+3][:-1].split("=")[1])
         dname="afo_"+str(i)+"_"+str(j)+"_"+str(k)
         if j==0:
            if not i in mons:
               print >>f_t1,dname
            mons.append(i)
         else:
            print >>f_t2,dname
         fname=dname+".inp"
         f=open(fname,'w')
         atoms=[]
         for ind2 in range(nat):
            pass
            v2=lines[ind+11+ind2].split()
            if j==0:
               a,c,x,y,z,e=v2[0][3:],float(v2[1]),float(v2[2]),float(v2[3]),float(v2[4]),float(v2[5])
               print >>f,a,c,x,y,z
            else:
               a,c,x,y,z=v2[0][3:],float(v2[1]),float(v2[2]),float(v2[3]),float(v2[4])
               print >>f,a,c,x,y,z
         f.close()
   f_t1.close()
   f_t2.close()

def step2(fname0):
   fname=fname0+".inp"
   f=open(fname,'r')
   text=f.read()
   f.close()

   mat1=re.search(r"indat\(1\)=([\,\d\s]+)[\$\w]",text,re.MULTILINE|re.DOTALL|re.IGNORECASE)
   mat2=re.search(r"fmobnd\s*(.+)\s*\$end",text,re.MULTILINE|re.DOTALL|re.IGNORECASE)

   if mat1: 
      indat=mat1.group(1)
      indat=indat.replace(" ","")
      indat=indat.split(",")
      indat=map(int,indat)

   if mat2:
      bnds0=mat2.group(1).split()
      bnds0=map(int,bnds0)
      bnds=[]
      flag=True
      for  b in bnds0:
         if flag:
            b1=-b
            flag=False
         else: 
            bnds.append((b1,b))
            flag=True

   infoc={}
   frag_len=[]
   ind0=0
   iold=0
   ind=0
   fcount={}
   for i in indat:
      ind0=ind0+1
      if not i==iold and ind0>1:
         fcount[iold]=ind
         ind=0
      ind+=1
      iold=i
      infoc[ind0]=(i,ind)
   fcount[iold]=ind

   f1=open("temp1",'r')
   lines1=f1.readlines()
   f1.close()

   for line1 in lines1:
      f0=open(line1[:-1]+".inp",'r')
      lines0=f0.readlines()
      f0.close()
      t1,i,t2,t3=line1.split('_')
      i=int(i)
      flag=0
      for b in bnds:
         j=infoc[b[0]][0]  
         k=infoc[b[0]][1]  
         l=infoc[b[1]][0]  
         m=infoc[b[1]][1]  
         if j==i:
            a,c,x,y,z=lines0[k-1].split()
            lines0[k-1]=a+" "+str(float(c)+1)+" "+x+" "+y+" "+z
            fn="afo_"+str(l)+"_0_0.inp"
            f=open(fn,'r') 
            lines=f.readlines()
            f.close()
            a,c,x,y,z,=lines[m-1].split()
            lines0.append("H 1.0 "+x+" "+y+" "+z)
         elif l==i:
            flag=1
      if flag:
         n=len(lines0)
         for ind in range(n):
            a,c,x,y,z=lines0[ind].split() 
            if int(float(c))==1:
               a="H "
               lines0[ind]=a+c+" "+x+" "+y+" "+z
      f=open(line1[:-1]+".inp",'w')
      for l0 in lines0:
         print >>f,l0[:-1]
      f.close()

   f2=open("temp2",'r')
   lines2=f2.readlines()
   f2.close()

   for line2 in lines2:
      f0=open(line2[:-1]+".inp",'r')
      lines0=f0.readlines()
      f0.close()
      t1,i1,i2,t2=line2.split('_')
      i1=int(i1)
      i2=int(i2)
      flag=0
      for b in bnds:
         j=infoc[b[0]][0]  
         k=infoc[b[0]][1]  
         l=infoc[b[1]][0]  
         m=infoc[b[1]][1]  
         if (j==i1 and (not l==i2)) or (j==i2 and (not l==i1)):
            if j==i2: k+=fcount[i1]
            a,c,x,y,z=lines0[k-1].split()
            lines0[k-1]=a+" "+str(float(c)+1)+" "+x+" "+y+" "+z
            fn="afo_"+str(l)+"_0_0.inp"
            f=open(fn,'r') 
            lines=f.readlines()
            f.close()
            a,c,x,y,z,=lines[m-1].split()
            lines0.append("H 1.0 "+x+" "+y+" "+z)
         elif (l==i1 and (not j==i1)) or (l==i2 and (not j==i1)):
            flag=1
      if flag:
         n=len(lines0)
         for ind in range(n):
            a,c,x,y,z=lines0[ind].split() 
            if int(float(c))==1:
               a="H "
               lines0[ind]=a+c+" "+x+" "+y+" "+z
      f=open(line2[:-1]+".inp",'w')
      for l0 in lines0:
         print >>f,l0[:-1]
      f.close()

def step3():

   f=open("templ_files/head")
   head_lines=f.readlines()
   f.close()

   f=open("templ_files/H.bas")
   Hbas_lines=f.readlines()
   f.close()

   f=open("templ_files/O.bas")
   Obas_lines=f.readlines()
   f.close()

   f=open("templ_files/C.bas")
   Cbas_lines=f.readlines()
   f.close()

   f=open("templ_files/N.bas")
   Nbas_lines=f.readlines()
   f.close()

   f=open("templ_files/H.ecp")
   Hecp_lines=f.readlines()
   f.close()

   f=open("templ_files/O.ecp")
   Oecp_lines=f.readlines()
   f.close()

   f=open("templ_files/C.ecp")
   Cecp_lines=f.readlines()
   f.close()

   f=open("templ_files/N.ecp")
   Necp_lines=f.readlines()
   f.close()

   f1=open("temp1",'r')
   lines1=f1.readlines()
   f1.close()

   f2=open("temp2",'r')
   lines2=f2.readlines()
   f2.close()

   for line1 in lines1:
      f=open(line1[:-1]+".inp",'r')
      lines=f.readlines()
      f.close()
      f=open(line1[:-1]+".inp",'w')
      for hline in head_lines:
         print >>f,hline[:-1]
      atoms=[]
      for line in lines:
         a,c,x,y,z=line[:-1].split()
         atoms.append(a)
         print >>f,line[:-1]
         if a=="H":
            for line in Hbas_lines:
               print >>f,line[:-1]
         if a=="O":
            for line in Obas_lines:
               print >>f,line[:-1]
         if a=="C":
            for line in Cbas_lines:
               print >>f,line[:-1] 
         if a=="N":
            for line in Nbas_lines:
               print >>f,line[:-1] 

      print >>f," $END"
      print >>f," $ECP"
      H_flag=0
      O_flag=0
      C_flag=0
      N_flag=0
      for a in atoms:
         if a=="H":
            if H_flag==0:
               for line in Hecp_lines:
                  print >>f,line[:-1]
               H_flag=1
            else:
               print >>f,"H-QMC GEN"
         elif a=="O":
            if O_flag==0:
               for line in Oecp_lines:
                  print >>f,line[:-1]
               O_flag=1
            else:
               print >>f,"O-QMC GEN"
         elif a=="C":
            if C_flag==0:
               for line in Cecp_lines:
                  print >>f,line[:-1]
               C_flag=1
            else:
               print >>f,"C-QMC GEN"
         elif a=="N":
            if N_flag==0:
               for line in Necp_lines:
                  print >>f,line[:-1]
               N_flag=1
            else:
               print >>f,"N-QMC GEN"
      print >>f," $END"
      f.close()

   for line2 in lines2:
      f=open(line2[:-1]+".inp",'r')
      lines=f.readlines()
      f.close()
      f=open(line2[:-1]+".inp",'w')
      for hline in head_lines:
         print >>f,hline[:-1]
      atoms=[]
      for line in lines:
         a,c,x,y,z=line[:-1].split()
         atoms.append(a)
         print >>f,line[:-1]
         if a=="H":
            for line in Hbas_lines:
               print >>f,line[:-1]
         if a=="O":
            for line in Obas_lines:
               print >>f,line[:-1]
         if a=="C":
            for line in Cbas_lines:
               print >>f,line[:-1] 
         if a=="N":
            for line in Nbas_lines:
               print >>f,line[:-1] 

      print >>f," $END"
      print >>f," $ECP"
      H_flag=0
      O_flag=0
      C_flag=0
      N_flag=0
      for a in atoms:
         if a=="H":
            if H_flag==0:
               for line in Hecp_lines:
                  print >>f,line[:-1]
               H_flag=1
            else:
               print >>f,"H-QMC GEN"
         elif a=="O":
            if O_flag==0:
               for line in Oecp_lines:
                  print >>f,line[:-1]
               O_flag=1
            else:
               print >>f,"O-QMC GEN"
         elif a=="C":
            if C_flag==0:
               for line in Cecp_lines:
                  print >>f,line[:-1]
               C_flag=1
            else:
               print >>f,"C-QMC GEN"
         elif a=="N":
            if N_flag==0:
               for line in Necp_lines:
                  print >>f,line[:-1]
               N_flag=1
            else:
               print >>f,"N-QMC GEN"
      print >>f," $END"
      f.close()

step1(fname0)
step2(fname0)
step3()
