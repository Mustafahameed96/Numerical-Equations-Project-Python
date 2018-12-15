# -*- coding: utf-8 -*-
"""
Created on Tue Apr 17 16:21:43 2018

@author: mmustafa

"""

#Mini Project Numerical Computing
#Name: Muhammad Mustafa Hameed
#Roll NO: CS151020
#Section: CS-6A
#Subject: Numerical Computing
#Notes: Change if condition on line 98 as per wish to see different results on basis of h and k max 8




from numpy.linalg import inv
import numpy as np
import matplotlib.pyplot as plt


print("Multi-Dimensional Unconstrained Optimization Problem")

#complete list for h and k varying from 0.1 to 1.0
newlist=[0.1,0.5,1.0]



j=0
while True:

#declaring h and k
    h=newlist[j]
    k=newlist[j]
    print ("\n******************************")
    print ("H=",h,"      &&      K =",k)
    print ("******************************\n\n")    
#Declaring initial variables
    x = 0.0 
    y = (1+5+1+0+2+0)
    startMatrix= np.matrix([[x],[y]])
    xlist = np.linspace(-2.75, 1.5, 100)
    ylist = np.linspace(-2.75, 1.5, 100)
    X, Y = np.meshgrid(xlist, ylist)
    Z = np.sqrt(X**2 + Y**2)
    plt.figure()
    cp = plt.contourf(X, Y, Z)
    plt.colorbar(cp)
    plt.title('Filled Contoured Plot of convergence')
    myListX=[]
    myListY=[]
    plt.xlabel('x (cm)')
    plt.ylabel('y (cm)')
   
    
    i = 0
    while True:
                  
            print ("Iteration Number :",i,"\n")
            print("Starting matrix = \n",startMatrix,"\n")     
            
#First Derivative i.e function'x
            
            var1 = ( ((x+h)**3) + (2*(x+h)*y) + (2*(y)**2) - (2*(x+h)) + (y) + 8 )
            var2 = ( ((x-h)**3) + (2*(x-h)*y) + (2*(y)**2) - (2*(x-h)) + (y) + 8 )
            
            fx = (var1 - var2)/(2 * h)
            
#First Derivative i.e function'y
            var1 = ( ((x)**3) + (2*x*(y+k)) + (2*(y+k)**2)  - (2*(x)) + (y+k) + 8)
            var2 = ( ((x)**3) + (2*x*(y-k)) + (2*(y-k)**2)  - (2*(x)) + (y-k) + 8)
            
            fy = (var1 - var2)/(2 * k)
            
#Second Derivative i.e function"x
            var1 = ( ((x+h)**3) + (2*(x+h)*y) + (2*(y)**2) - (2*(x+h)) + (y) + 8)
            var2 = ( ((x)**3)   + (2*(x)*y)   + (2*(y)**2) - (2*(x))   + (y) + 8)
            var3 = ( ((x-h)**3) + (2*(x-h)*y) + (2*(y)**2) - (2*(x-h)) + (y) + 8)
            
            fxx = (var1 - (2 * var2) + var3)/((h)**2)
            
#Second Derivative i.e function"y
            var1 = ( ((x)**3) + (2*x*(y+k)) + (2*(y+k)**2) - (2*(x)) + (y+k) + 8)
            var2 = ( ((x)**3) + (2*(x)*y)   + (2*(y)**2)   - (2*(x)) + (y)   + 8)
            var3 = ( ((x)**3) + (2*x*(y-k)) + (2*(y-k)**2) - (2*(x)) + (y-k) + 8)
            
            fyy = (var1 - (2 * var2) + var3)/((k)**2)
            
#Second Derivative i.e function"xy
            var1 = ( ((x+h)**3)  + (2*(x+h)*(y+k)) + (2*(y+k)**2) - (2*(x+h)) + (y+k) + 8)
            var2 = ( ((x+h)**3)  + (2*(x+h)*(y-k)) + (2*(y-k)**2) - (2*(x+h)) + (y-k) + 8)            
            var3 = ( ((x-h)**3)  + (2*(x-h)*(y+k)) + (2*(y+k)**2) - (2*(x-h)) + (y+k) + 8)         
            var4 = ( ((x-h)**3)  + (2*(x-h)*(y-k)) + (2*(y-k)**2) - (2*(x-h)) + (y-k) + 8)
   
            fxy = (var1 - var2 - var3 + var4)/(4*h*k)
          
            var1 = fxx
            var2 = fxy
            var3 = fxy
            var4 = fyy
            vara = fx
            varb = fy
            
#   
            
            naabla= np.array([[vara],[ varb]])
    
#            print("delta F Nable = \n",naabla,"\n")
        
            Hessian= np.array([[var1,var2 ], [var3, var4]])
#            print("Hessian = \n",Hessian,"\n")
            
            Hinv = inv(Hessian)
#            print("Hessian inverse = \n",Hinv,"\n")
   
            startMatrix=startMatrix-(np.matmul(Hinv,naabla))
            
            print("Xi+1 = Xi - (Hessian Inverse * Naabla) \n")
            print(startMatrix,"\n")
            
            x=startMatrix[0,0]
            x=float(format(x,".5f"))
        
            y=startMatrix[1,0]
            myListX.append([x])
            myListY.append([y])
                  
#breaking 2nd Loop            
            if i >5:
                break
            i = i + 1
    plt.scatter(myListX,myListY)
    plt.show()             
#breaking 1st loop            
    if j >1:
        break
    j = j + 1
