#This is a Gene Regulatory Network describing the interaction between 
#Cdx2 and Oct4
#S = 0.6 for Outer Cell
#S = -0.2 for Inner Cell

import tellurium as te
import matplotlib.pyplot as plt
import numpy as np

Model = '''
    //Reactions
    J0: -> x ; k * ( b + S + x^h/(x^h+0.5^h)) * ((1-I) + I * 0.5^h/(y^h + 0.5^h)) + x * rando ;
    J1: x -> ; x ; 
    J2: -> y ;  k * ( b + y^h/(y^h+0.5^h)) * ((1-I) + I * 0.5^h/(x^h + 0.5^h)) + y * rando ;
    J3: y -> ; y ;
    
    //Species
    x = 1 ;
    y = 0 ;
    
    //Parameters
    k = 0.7 ;
    b = 0.7 ;
    S = - 0.2 ; 
    h = 4 ;
    I = 0.6 ;
    d = 0.5 ;
    rando = 0 ;
    
'''
m = te.loada(Model)
#for k in np.arange(1.0,1.6,0.05):
    #corrections = 0
t = 140
t = np.arange(0,140)
sigma =4.0
results = [[0,1,0]]
corrections = 0
for j in range(0,100): #simulating 100 cells
    m = te.loada(Model)
    for i in range(1,len(t)):
        if i>0 and i%10==0:
            m.rando = np.random.normal(0,sigma)
        s = m.simulate(t[i-1],t[i],2)
        results.append([s[1][0],s[1][1],s[1][2]])
    if results[len(results)-1][1]<results[len(results)-1][2]:
        corrections+=1
print corrections/100.0
#results = np.array(results)
#plt.plot(t,results[:,1])
#plt.plot(t,results[:,2])
#plt.axis([0,140,0,1.5])
#plt.show()



#REPLICATING PART
#time_length = 200.0
#points = time_length*100
#sigma = 3
#t = np.arange(0,time_length,time_length/points)
#R0 = [[0, 1.0, 0.0], [0, 1.05, 0.0], [0, 1.1, 0.0], [0, 1.15, 0.0], [0, 1.2, 0.0], [0, 1.25, 0.0], [0, 1.3, 0.0], [0, 1.35, 0.0], [0, 1.4, 1.0], [0, 1.45, 1.0], [0, 1.5, 1.0], [0, 1.55, 1.0], [0, 1.6, 1.0]]
#R1 = [[1, 1.0, 0.0], [1, 1.05, 0.0], [1, 1.1, 0.0], [1, 1.15, 0.0], [1, 1.2, 0.0], [1, 1.25, 0.0], [1, 1.3, 0.08], [1, 1.35, 0.71], [1, 1.4,1.00], [1, 1.45, 1.00], [1, 1.5, 1.00], [1, 1.55, 1.00], [1, 1.6, 1.00]]
#R2 = [[2, 1.0, 0.0], [2, 1.05, 0.0], [2, 1.1, 0.02], [2, 1.15, 0.17], [2, 1.2, 0.28], [2, 1.25, 0.69], [2, 1.3, 0.93], [2, 1.35, 1.0], [2, 1.4, 1.0], [2, 1.45, 1.0], [2, 1.5, 1.0], [2, 1.55, 1.0], [2, 1.6, 1.0]]
#R3 = [[3, 1.0, 0.17], [3, 1.05, 0.35], [3, 1.1, 0.66], [3, 1.15, 0.84], [3, 1.2, 0.97], [3, 1.25, 0.98], [3, 1.3, 1.0], [3, 1.35, 1.0], [3, 1.4, 1.0], [3, 1.45, 1.0], [3, 1.5, 1.0], [3, 1.55, 1.0], [3, 1.6, 1.0]]
#Rx = [[3, 1.0, 0.0], [3, 1.05, 0.04], [3, 1.1, 0.03], [3, 1.15, 0.18], [3, 1.2, 0.42], [3, 1.25, 0.81], [3, 1.3, 0.98], [3, 1.35, 1.0], [3, 1.4, 1.0], [3, 1.45, 1.0], [3, 1.5, 1.0], [3, 1.55, 1.0], [3, 1.6, 1.0]]
#Ry = [[3, 1.0, 0.06], [3, 1.05, 0.09], [3, 1.1, 0.21], [3, 1.15, 0.38], [3, 1.2, 0.66], [3, 1.25, 0.87], [3, 1.3, 0.98], [3, 1.35, 1.0], [3, 1.4, 1.0], [3, 1.45, 1.0], [3, 1.5, 1.0], [3, 1.55, 1.0], [3, 1.6, 1.0]]
##for k in np.arange(1.0,1.6,0.05):
##    corrections = 0
##    for j in range(0,100):
##        m = te.loada(Model)
##        results = []
##        m.c = k
##        for i in range(1,len(t)):
##            m.xsigma = np.random.normal(0,0)
##            m.ysigma = np.random.normal(0,sigma)
##            s = m.simulate(t[i-1],t[i],2)
##            results.append([s[1][0],s[1][1],s[1][2]])
##        if results[len(results)-1][1]<results[len(results)-1][2]:
##            corrections+=1
##    Ry.append([sigma,m.d,corrections/100.0])
##    print Ry
#
##Plotting
#R0 = np.array(R0)
#R1 = np.array(R1)
#R2 = np.array(R2)
#R3 = np.array(R3)
#Rx = np.array(Rx)
#Ry = np.array(Ry)
#plt.plot(R0[:,1],R0[:,2])
#plt.plot(R1[:,1],R1[:,2])
#plt.plot(R2[:,1],R2[:,2])
#plt.plot(R3[:,1],R3[:,2])
#plt.plot(Rx[:,1],Rx[:,2],'--')
#plt.plot(Ry[:,1],Ry[:,2],'--')
#plt.show()
