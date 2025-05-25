# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:12:25 2020

@author: LIPO1
"""
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
import math
import scipy.integrate as integrate

inf = math.inf
exp = math.exp
sin= np.sin
h2=(6.626*10**-34)**2 #Js
pi2=3.14**2
pi=3.14
e=1.6*10**-19
eps0=8.854*10**-12
#CdSe
eps2=5.7 # 5,7
me=0.19*9.1*10**-31
mh=0.8*9.1*10**-31
Eg=2.58
esum=0.0
##ZnO
#eps2=3.7 # 5,7
#me=0.24*9.1*10**-31
#mh=0.45*9.1*10**-31
#Eg=3.44
##InSb
eps2=15.6
me=0.015*9.1*10**-31
mh=0.39*9.1*10**-31
Eg=0.24
##GaAs
#eps2=10.9
#me=0.07*9.1*10**-31
#mh=0.68*9.1*10**-31
#Eg=1.52
E=np.array([])
Ekin=np.array([])
Ekul=np.array([])
Epol=np.array([])
X=np.array([])

for r in range (30,100):
    R=r*10**-10
    E1=(((h2)/(8*R**2))*((1/me) + (1/mh)))/e
    E2= ((-1.8*e**2)/(4*pi*eps2*eps0*R))/e
    E3=(((e**2)*(eps2-1)/(4*pi*eps2*eps0*(R**2)))*integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((eps2+1)*(y+1)), [0, inf]), 0, R)[0])/e
    E=np.append(E,[Eg+E1+E2+E3])
    Ekin=np.append(Ekin,[E1])
    Ekul=np.append(Ekul,[E2])
    Epol=np.append(Epol,[E3])
    X=np.append(X,[r])

plt.figure()
plt.plot(X,E)
plt.xlabel('Promień (A)')
plt.ylabel('Energia (eV)')
plt.title('Energia całkowita')

plt.figure()
plt.plot(X,Ekin)
plt.xlabel('Promień (A)')
plt.ylabel('Energia (eV)')
plt.title('Energia kinetyczna nosnikow')

plt.figure()
plt.plot(X,Ekul)
plt.xlabel('Promień (A)')
plt.ylabel('Energia (eV)')
plt.title('Energia odpychania kulombowskiego')

plt.figure()
plt.plot(X,Epol)
plt.xlabel('Promień (A)')
plt.ylabel('Energia (eV)')
plt.title('Energia polaryzacji')

#Epol=0.0
#R=30*10**-10
#print(eint)
#eint = integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((eps2+1)*(y+1)), [0, inf]), 0, R)
#print(eint[0])
#żeby uzyskac Epol taka z korelacja elektron dziura nalezy dzielic pierwszy czlon przez 2, nie przez 4
#Epol = (((e**2)*(eps2-1)/(4*pi*eps2*eps0*(R**2)))*integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((eps2+1)*(y+1)), [0, inf]), 0, R)[0])/e
#print(Epol)