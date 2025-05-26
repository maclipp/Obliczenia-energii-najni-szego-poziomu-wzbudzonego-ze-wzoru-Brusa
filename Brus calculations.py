# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 15:12:25 2024

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
#phsycial constants
h2=(6.626*10**-34)**2 #Js
pi2=3.14**2
pi=3.14
e=1.6*10**-19
eps0=8.854*10**-12

class MaterialParameters:
    def __init__(self, eps2, me, mh, Eg):
        self.eps2=eps2
        self.me=me
        self.mh=mh
        self.Eg=Eg
    
CdSe=MaterialParameters(5.7, 0.19*9.1*10**-31, 0.8*9.1*10**-31, 2.58)
ZnO=MaterialParameters(3.7, 0.24*9.1*10**-31, 0.45*9.1*10**-31, 3.44)
InSb=MaterialParameters(15.6, 0.015*9.1*10**-31, 0.39*9.1*10**-31, 0.24)
GaAs=MaterialParameters(10.9, 0.07*9.1*10**-31, 0.68*9.1*10**-31, 1.52)

#CdSe
# eps2=5.7 # 5,7
# me=0.19*9.1*10**-31
# mh=0.8*9.1*10**-31
# Eg=2.58

##ZnO
#eps2=3.7 # 5,7
#me=0.24*9.1*10**-31
#mh=0.45*9.1*10**-31
#Eg=3.44

##InSb
# eps2=15.6
# me=0.015*9.1*10**-31
# mh=0.39*9.1*10**-31
# Eg=0.24

##GaAs
#eps2=10.9
#me=0.07*9.1*10**-31
#mh=0.68*9.1*10**-31
#Eg=1.52

esum=0.0
E=np.array([])
Ekin=np.array([])
Ekul=np.array([])
Epol=np.array([])
X=np.array([])

Material_Input=input("Choose the material to calculate the energy parameters \n List: CdSe, ZnO, InSb, GaAs \n Selected material: ")

if Material_Input=="CdSe":
    for r in range (30,100):
        R=r*10**-10
        E1=(((h2)/(8*R**2))*((1/CdSe.me) + (1/CdSe.mh)))/e
        E2= ((-1.8*e**2)/(4*pi*CdSe.eps2*eps0*R))/e
        E3=(((e**2)*(CdSe.eps2-1)/(4*pi*CdSe.eps2*eps0*(R**2)))*integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((CdSe.eps2+1)*(y+1)), [0, inf]), 0, R)[0])/e
        E=np.append(E,[CdSe.Eg+E1+E2+E3])
        Ekin=np.append(Ekin,[E1])
        Ekul=np.append(Ekul,[E2])
        Epol=np.append(Epol,[E3])
        X=np.append(X,[r])
        
if Material_Input=="ZnO":
    for r in range (30,100):
        R=r*10**-10
        E1=(((h2)/(8*R**2))*((1/ZnO.me) + (1/ZnO.mh)))/e
        E2= ((-1.8*e**2)/(4*pi*ZnO.eps2*eps0*R))/e
        E3=(((e**2)*(ZnO.eps2-1)/(4*pi*ZnO.eps2*eps0*(R**2)))*integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((ZnO.eps2+1)*(y+1)), [0, inf]), 0, R)[0])/e
        E=np.append(E,[ZnO.Eg+E1+E2+E3])
        Ekin=np.append(Ekin,[E1])
        Ekul=np.append(Ekul,[E2])
        Epol=np.append(Epol,[E3])
        X=np.append(X,[r])
        
if Material_Input=="InSb":
    for r in range (30,100):
        R=r*10**-10
        E1=(((h2)/(8*R**2))*((1/InSb.me) + (1/InSb.mh)))/e
        E2= ((-1.8*e**2)/(4*pi*ZnO.eps2*eps0*R))/e
        E3=(((e**2)*(InSb.eps2-1)/(4*pi*InSb.eps2*eps0*(R**2)))*integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((InSb.eps2+1)*(y+1)), [0, inf]), 0, R)[0])/e
        E=np.append(E,[InSb.Eg+E1+E2+E3])
        Ekin=np.append(Ekin,[E1])
        Ekul=np.append(Ekul,[E2])
        Epol=np.append(Epol,[E3])
        X=np.append(X,[r])
        
if Material_Input=="GaAs":
    for r in range (30,100):
        R=r*10**-10
        E1=(((h2)/(8*R**2))*((1/GaAs.me) + (1/GaAs.mh)))/e
        E2= ((-1.8*e**2)/(4*pi*GaAs.eps2*eps0*R))/e
        E3=(((e**2)*(GaAs.eps2-1)/(4*pi*GaAs.eps2*eps0*(R**2)))*integrate.quad(lambda x: (sin(pi*x/R)**2)*mp.nsum(lambda y: ((x/R)**(2*y))*(y+1)/((GaAs.eps2+1)*(y+1)), [0, inf]), 0, R)[0])/e
        E=np.append(E,[GaAs.Eg+E1+E2+E3])
        Ekin=np.append(Ekin,[E1])
        Ekul=np.append(Ekul,[E2])
        Epol=np.append(Epol,[E3])
        X=np.append(X,[r])

plt.figure()
plt.plot(X,E)
plt.xlabel('Radius (A)')
plt.ylabel('Energy (eV)')
plt.title('Total energy')

plt.figure()
plt.plot(X,Ekin)
plt.xlabel('Radius (A)')
plt.ylabel('Energy (eV)')
plt.title('Carriers kinetic energy')

plt.figure()
plt.plot(X,Ekul)
plt.xlabel('Radius (A)')
plt.ylabel('Energy (eV)')
plt.title('Coulomb repulsion energy')

plt.figure()
plt.plot(X,Epol)
plt.xlabel('Radius (A)')
plt.ylabel('Energy (eV)')
plt.title('Polarization energy')
