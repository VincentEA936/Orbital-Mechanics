# -*- coding: utf-8 -*-
"""
Vincent Andrews
Orbital Mechanics - Final Exam
"""
"""
functions 
"""
def calculate_stumpff_c(z):
    
    if z>0:
        c = (1 - np.sqrt(z))/z
    elif z < 0:
        c = (np.cosh(np.sqrt(-z)) - 1)/(-z)
    else:
        c = 1/2
    return c
def calculate_stumpff_s(z):
    if z>0:
        s = (np.sqrt(z) - np.sin(np.sqrt(z)))/pow((np.sqrt(z)),3)
    elif z<0:
        s = (np.sinh(np.sqrt(-z)) - np.sqrt(-z))/pow((np.sqrt(-z)),3)
    else:
        s = 1/6
    return s

"""
Question 4 - solve for the initial velocity of the object at first sight by 
solving lamberts problem. Then compute the 6 orbital elements of the orbit
"""
import numpy as np
import math
pi = math.pi 
mu_e = 398600

"""
Part a- use lambert's problem to solve for the initial velocity of the satellite
"""
#define position vectors and delta t
r1 = np.array([-29250, -32000, 16540])
r2 = np.array([-11040, -42280, 15620])
delta_t = 7200 # seconds

#step 1: find magnitudes of r vectors
r1_mag = np.linalg.norm(r1)
r2_mag = np.linalg.norm(r2)

#step 2: find delta_theta 
r1dotr2 = np.dot(r1,r2)
delta_theta_rad = np.arccos(r1dotr2/(r1_mag*r2_mag))
delta_theta_deg = np.rad2deg(delta_theta_rad)

#step 3: find A
A = np.sin(delta_theta_deg)*np.sqrt((r1_mag*r2_mag)/(1 - np.cos(delta_theta_deg)))

#step 4: solve for z to evaluate functions 
z = 0.2066 # solved by plotting F and F prime on calculator and finding intercept

#step 5: calculate Stumpff functions 

C = calculate_stumpff_c(z)
S = calculate_stumpff_s(z)

#step 6: calculate y
y = r1_mag + r2_mag + A * ((z*S - 1)/np.sqrt(C))

#step 7: compute lagrange coefficients 
kai = np.sqrt(y/C)
f = 1 - (kai**2/r1_mag) * C
g = delta_t - (1/np.sqrt(mu_e)) * pow(kai,3) * S

#step 8: compute v1
v1 = (1/g) * (r2 - (f*r1))

"""
compute the orbital elements 
"""
# compute angular momentum 
h = np.cross(r1, v1)
h_mag = np.linalg.norm(h)
print("Angular momentum = ", h_mag)

#compute inclination 
i = (np.arccos(0) * 180 )/pi
print("inclination = ", i)

#Compute N vector and magnitude 
k = np.array([0, 0, 1], dtype=np.int64)
N = np.cross(k, h)
N_mag = np.linalg.norm(N)

#computer right ascension of ascending node
ra = (np.arccos(52000/h_mag) * 180 )/pi
print("right ascension of ascending node = ", ra)

#compute eccentricity vector 
#left side 
v_mag = np.sqrt(np.dot(v1,v1))
v_squared = v_mag**2
r_mag = np.linalg.norm(r1)
left_scalar = (v_squared - mu_e/r_mag)/mu_e
left_vector = np.array([left_scalar * x for x in r2])

#right side 
rdotv = np.dot(r1,v1)
right_scalar = rdotv /mu_e
right_vector = np.array([right_scalar * x for x in v1])

e_vector = np.subtract(left_vector, right_vector)
e_mag = np.linalg.norm(e_vector)
print("eccentricity = ", e_mag)


