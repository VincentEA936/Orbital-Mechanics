"""
Vincent Andrews
Homework 9
"""
import numpy as np
import math
pi = math.pi 
mu_e = 398600
mu_sun = 132.71*pow(10,9)
R_etosun = 149.6*pow(10,6)
R_mtosun = 227.9*pow(10,6)
"""
Question 1 - What is the total ∆v required for a Hohmann transfer between Earth’s
orbit and Venus’ orbit? Ignore the velocity changes within each planet’s SOI, 
and focus just on the heliocentric portion of the transfer
"""
# only need to consider the delta v required for departure and arrival. 
# at Earth's SOI, v = v_D
# v_D = v_p - V_e
R_vtosun = 108.2*pow(10,6)
h = np.sqrt(2*mu_sun*((R_etosun*R_vtosun)/(R_etosun + R_vtosun)))
v_p = h/R_etosun

v_e  = np.sqrt(mu_sun/R_etosun)
v_D = v_p - v_e

# compute arrival speed into Venus orbit 
v_A = np.sqrt(mu_sun/R_vtosun)*(1 - np.sqrt((2*R_etosun)/(R_vtosun + R_etosun)))

delta_v = v_A + v_D
"""
Question 1 Results:
    The delta v for an Earth to Venus transfer focusing only on the helicoentric 
    velocities is -5.2035 km/s
"""

"""
Question 2 - Compute the synodic period of Jupiter relative to Saturn
"""
T_J = 11.86
T_S = 29.46
T_syn = (T_J*T_S)/abs(T_J-T_S)
print("The launch window opens every ", T_syn,"years ")
"""
Question 2 Results:
    The launch window from Jupiter to Saturn opens once every T_syn years which
    is Y_syn = 19.852 years 
"""

"""
Question 3 - Compute the radius of Trappist 1b's sphere of influence
"""
M_sun = 1.988*pow(10,30)
M_earth = 5.974*pow(10,24)
M_star = 0.0898*M_sun
M_trap1b = 1.374*M_earth
d = 1.726*pow(10,6)

mass_ratio = M_trap1b/M_star
R_SOI_trap1b = d*pow(mass_ratio,2/5)
"""
Question 3 Results:
    The radius of the sphere of influence of Trappist-1b is 31,773.458 km
"""

"""
Question 4 - Consider a mission from Earth to Mars. The Earth mission starts in a circular parking orbit of
radius (not altitude) 7000 km.
"""
# a) compute the delta v required as measured in the heliocentric frame 
#step 1 - find the angular momentum of the Hohmann ellipse 

h = np.sqrt(2*mu_sun*((R_etosun*R_mtosun)/(R_etosun + R_mtosun)))

#step 2 - compute Earth's orbital speed
v_e = np.sqrt(mu_e/R_etosun)
#step 3 - compute the perihelion velocity 
v_ph = h/R_etosun
#step 4 - calculate V_infinity 
v_infin = v_ph - v_e


# b)  compute the total delta v from circular orbit to the hyperbolic trajectory 
# matching the Hohmann transfer ellipse 
# v_p_hyp^2 = v_infin^2 + v_esc^2 
#solve for v_esc 
r_park = 7000 # parking orbit radius 
v_esc = np.sqrt((2*mu_sun)/r_park)
vp_hyp = np.sqrt(v_infin**2 + v_esc**2)

v_i = np.sqrt(mu_e/r_park) # initial velocity in the  circular parking orbit 
delta_v = vp_hyp - v_i
"""
Question 4 Results:
    part a) v = 2.943 km/s
    part b) delta_v = 3.5240 km/s
"""