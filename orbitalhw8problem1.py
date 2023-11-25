
import numpy as np
pi = np.pi
mu = 398600
R_e = 6378

"""
Question 1- Find the right ascension α and declination δ after 30 minutes have passed.
"""
z = 500
r = z + R_e
RA_0 = 120
DEC_0 = -20
dt = 1800

#GEF trajectory 
l = (np.cos(RA_0) * 180/pi) * (np.cos(DEC_0) * 180/pi)
m = (np.sin(RA_0) * 180/pi) * (np.cos(DEC_0) * 180/pi)
n = np.sin(DEC_0) * 180/pi

u_r = np.array([l, m, n])
r_0 = np.array([r * x for x in u_r])
v_0 = np.array([0, 0, 10])

v_0_mag = np.sqrt(np.dot(v_0, v_0))
h1 = np.cross(r_0,v_0)
h1_mag = np.sqrt(np.dot(h1,h1))

Vr_0 = (np.dot(r_0, v_0))/r

#apply conservation of energy 
alpha = 2/r - np.square(v_0_mag)/mu
X = np.sqrt(mu) * abs(alpha) * dt
z = alpha * np.square(X)
print("initiasl guess: ", X)
print("\n")

#need to find stumpff functions 
if z > 0:
    s = (np.sqrt(z) - (np.sin(np.sqrt(z)) * 180/pi))/np.power(np.sqrt(z), 3)
elif alpha < 0:
    s = (np.sinh(np.sqrt(-1 * z)) - np.sqrt(-1 * z))/np.power(np.sqrt(-1 * z), 3)
else:
    s = 1/6

if z > 0:
    c = (1 - (np.cos(np.sqrt(z)) * 180/pi))/z
elif z < 0:
    c = (np.cosh(np.sqrt(-1 * z)) - 1)/(-1 * z)
else:
    c = 1/2

#Solving kepler's equation numerically
err = 1e-8
Max_n = 100

n=0
ratio=1
while abs(ratio) > err and n <= Max_n:
    n = n + 1
    F = r*Vr_0/np.sqrt(mu)*(np.square(X)*c) + ((1-alpha*r)*np.power(X, 3)*s + r*X - np.sqrt(mu)*dt)
    dFdX = (r*Vr_0/np.sqrt(mu))*X*(1 - alpha*np.square(X)*s) + ((1 - alpha*r)*np.square(X)*c + r)
    ratio = F/dFdX
    X = X - ratio
    print(X)

print("\n")

# final stumpff functions and f, g computed in calculator 
l2 = 1
m2 = 3
n2 = 2

DEC_final = np.sin(n2) * 180/pi

if m2 > 0:
    RA_final = np.cos(11/np.cos(DEC_final)) * 180/pi
else:
    RA_final = 360 - np.cos(11/np.cos(DEC_final)) * 180/pi
print("The final declination is: ", DEC_final)
print("The final right ascension is: ", RA_final)
print("\n")

"""
z_f = X * np.square(alpha)
c_f = (1 - (np.cos(np.sqrt(z_f)) * 180/pi)/z_f

       
#f and g components using X
f = (1 - np.square(X)/r )* c_f
g = dt - (1/np.sqrt(mu))*np.power(X,3)*s_f

rfg1 = np.array([f * x for x in r_0])
rfg2 = np.array([g * x for x in v_0])
r_fg = np.add(rfg1, rfg2)
r_fg_mag = np.sqrt(np.dot(r_fg, r_fg))

finding later RA and DEC
l2 = r_fg[0]/r_fg_mag
m2 = r_fg[1]/r_fg_mag
"""