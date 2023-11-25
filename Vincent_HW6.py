"""
Vincent Andrews 
Orbital Mechanics Homework 6
"""
import numpy as np
pi = np.pi
mu = 398600
R_e = 6378

"""
Question 1- Find the right ascension α and declination δ after 30 minutes have passed.
"""
print("QUESTION 1 RESULTS:")
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
print("initial guess: ", X)
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
    
# final stumpff functions and f, g computed in calculator 
l2 = -0.1990
m2 = 3
n2 = 0.9174
DEC_final = np.arcsin(n2) * 180/pi

cos_df = np.cos(DEC_final)
RA_final = 120

print("The final declination is: ", DEC_final)
print("The final right ascension is: ", RA_final)
"""
QUESTION 1 RESULTS:
initial guess:  45.348090892279906
The final declination is:  66.54889148795415
The final right ascension is:  120
"""

"""
Question 2 - Find the six orbital elements of a geocentric satellite whose inertial position and velocity vectors
in the geocentric equatorial frame are r and v respectively 
"""
pi = np.pi
r2 = np.array([0, 0, -13000], dtype=np.int64)
v2 = np.array([4, 5, 6], dtype=np.int64)
mu = 398600
print("\n")
print("QUESTION 2 RESULTS:")

print("The six Orbital elements are:\n")
#computing angular momentum vector and magnitude
h2 = np.cross(r2, v2)

hdoth = np.int64(0)
for x in h2:
    hdoth = np.int64(x*x + np.int64(hdoth))
h_mag2 = np.sqrt(hdoth)
print("Angular momentum = ", h_mag2)

#compute inclination 
i2 = (np.arccos(0) * 180 )/pi
print("inclination = ", i2)

#Compute N vector and magnitude 
k = np.array([0, 0, 1], dtype=np.int64)
N = np.cross(k, h2)

NdotN = np.int64(0)
for x in N:
    NdotN = np.int64(x*x + np.int64(NdotN))
N_mag = np.sqrt(NdotN)

#computer right ascension of ascending node
ra2 = (np.arccos(52000/h_mag2) * 180 )/pi
print("right ascension of ascending node = ", ra2)

#compute eccentricity vector 
#left side 
v_mag2 = np.sqrt(np.dot(v2,v2))
v_squared2 = np.square(v_mag2)
r_mag2 = np.sqrt(np.dot(r2,r2))
left_scalar = (v_squared2 - mu/r_mag2)/mu
left_vector = np.array([left_scalar * x for x in r2])

#right side 
rdotv2 = np.dot(r2,v2)
right_scalar = rdotv2 /mu
right_vector = np.array([right_scalar * x for x in v2])

e_vector2 = np.subtract(left_vector, right_vector)
e_mag2 = np.sqrt(np.dot(e_vector2, e_vector2))
print("eccentricity = ", e_mag2)

#argument of perigee 
Ndote = np.dot(N,e_vector2)
Ne = N_mag * e_mag2
w2 = 360 - (np.arccos(Ndote/Ne) * 180)/pi
print("argument of perigee = ", w2)

#Compute true anomaly 
edotr2 = np.dot(e_vector2,r2)
er = e_mag2 * r_mag2
true_anomaly2 = ((2 * pi - np.arccos(edotr2/er)) * 180)/pi
print("true anomaly = ", true_anomaly2)
"""
QUESTION 2 RESULTS:
The six Orbital elements are:

Angular momentum =  83240.61508662703
inclination =  90.0
right ascension of ascending node =  51.34019174590991
eccentricity =  1.2975693345987163
argument of perigee =  344.93852998712333
true anomaly =  285.06147001287667
"""

"""
Question 3 - Given r and e, as well as the additional information that the satellite is approaching
perigee, compute the following:
"""
#calculate true anomaly
r3 = np.array([-6600, -1300, -5200], dtype=np.int64)
e3 = np.array([-0.4, -0.5, -0.6])
edotr3 = np.dot(e3,r3)
r_mag3 = np.sqrt(np.dot(r3,r3))
e_mag3 = np.sqrt(np.dot(e3,e3))
er3 = r_mag3 * e_mag3
true_anomaly3 = ((2 * pi - np.arccos(edotr3/er3)) * 180)/pi
print("\n")
print("QUESTION 3 RESULTS:")

print("true anomaly = ", true_anomaly3)

#calculate angular momentum 
#solved for h using orbit equation
b3 = r_mag3 * mu
c = (1 + e_mag3)
d = b3 * c
h3 = np.sqrt(d)
print("angular momentum = ", h3)

#compute speed of satellite 
#first find velcoity vector 
v_r3 = (mu/h3) * e_mag3 * np.sin(true_anomaly3)
v_p3 = (mu/h3) * (1 + e_mag3 * np.cos(true_anomaly3))
speed3 = np.sqrt(np.square(v_r3) + np.square(v_p3))
print("speed  = ", speed3)
"""
QUESTION 3 RESULTS:
true anomaly =  329.2222885713851
angular momentum =  79767.82754480801
speed  =  3.029363656733253
"""

"""
Question 4 - Consider the classical Euler sequence α “ 350 ̋, β “ 170 ̋, γ “ 300 ̋.
"""
a = (350 * pi/180)
b = (170 * pi/180)
y = (300 * pi/180)
Q11 = ((-np.sin(a) * np.cos(b) * np.sin(y)) + (np.cos(a) * np.cos(y))) 
Q12 = ((np.cos(a) * np.cos(b) * np.sin(y)) + (np.sin(a) * np.cos(y))) 
Q13 = (np.sin(b) * np.sin(y)) 
Q21 = ((-np.sin(a) * np.cos(b) * np.sin(y)) - (np.cos(a) * np.cos(y))) 
Q22 = ((np.cos(a) * np.cos(b) * np.sin(y)) - (np.sin(a) * np.cos(y))) 
Q23 = (np.sin(b) * np.cos(y)) 
Q31 = (np.sin(a) * np.sin(b)) 
Q32 = (-np.cos(a) * np.sin(b)) 
Q33 = np.cos(b) 
DCM = np.array([[[Q11, Q12, Q13], [Q21,Q22,Q23], [Q31,Q32,Q33]]])
print("\n")
print("QUESTION 4 RESULTS:")

print("The direction cosdine matrix is: \n", DCM)
print("\n")

#Compute yaw,pitch,roll sequence 
yaw = (np.arctan(Q12/Q11) * 180)/pi
pitch = (np.arcsin(-Q13) * 180)/pi
roll = (np.arctan(Q23/Q33) * 180)/pi + 180
print("the yaw,pitch,roll sequence is:")
print("yaw = ", yaw)
print("pitch = ", pitch)
print("roll = ", roll)
"""
QUESTION 4 RESULTS:
The direction cosdine matrix is: 
 [[[ 0.64050294  0.75308745 -0.15038373]
  [-0.34430481  0.92673563  0.08682409]
  [-0.03015369 -0.17101007 -0.98480775]]]


the yaw,pitch,roll sequence is:
yaw =  49.61874485752944
pitch =  8.649165105287574
roll =  174.9616312267025
"""
"""
Question 5 - The hyperbolic trajectory of a spacecraft passing Earth has the following orbital parameters:
e = 1.5, perigee altitude of 3300 km, i = 35 , Ω = 130 and ω = 115 :Find r and v at perigee in...

"""
#find r and v at perigee in the perifocal frame 
e5 = 1.5
z_p5 = 3300
i5 = 35
W = 130
w5 = 115
print("\n")
print("QUESTION 5 RESULTS:")

r_p5 = R_e + z_p5

#theta = 0 at perigee
theta5 = 0 
h5 = np.sqrt(r_p5 * mu * (1+e5))
print("In the perifocal frame: ")
pf_r_unitvector = np.array([np.cos(theta5), np.sin(theta5), 0])
r_pf = np.array([r_p5 * x for x in pf_r_unitvector])

pf_v_unitvector = np.array([np.sin(theta5), e5 + np.cos(theta5), 0])
v_pf = np.array([(mu/h5) * x for x in pf_v_unitvector])
print("the r vector is: ", r_pf)
print("the v vector is: ", v_pf)
print("\n")

# calculate direction cosine matrix, P 
# i5 = 35, w5 = 115, W=130
P11 = (np.cos(W) * np.cos(w5)) - (np.sin(W) * np.sin(w5) * np.cos(i5))
P12 = (-np.cos(W) * np.sin(w5)) - (np.cos(W) * np.cos(i5) * np.cos(w5))
P13 = np.sin(W) * np.sin(i5)
P21 = (np.sin(W) * np.cos(w5)) + (np.cos(W) * np.cos(i5) * np.sin(w5))
P22 = (-np.sin(W) * np.sin(w5)) + (np.cos(W) * np.cos(i5) * np.cos(w5))
P23 = -np.cos(W) * np.sin(i5)
P31 = np.sin(i5) * np.sin(w5)
P32 = np.sin(i5) * np.cos(w5)
P33 = np.cos(i5)
P = np.array([[[P11,P12,P13], [P21,P22,P23], [P31,P32,P33]]])
print("In the geocentric frame: ")

print("the direction cosine matrix for transformation is:\n", P)

R_GEO = np.array([-6532.64196, 5969.82532, -3917.83857])
V_GEO = np.array([4.62095, 7.82565, 1.41559])
print("\n")
print("the r vector is: ", R_GEO)
print("the v vector is: ", V_GEO)
"""
QUESTION 5 RESULTS:
In the perifocal frame: 
the r vector is:  [9678.    0.    0.]
the v vector is:  [ 0.         10.14719117  0.        ]


In the geocentric frame: 
the direction cosine matrix for transformation is:
 [[[-0.67499917  0.45539244  0.39825525]
  [ 0.61684494  0.77121279 -0.15726778]
  [-0.40481903  0.13950611 -0.90369221]]]


the r vector is:  [-6532.64196  5969.82532 -3917.83857]
the v vector is:  [4.62095 7.82565 1.41559]
"""

"""
Question 6
"""
# because r, h, and e are all orthonormal to each other, they form the orbital plane 
# the normal to this plane is defined as w
print("\n")
print("QUESTION 6 RESULTS:")
print("Because r, h, and e are orthonomral with each other, the normal vector can be expressed as w(angular velocity)")

r6 = np.array([-6600, -1300, 5200])
e6 = np.array([-0.4, -0.5, -0.6])

rcrosse6 = np.array([-1820, -1880, 2780])
rcrosse6_mag = np.sqrt(np.dot(rcrosse6, rcrosse6))

w6 = np.array([(1/rcrosse6_mag) * x for x in rcrosse6])
print("w = ", w6)
i6 = np.arccos(0.728178) * 180/pi
print("the inclination angle is: ", i6)

"""
QUESTION 6 RESULTS:
Because r, h, and e are orthonomral with each other, the normal vector can be expressed as w(angular velocity)
w =  [-0.47672083 -0.4924369   0.72817797]
the inclination angle is:  43.26613394608721
"""


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
