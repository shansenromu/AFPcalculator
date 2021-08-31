#!/usr/bin/python

# Simulates AFP-NMR sequence by integrating the Bloch equations
# by Jeff

# Editied by Sean to calcualte AFP of neutrons in TCN18-180 using course hand measured feilds and calucuated spinflipper coil. 

import scipy
print(scipy.__version__)

import sys

from math import *
from math import sqrt
import numpy as np
from scipy.integrate import ode
import scipy.constants as c

import matplotlib.pyplot as plt
from matplotlib.pyplot import figure



from b_grid import bField
from b_grid import b1Field
from b_grid import ReadB
from b_grid import ReadB1

print("Reading files")
X0D,Y0D,Z0D, B0xField, B0yField, B0zField = ReadB("CourseFieldMap-SCM18-070.txt", 13,2,2)
X1D,Y1D,Z1D, B1xField, B1yField, B1zField = ReadB1("B1grid.tsv", 44,44,401)
print("Files read")


# python afp_UCN18-180.py [Y value in cm] [Offset in the SF position]
Yvalue_arg = sys.argv[1] 
SF_Offset_arg = sys.argv[2]

# global parameters of the calculation
# constants of nature
pi=c.pi             # pi
hbar=c.hbar         # hbar
k=c.k               # Boltzmann's constant (SI)
go2pi = c.physical_constants["neutron gyromagnetic ratio over 2 pi"][0]

                    # gamma/(2*pi) in Hz/uT
gamma=go2pi*2*pi    # gamma (radians/s)/uT
print(gamma)
mu_p=c.physical_constants['neutron magnetic moment'][0]
                    # neutron magnetic moment (J/T)

# Boltzmann polarization
T=295               # Kelvin
c=mu_p/(k*T)/1e6    # (1/uT) Boltzmann polarization is Bz times this.

print(c)

# relaxation times
T1=100.0     # seconds
#T2=100.0        # seconds

#neutron stuff
velocity = -800 #cm/s
x0 = 130 #cm
xf =  10 #cm

# afp parameters
freq=7.4e3         # b1 frequency (Hz)
omega=freq*2*pi     # b1 frequency (rad/s)
centralb0=omega/gamma
                    # central b0 field (uT)
b1=100#0.1*centralb0 # b1 amplitude (uT)     #10 uT was good
b1scale = 1.16 #B1 gird was calculated with 1 A but actual Amps was 1.16 A
#deltab=b1*10.       # range to sweep (uT)
#deltat=2.           # sweep time (seconds)
#b00=centralb0-deltab/2
                    # starting b0
#print("starting b = {0}".format(b00))
#db0dt=deltab/deltat # b0 ramp rate (uT/s)


t0= 0.               # starting time
#t1=(xf-x0)/velocity # ending time (seconds)
t1 = 0.1499
dt=.00001           # time step (seconds)

#Starting position
#yAxis = -4.0
yAxis = float(Yvalue_arg)

#Center point of SF is  cm away from the 1st flange of the analyzer foil
SFcentmin = 28.1327 #cm
#SFoffset = 20 #cm
SFoffset = float(SF_Offset_arg)
SFcent = SFcentmin + SFoffset
print(SFcent)
# lock-in parameters
tau=.02               # lock-in time constant (seconds)

#Calculate the unit vector in the direction of B_total(t=0)
B0start = bField(X0D,Y0D,Z0D, B0xField, B0yField, B0zField,   x0,yAxis,0,  13,2,2 ) 
B1start = b1Field(X1D,Y1D,Z1D, B1xField, B1yField, B1zField,  x0-SFcent,yAxis,0, 44,44,401)  
B1scale = [b1scale*B1start[0],b1scale*B1start[1],b1scale*B1start[2]]
vec = B0start+B1scale
coeff_vec = 1/sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2])
unitvec =[coeff_vec*vec[0],coeff_vec*vec[1],coeff_vec*vec[2]] 
s0=unitvec      # starting polarization vector

# sparsification factor
sparse=1000000
i=0

def b(t):
    B0calc = bField(X0D,Y0D,Z0D, B0xField, B0yField, B0zField, x0+velocity*t,yAxis,0,  13,2,2 )  #bField("CourseFieldMap-SCM18-070.txt", x0+velocity*t,-5,0, 13,2,2)
    B1calc = b1Field(X1D,Y1D,Z1D, B1xField, B1yField, B1zField,x0+velocity*t-SFcent,yAxis,0, 44,44,401)  #b1Field("B1grid.tsv", 0,-5.0,x0+velocity*t-30,44,44,401 )
    # need to scale B1 by 1.16
    bx0 = B0calc[0]
    by0 = B0calc[1]
    bz0 = B0calc[2]
    # trying to fix the mismatch of coordinate systems
    bx1_m = b1scale*B1calc[2]
    by1_m = b1scale*B1calc[1]
    bz1_m = b1scale*B1calc[0]

    bx1 = cos(omega*t)*bx1_m
    by1 = cos(omega*t)*by1_m
    bz1 = cos(omega*t)*bz1_m



    bx= bx0 + bx1
    by= by0 + by1 
    bz= bz0 + bz1
    
  

    return [bx,by,bz,bx0,by0,bz0,bx1,by1,bz1,bx1_m,by1_m,bz1_m]  # uT

def dsdt(t,s):
    # The Bloch equations in the non-rotating frame
    # f = dS/dt = gamma S x B + relaxation
    sx,sy,sz=s
    bx,by,bz,bx0,by0,bz0,bx1,by1,bz1,bx1_m,by1_m,bz1_m=b(t)
    dsxdt=gamma*(sy*bz-sz*by)# -(sx-c*bx)/T2
    dsydt=gamma*(sz*bx-sx*bz)# -(sy-c*by)/T2
    dszdt=gamma*(sx*by-sy*bx)# -(sz-c*bz)/T1
    return [dsxdt,dsydt,dszdt]

r=ode(dsdt).set_integrator('dop853',atol=s0[2]*1e-8)
# Seems necessary to crank down atol for small initial polarizations?
r.set_initial_value(s0,t0)

#s_demod_x=np.array([0,0,0]) # components of polarization demodulated with cos
#s_demod_y=np.array([0,0,0]) # components of polarization demodulated with sin
sol=[]
#sol_demod_x=[]
#sol_demod_y=[]
while r.successful() and r.t < t1:
    r.integrate(r.t+dt) # do an integration step
    t=r.t               # time after the step
    s=r.y               # polarization vector after the step
    sx,sy,sz=s          # polarization vector after the step
    s=np.array(s)
    #s_demod_x=s_demod_x+(-s_demod_x+s*cos(omega*t))*dt/tau
    #s_demod_y=s_demod_y+(-s_demod_y+s*sin(omega*t))*dt/tau
    i+=1
  
    
    if(i%sparse):
        
        bxdt,bydt,bzdt,bx0dt,by0dt,bz0dt,bx1dt,by1dt,bz1dt,bx1_mdt,by1_mdt,bz1_mdt= b(t+10*dt)
        bx,by,bz,bx0,by0,bz0,bx1,by1,bz1,bx1_m,by1_m,bz1_m= b(t)
        b_coeff=1.0/(sqrt(bx*bx+by*by+bz*bz))
        b0_coeff=1.0/(sqrt(bx0*bx0+by0*by0+bz0*bz0))
        dx = dt/velocity
        #k =(gamma*( bx1_m*bx1_m + by1_m*by1_m + bz1_m*bz1_m ))/((sqrt(bx0dt*bx0dt+by0dt*by0dt+bz0dt*bz0dt)-sqrt(bx0*bx0+by0*by0+bz0*bz0))/(dt))
        #print ( gamma*(bx1_m*bx1_m + by1_m*by1_m + bz1_m*bz1_m),(sqrt(bx0dt*bx0dt+by0dt*by0dt+bz0dt*bz0dt)-sqrt(bx0*bx0+by0*by0+bz0*bz0))/(10*dt))#,k )



        if((sqrt(bx0dt*bx0dt+by0dt*by0dt+bz0dt*bz0dt)-sqrt(bx0*bx0+by0*by0+bz0*bz0))/(10*dt) == 0.0):
            k = 500.0
            k0 = 500.0
        else:
            k =abs((gamma*( bx1_m*bx1_m + by1_m*by1_m + bz1_m*bz1_m ))/((sqrt(bx0dt*bx0dt+by0dt*by0dt+bz0dt*bz0dt)-sqrt(bx0*bx0+by0*by0+bz0*bz0))/(10*dt)))
            k0 =abs((gamma*( bx0*bx0 + by0*by0 + bz0*bz0 ))/((sqrt(bx0dt*bx0dt+by0dt*by0dt+bz0dt*bz0dt)-sqrt(bx0*bx0+by0*by0+bz0*bz0))/(10*dt)))

        sol.append([t,sx,sy,sz,bx,by,bz,sx-b_coeff*bx,sy-b_coeff*by,sz-b_coeff*bz,sqrt(bx*bx+by*by+bz*bz),b_coeff*(sx*bx+sy*by+sz*bz),b0_coeff*(sx*bx0+sy*by0+sz*bz0),k,centralb0,k0])
        
        #sol_demod_x.append([s_demod_x[0],s_demod_x[1],s_demod_x[2]])
        #sol_demod_y.append([s_demod_y[0],s_demod_y[1],s_demod_y[2]])
        PerDone = t/0.1499 *100.00 
        print ("%.2f"% round(PerDone,2))
        i=0

sol=np.array(sol)

figure(figsize=(15,15),dpi=120)

plt.plot(sol[:,0],sol[:,1],color="red",label="Sx")
plt.plot(sol[:,0],sol[:,2],color="green",label="Sy")
plt.plot(sol[:,0],sol[:,3],color="blue",label="Sz")
plt.plot(sol[:,0],sol[:,11],color="black",label="S*B/|B|")
plt.plot(sol[:,0],sol[:,12],color="goldenrod",label="S*B0/|B0|")
#plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('Spin components and demodulation')

save1 = "SpinCompnents-AlongY{AlongY:.0f}cmAt{SFOffset:.0f}cmSFOffsent.png"
plt.savefig(save1.format(AlongY = yAxis, SFOffset = SFoffset))
plt.close()

figure(figsize=(15,15),dpi=120)

plt.plot(sol[:,0],sol[:,4],color="red",label="Bx")
plt.plot(sol[:,0],sol[:,5],color="green",label="By")
plt.plot(sol[:,0],sol[:,6],color="blue",label="Bz")
plt.plot(sol[:,0],sol[:,10],color="black",label="B tot")
#plt.plot(sol[:,0],sol[:,13],color="goldenrod",label="k")
plt.plot(sol[:,0],sol[:,14],color="hotpink",label="Central B0")
#plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('B ($\mu$T)')
#plt.show()

save2 = "MagneticFieldComponents-AlongY{AlongY:.0f}cmAt{SFOffset:.0f}cmSFOffsent.png"
plt.savefig(save2.format(AlongY = yAxis, SFOffset = SFoffset))
plt.close()

figure(figsize=(15,15),dpi=120)
plt.plot(sol[:,0],sol[:,13],color="goldenrod",label="k afp")
plt.plot(sol[:,0],sol[:,15],color="hotpink",label="k asp")
plt.xlabel('Time (s)')
plt.ylabel('K values')
plt.yscale("log")
save3 = "AdiabaticConditions-AlongY{AlongY:.0f}cmAt{SFOffset:.0f}cmSFOffsent.png"
plt.savefig(save3.format(AlongY = yAxis, SFOffset = SFoffset))
plt.close()
