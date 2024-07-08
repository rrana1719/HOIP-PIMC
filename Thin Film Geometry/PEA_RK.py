import numba, sys
from numpy import random
import scipy as scipy
import numpy as np
import hankel 
from hankel import HankelTransform 
from scipy.special import kn
from scipy import interpolate
from array import array
from statistics import stdev
from statistics import mean

#Here, we Fourier transform the quantum-confined Rytova-Keldysh potential into r-space
#We interpolate the potential and its derivative for use in the PIMC code
def ew(r,u):
    v = u*Lz/np.sqrt(r**2 + rc**2)
    num = 1-q*(1+v) + np.exp(v)*(q+v-1)
    den = (np.exp(v)-q)*(v**2)*(v**2 + 4*(np.pi**2))**2
    return (num/den)
 
def dumEw(r):
    V = -1/(Lz*ep)
    t = np.sqrt(r**2 + rc**2)
    arg = 2*np.pi*t/Lz
    func = lambda k: ew(r,k)
    nas = 3*kn(0,arg) + arg*kn(1,arg) + 32*(1/t)*(np.pi**4)*Lz*ht.integrate(func)[0]
    return( V*(nas))




ep= 6.1
ep_o = 3.3
rc = 0.898
q = (ep-ep_o)/(ep+ep_o)

s = 1*6.39*1.8897259886
Lz = 1*s

ht = HankelTransform(
    nu= 0,     # The order of the bessel function
    N = 15,   # Number of steps in the integration
    h = 0.008   # Proxy for "size" of steps in integration
)


r_s  = np.linspace(0,500,10000,endpoint=False)

sums = []

for i in range(len(r_s)):
    sums.append(dumEw(r_s[i])) 
    
    
tck = interpolate.splrep(r_s,sums)
der = interpolate.splev(r_s,tck,der=1)
Est_ar = []
ynew = interpolate.splev(r_s,tck,der=0)
for i in range(len(r_s)):
    Est_ar.append(ynew[i] + 0.5*r_s[i]*der[i])
   

#We randomly initialize the electron and hole ring polymers
@numba.jit(nopython=True)
def init(x,h):

    for i in range(0,N):
        x[i,0] = ((Lx)/(np.random.uniform(1,2*N)))
        x[i,1] = ((Ly)/(np.random.uniform(1,2*N)))
        h[i,0] = ((Lx)/(np.random.uniform(1,2*N)))
        h[i,1] = ((Ly)/(np.random.uniform(1,2*N)))
    
    return x,h

#This returns the quantum-confined Rytova-Keldysh potential at a given value of r
@numba.jit(nopython=True)
def SN(xn0,yn0,hn0,hn1):
    xu = (xn0-hn0) - Lx*np.round((xn0-hn0)/Lx)
    yu = (yn0-hn1) - Ly*np.round((yn0-hn1)/Ly)
    r = np.sqrt(xu**2 + yu**2)
    return( np.interp(r,r_s, potSums)    )
    
#This is the Virial estimator for the quantum-confined Rytova-Keldysh potential    
@numba.jit(nopython=True)
def Est(xn0,yn0,hn0,hn1):
    xu = (xn0-hn0) - Lx*np.round((xn0-hn0)/Lx)
    yu = (yn0-hn1) - Ly*np.round((yn0-hn1)/Ly)
    r = np.sqrt(xu**2 + yu**2)
    return( np.interp(r,r_s, derivSums)   )

#This computes the self-attractive piece from the charge-phonon interaction potential
@numba.jit(nopython=True)
def polPot_self(r):
    rnorm = r
    V = (2*ep_o/ref)*rnorm
    arg = (ep_o/ref)*rnorm
    if (arg < 3):
        S = 1.909859164*(arg/3) -1.909855001*(arg/3)**3 +0.687514637*(arg/3)**5 -0.126164557*(arg/3)**7 +0.013828813*(arg/3)**9 -0.000876918*(arg/3)**11
        J = 0.9999999 -2.24999239*(arg/3)**2 +1.26553572*(arg/3)**4 -0.31602189*(arg/3)**6 +0.04374224*(arg/3)**8 -0.00331563*(arg/3)**10
        N = (2/np.pi)*np.log(arg/2)*J +0.36746703 +0.60558498*(arg/3)**2 -0.74340225*(arg/3)**4 +0.25256673*(arg/3)**6 -0.04177345*(arg/3)**8 +0.00353354*(arg/3)**10
        S1 = 1.909859286*(arg/3)**2 -1.145914713*(arg/3)**4 +0.294656958*(arg/3)**6 -0.042070508*(arg/3)**8 +0.003785727*(arg/3)**10 -0.000207183*(arg/3)**12
        J1 = 0.49999999 -0.56249945*(arg/3)**2 +0.21093101*(arg/3)**4 -0.03952287*(arg/3)**6 +0.00439494*(arg/3)**8 -0.00028397*(arg/3)**10
        N1 = (2/np.pi)*(np.log(arg/2)*J1*arg - 1/arg) +0.07373571*(arg/3) +0.72276433*(arg/3)**3 -0.43885620*(arg/3)**5 +0.10418264*(arg/3)**7 -0.01340825*(arg/3)**9 +0.00094249*(arg/3)**11
        return -(V + np.pi*(S-N) - (np.pi*ep_o/ref)*rnorm*(S1-N1) )
    
    elif (arg >= 3):
        num = 0.99999906 + 4.77228920*(3/arg)**2 + 3.85542044*(3/arg)**4 + 0.32303607*(3/arg)**6
        den = arg*(1+ 4.88331068*(3/arg)**2 +4.28957333*(3/arg)**4 + 0.52120508*(3/arg)**6)
        S_N= (2/np.pi)*(num/den)
        num1 = 1.00000004 + 3.92205313*(3/arg)**2 + 2.64893033*(3/arg)**4 + 0.27450895*(3/arg)**6
        den1 = 1+ 3.81095112*(3/arg)**2 +2.26216956*(3/arg)**4 + 0.10885141*(3/arg)**6
        S_N1= (2/np.pi)*(num1/den1)
        return -(V + np.pi*(S_N) - (np.pi*ep_o/ref)*rnorm*(S_N1) )

#This computes the repulsive piece from the charge-phonon interaction potential   
@numba.jit(nopython=True)
def polPot_cross(r):
    rnorm = r
    V = (2*ep_o/ref)*rnorm
    arg = (ep_o/ref)*rnorm
    if (arg < 3):
        S = 1.909859164*(arg/3) -1.909855001*(arg/3)**3 +0.687514637*(arg/3)**5 -0.126164557*(arg/3)**7 +0.013828813*(arg/3)**9 -0.000876918*(arg/3)**11
        J = 0.9999999 -2.24999239*(arg/3)**2 +1.26553572*(arg/3)**4 -0.31602189*(arg/3)**6 +0.04374224*(arg/3)**8 -0.00331563*(arg/3)**10
        N = (2/np.pi)*np.log(arg/2)*J +0.36746703 +0.60558498*(arg/3)**2 -0.74340225*(arg/3)**4 +0.25256673*(arg/3)**6 -0.04177345*(arg/3)**8 +0.00353354*(arg/3)**10
        S1 = 1.909859286*(arg/3)**2 -1.145914713*(arg/3)**4 +0.294656958*(arg/3)**6 -0.042070508*(arg/3)**8 +0.003785727*(arg/3)**10 -0.000207183*(arg/3)**12
        J1 = 0.49999999 -0.56249945*(arg/3)**2 +0.21093101*(arg/3)**4 -0.03952287*(arg/3)**6 +0.00439494*(arg/3)**8 -0.00028397*(arg/3)**10
        N1 = (2/np.pi)*(np.log(arg/2)*J1*arg - 1/arg) +0.07373571*(arg/3) +0.72276433*(arg/3)**3 -0.43885620*(arg/3)**5 +0.10418264*(arg/3)**7 -0.01340825*(arg/3)**9 +0.00094249*(arg/3)**11
        return (V + np.pi*(S-N) - (np.pi*ep_o/ref)*rnorm*(S1-N1) )
    
    elif (arg >= 3):
        num = 0.99999906 + 4.77228920*(3/arg)**2 + 3.85542044*(3/arg)**4 + 0.32303607*(3/arg)**6
        den = arg*(1+ 4.88331068*(3/arg)**2 +4.28957333*(3/arg)**4 + 0.52120508*(3/arg)**6)
        S_N= (2/np.pi)*(num/den)
        num1 = 1.00000004 + 3.92205313*(3/arg)**2 + 2.64893033*(3/arg)**4 + 0.27450895*(3/arg)**6
        den1 = 1+ 3.81095112*(3/arg)**2 +2.26216956*(3/arg)**4 + 0.10885141*(3/arg)**6
        S_N1= (2/np.pi)*(num1/den1)
        return (V + np.pi*(S_N) - (np.pi*ep_o/ref)*rnorm*(S_N1) )    
    
#Helper function to construct the Virial estimator for the self-attractice piece    
@numba.jit(nopython=True)
def Vir_self(rm,rd):
    rnorm = rm
    rnast = rd
    V = (2*ep_o/ref)*(rnast/rnorm)
    arg = (ep_o/ref)*rnorm
    if (arg < 3):
        S = 1.909859164*(arg/3) -1.909855001*(arg/3)**3 +0.687514637*(arg/3)**5 -0.126164557*(arg/3)**7 +0.013828813*(arg/3)**9 -0.000876918*(arg/3)**11
        J = 0.9999999 -2.24999239*(arg/3)**2 +1.26553572*(arg/3)**4 -0.31602189*(arg/3)**6 +0.04374224*(arg/3)**8 -0.00331563*(arg/3)**10
        N = (2/np.pi)*np.log(arg/2)*J +0.36746703 +0.60558498*(arg/3)**2 -0.74340225*(arg/3)**4 +0.25256673*(arg/3)**6 -0.04177345*(arg/3)**8 +0.00353354*(arg/3)**10
        S1 = 1.909859286*(arg/3)**2 -1.145914713*(arg/3)**4 +0.294656958*(arg/3)**6 -0.042070508*(arg/3)**8 +0.003785727*(arg/3)**10 -0.000207183*(arg/3)**12
        J1 = 0.49999999 -0.56249945*(arg/3)**2 +0.21093101*(arg/3)**4 -0.03952287*(arg/3)**6 +0.00439494*(arg/3)**8 -0.00028397*(arg/3)**10
        N1 = (2/np.pi)*(np.log(arg/2)*J1*arg - 1/arg) +0.07373571*(arg/3) +0.72276433*(arg/3)**3 -0.43885620*(arg/3)**5 +0.10418264*(arg/3)**7 -0.01340825*(arg/3)**9 +0.00094249*(arg/3)**11
        return ( (2)*polPot_self(rm) -  (V/2 - (np.pi* ((ep_o)**2)/(2* (ref)**2))*rnast*(S-N) - (np.pi* ((ep_o))/(2* (ref)))*(rnast/rnorm)*(S1-N1 - 2/np.pi) )  )
    
    elif (arg >= 3):
        num = 0.99999906 + 4.77228920*(3/arg)**2 + 3.85542044*(3/arg)**4 + 0.32303607*(3/arg)**6
        den = arg*(1+ 4.88331068*(3/arg)**2 +4.28957333*(3/arg)**4 + 0.52120508*(3/arg)**6)
        S_N= (2/np.pi)*(num/den)
        num1 = 1.00000004 + 3.92205313*(3/arg)**2 + 2.64893033*(3/arg)**4 + 0.27450895*(3/arg)**6
        den1 = 1+ 3.81095112*(3/arg)**2 +2.26216956*(3/arg)**4 + 0.10885141*(3/arg)**6
        S_N1= (2/np.pi)*(num1/den1)
        return ( (2 )*polPot_self(rm)  - (V/2 - (np.pi* ((ep_o)**2)/(2* (ref)**2))*rnast*(S_N) - (np.pi* ((ep_o))/(2* (ref)))*(rnast/rnorm)*(S_N1- 2/np.pi) )   )
    

    
#Helper function to construct the Virial estimator for the repulsive piece   
@numba.jit(nopython=True)
def Vir_cross(r):
    rnorm = r
    
    V = (2*ep_o/ref)*rnorm
    arg = (ep_o/ref)*rnorm
    if (arg < 3):
        S = 1.909859164*(arg/3) -1.909855001*(arg/3)**3 +0.687514637*(arg/3)**5 -0.126164557*(arg/3)**7 +0.013828813*(arg/3)**9 -0.000876918*(arg/3)**11
        J = 0.9999999 -2.24999239*(arg/3)**2 +1.26553572*(arg/3)**4 -0.31602189*(arg/3)**6 +0.04374224*(arg/3)**8 -0.00331563*(arg/3)**10
        N = (2/np.pi)*np.log(arg/2)*J +0.36746703 +0.60558498*(arg/3)**2 -0.74340225*(arg/3)**4 +0.25256673*(arg/3)**6 -0.04177345*(arg/3)**8 +0.00353354*(arg/3)**10
        S1 = 1.909859286*(arg/3)**2 -1.145914713*(arg/3)**4 +0.294656958*(arg/3)**6 -0.042070508*(arg/3)**8 +0.003785727*(arg/3)**10 -0.000207183*(arg/3)**12
        J1 = 0.49999999 -0.56249945*(arg/3)**2 +0.21093101*(arg/3)**4 -0.03952287*(arg/3)**6 +0.00439494*(arg/3)**8 -0.00028397*(arg/3)**10
        N1 = (2/np.pi)*(np.log(arg/2)*J1*arg - 1/arg) +0.07373571*(arg/3) +0.72276433*(arg/3)**3 -0.43885620*(arg/3)**5 +0.10418264*(arg/3)**7 -0.01340825*(arg/3)**9 +0.00094249*(arg/3)**11
        return ( (4)*polPot_cross(r) +  (V - (np.pi* ((ep_o)**2)/( (ref)**2))*(rnorm**2)*(S-N) - (np.pi* ((ep_o))/( (ref)))*(rnorm)*(S1-N1- 2/np.pi) )  )
    
    elif (arg >= 3):
        num = 0.99999906 + 4.77228920*(3/arg)**2 + 3.85542044*(3/arg)**4 + 0.32303607*(3/arg)**6
        den = arg*(1+ 4.88331068*(3/arg)**2 +4.28957333*(3/arg)**4 + 0.52120508*(3/arg)**6)
        S_N= (2/np.pi)*(num/den)
        num1 = 1.00000004 + 3.92205313*(3/arg)**2 + 2.64893033*(3/arg)**4 + 0.27450895*(3/arg)**6
        den1 = 1+ 3.81095112*(3/arg)**2 +2.26216956*(3/arg)**4 + 0.10885141*(3/arg)**6
        S_N1= (2/np.pi)*(num1/den1)
        return ( (4)*polPot_cross(r) +  (V - (np.pi* ((ep_o)**2)/( (ref)**2))*(rnorm**2)*(S_N) - (np.pi* ((ep_o))/( (ref)))*(rnorm)*(S_N1- 2/np.pi) )  )

#Virial estimator for the self-attractive piece    
@numba.jit(nopython=True)
def Est_self(x):
    energ = 0.
    for i in range(N):
        for j in range(N):
            if (i != j):
                ix = x[i,0] - Lx*np.round(x[i,0] /Lx)
                iy = x[i,1] - Ly*np.round(x[i,1] /Ly)
                jx = x[j,0] - Lx*np.round(x[j,0] /Lx)
                jy = x[j,1] - Ly*np.round(x[j,1] /Ly)
            
                rmag = np.sqrt( (ix-jx)**2 + (iy-jy)**2  )
                rdot = (ix**2)+(iy**2) - ((ix*jx) +(iy*jy) )
                energ = energ + np.exp(-wlo*abs(i-j)/T/N)*(Vir_self(rmag,rdot) - (wlo/T/N)*abs(i-j)*polPot_self(rmag) )
    
    
    return ( (garb/T/N/N)*energ )

#Virial estimator for the repulsive piece
@numba.jit(nopython=True)
def Est_cross(x,h):
    energ = 0.
    for i in range(N):
        for j in range(N):
            if(i != j):
                ix = x[i,0] - Lx*np.round(x[i,0] /Lx)
                iy = x[i,1] - Ly*np.round(x[i,1] /Ly)
                jx = h[j,0] - Lx*np.round(h[j,0] /Lx)
                jy = h[j,1] - Ly*np.round(h[j,1] /Ly)
            
                rmag = np.sqrt( (ix-jx)**2 + (iy-jy)**2  )
                energ = energ + np.exp(-wlo*abs(i-j)/T/N)*(Vir_cross(rmag) - 2*(wlo/T/N)*abs(i-j)*polPot_cross(rmag) )
    
    
    return ( (garb/T/N/N)*energ )
        

    
  
#Monte Carlo Sweep
@numba.jit(nopython=True)
def Monte_Carlo(x,h,T,count):

 # Maximum size of trial displacement
 
 #Acceptance rate is kept to around 40%
   
    delta_xy = 6.2    
 

 # Preform 2N+2 random particle displacements
    for step in range(2*N + 2):
     
  # Draw random bead
        nn=np.random.randint(0,2*N+2)
        track = count
        
        #Propose to move single electron bead
        if (nn < N):
            # Draw a random displacement in x,y
            xn0=x[nn,0]*1.
            xn1=x[nn,0]+2*delta_xy*(np.random.rand()-0.5)
            
            
            yn0 = x[nn,1]*1. 
            yn1=x[nn,1]+2*delta_xy*(np.random.rand()-0.5)
            
            
            
            # Enforce MC acceptance
            

            energy=MetropolisB(x,xn0,xn1,yn0, yn1,nn,N,T,mass_e,h)
            if np.random.random()<=np.exp(-energy/T):
                x[nn,0]= (xn1 + 2*Lx)%Lx
                x[nn,1]=(yn1 + 2*Ly)%Ly
                track = track + 1
        
        #Propose to move single hole bead
        elif (nn > N and nn < 2*N):
            
            xn0=h[nn%N,0]*1.
            xn1=h[nn%N,0]+2*delta_xy*(np.random.rand()-0.5)
            
            
            yn0 = h[nn%N,1]*1. 
            yn1=h[nn%N,1]+2*delta_xy*(np.random.rand()-0.5)
            
            
            
            
            

            energy=MetropolisB(h,xn0,xn1,yn0, yn1,nn%N,N,T,mass_h,x)
            if np.random.random()<=np.exp(-energy/T):
                h[nn%N,0]= (xn1 + 2*Lx)%Lx
                h[nn%N,1]=(yn1 +2*Ly)%Ly
                track = track + 1
        
        #Propose to move whole electron ring polymer
        elif (nn == 2*N):
            x_New = x.copy()
            disp_x = 2*delta_xy*(np.random.rand()-0.5)
            disp_y = 2*delta_xy*(np.random.rand()-0.5)
            for i in range(N):
                x_New[i,0] = x_New[i,0] + disp_x 

                
                x_New[i,1] = x_New[i,1] + disp_y

            energy = MetropolisW(x,x_New,N,T,mass_e,h)
            if np.random.random()<=np.exp(-energy/T):
                x = x_New.copy()
                for i in range(N):
                    x[i,0] = (x[i,0]+2*Lx)%Lx
                    x[i,1] = (x[i,1]+2*Ly)%Ly
                    
                track = track + 1
                
        #Propose to move whole hole ring polymer            
        elif (nn == (2*N + 1) ):
            h_New = h.copy()
            disp_x = 2*delta_xy*(np.random.rand()-0.5)
            disp_y = 2*delta_xy*(np.random.rand()-0.5)
            for i in range(N):
                h_New[i,0] = h_New[i,0] + disp_x 

                
                h_New[i,1] = h_New[i,1] + disp_y


            energy = MetropolisW(h,h_New,N,T,mass_h,x)
            if np.random.random()<=np.exp(-energy/T):
                h = h_New.copy()
                for i in range(N):
                    h[i,0] = (h[i,0]+2*Lx)%Lx
                    h[i,1] = (h[i,1]+2*Ly)%Ly
                track = track + 1
            
    return x, h, track 

# Metropolis acceptance for single bead move
@numba.jit(nopython=True)
def MetropolisB(x,xn0,xn1,yn0, yn1, n,N,T,mass,h):
    fac = mass*N*T/2

    xr = (x[(n+1)%N,0]-xn0) - Lx*np.round((x[(n+1)%N,0]-xn0)/Lx)
    xl = (x[(n-1)%N,0]-xn0) - Lx*np.round((x[(n-1)%N,0]-xn0)/Lx)
    yr = (x[(n+1)%N,1]-yn0) - Ly*np.round((x[(n+1)%N,1]-yn0)/Ly)
    yl = (x[(n-1)%N,1]-yn0) - Ly*np.round((x[(n-1)%N,1]-yn0)/Ly)
    it =  (xr)**2 + (xl)**2 + (yl)**2 + (yr)**2 

  
    KEi=fac*T*(it)
    
    PEi = SN(xn0,yn0,h[n,0],h[n,1])/N 
    

    
    nx = xn0 - Lx*np.round(xn0 /Lx)  
    ny = yn0 - Ly*np.round(yn0 /Ly)
    
    PEi_self = 0.
    PEi_cross = 0.
    for j in range(N):
        if (j != n%N):
            ix = x[j,0] - Lx*np.round(x[j,0] /Lx) 
            iy = x[j,1] - Ly*np.round(x[j,1] /Ly)
            rtot = np.sqrt( (ix-nx)**2 + (iy-ny)**2  )
            PEi_self = PEi_self + np.exp(-wlo*abs(n-j)/T/N)*polPot_self(rtot)
            
            jx = h[j,0] - Lx*np.round(h[j,0] /Lx) 
            jy = h[j,1] - Ly*np.round(h[j,1] /Ly)
            rtot = np.sqrt( (jx-nx)**2 + (jy-ny)**2  )
            PEi_cross = PEi_cross + np.exp(-wlo*abs(n-j)/T/N)*polPot_cross(rtot)

   
    xr = (x[(n+1)%N,0]-xn1) - Lx*np.round((x[(n+1)%N,0]-xn1)/Lx)
    xl = (x[(n-1)%N,0]-xn1) - Lx*np.round((x[(n-1)%N,0]-xn1)/Lx)
    yr = (x[(n+1)%N,1]-yn1) - Ly*np.round((x[(n+1)%N,1]-yn1)/Ly)
    yl = (x[(n-1)%N,1]-yn1) - Ly*np.round((x[(n-1)%N,1]-yn1)/Ly)
    fi =  (xr)**2 + (xl)**2 + (yl)**2 + (yr)**2 
    

 
    KEf=fac*T*(fi)
    PEf = SN(xn1,yn1,h[n,0],h[n,1])/N 
    

    
  
    nx = xn1 - Lx*np.round(xn1 /Lx)  
    ny = yn1 - Ly*np.round(yn1 /Ly)
    
     
    PEf_self = 0.
    PEf_cross = 0.
    for j in range(N):
        if (j != n%N):
            ix = x[j,0] - Lx*np.round(x[j,0] /Lx) 
            iy = x[j,1] - Ly*np.round(x[j,1] /Ly)
            rtot = np.sqrt( (ix-nx)**2 + (iy-ny)**2  )
            PEf_self = PEf_self + np.exp(-wlo*abs(n-j)/T/N)*polPot_self(rtot)
            
            
            jx = h[j,0] - Lx*np.round(h[j,0] /Lx) 
            jy = h[j,1] - Ly*np.round(h[j,1] /Ly)
            rtot = np.sqrt( (jx-nx)**2 + (jy-ny)**2  )
            PEf_cross = PEf_cross + np.exp(-wlo*abs(n-j)/T/N)*polPot_cross(rtot)
        
    
    sumE=KEf+PEf -KEi - PEi + (garb/T/N/N)*2*(PEf_self + PEf_cross - PEi_self - PEi_cross)
    

    return (sumE)

# Metropolis acceptance for ring polymer move
@numba.jit(nopython=True)
def MetropolisW(x,x_New,N,T,mass,h):
    sumPE = 0
    
    for i in range(N):
        PEi=SN(x[i,0],x[i,1],h[i,0],h[i,1])/N
        PEf=SN(x_New[i,0],x_New[i,1],h[i,0],h[i,1])/N
        sumPE = sumPE + (PEf-PEi)
        
    

            
            
    sumPE_cross = 0        
    for i in range(N):
        for j in range(N):
            if (i!=j):
                ix = x[i,0] - Lx*np.round(x[i,0] /Lx)
                iy = x[i,1] - Ly*np.round(x[i,1] /Ly)
                jx = h[j,0] - Lx*np.round(h[j,0] /Lx)
                jy = h[j,1] - Ly*np.round(h[j,1] /Ly)
                rmag = np.sqrt( (ix-jx)**2 + (iy-jy)**2  )
                PEi = np.exp(-wlo*abs(i-j)/T/N)*polPot_cross(rmag)
            
                ix = x_New[i,0] - Lx*np.round(x_New[i,0] /Lx)
                iy = x_New[i,1] - Ly*np.round(x_New[i,1] /Ly)
                rmag = np.sqrt( (ix-jx)**2 + (iy-jy)**2  )
                PEf = np.exp(-wlo*abs(i-j)/T/N)*polPot_cross(rmag)
            
                sumPE_cross = sumPE_cross + (PEf-PEi)      
        
        
        
        
        
    return (sumPE + (garb/T/N/N)*( 2*sumPE_cross )  )


# Parameters #
##
##

# Total number of beads
Nx = 20
Ny = 10
N = Nx*Ny

#Size of Periodic Box
Lx = 500
Ly = 500

#Effective electron and hole masses
mass_e=0.20
mass_h = 0.20

#Dielectric constants of inorganic and organic layers
ep= 6.1
ep_o = 3.3



#Variables needed for the charge-phonon interaction potential and Virial estimator
ro = ep*Lz/(2*ep_o)
ref = Lz*ep_o*( (1/3)*(ep_o/ep) + (0.5*ep/ep_o)*(1-  (ep_o/ep)**2  )  )

ep_st=13.0
ep_in = ep
wlo = 0.0005144902693169139 #wlo in atomic units
garb = Lz*wlo*(ep_st-ep_in)/(16* (ref**2) )

# Number of sweeps
steps = 500000
# Averaging frequency
Nsamp = 500


r_s  = np.linspace(0,500,10000,endpoint=False)
potSums = np.array(sums)
derivSums = np.array(Est_ar)



print ("MC Run with N=", N, "Number of MC steps=", steps)
print("Lx and Ly = ", Lx)


# Container for positions
x = np.zeros([N,2])
h = np.zeros([N,2])

# Initialize ring polymers
x,h =init(x,h)

#Temperature in atomic units
T = 0.0003167*4.0
print("Temperature: ", T*100/0.0003167)

#Store r's for P(r)
r_m1 = []

#Store energies
Ener1 = []

count = 0 
    # Perform Nsweeps 
for step in range(steps):

    # Perform 2N+2 random displacements
    x, h, count =Monte_Carlo(x,h,T,count)
        

    # Every Nsamp steps, compute expectations
    if(step%Nsamp==0):
        en = 0
        r_sim = np.zeros(N)

        for i in range(N):
            en = en + Est(x[i,0],x[i,1],h[i,0],h[i,1])
            xu = (x[i,0]-h[i,0]) - Lx*np.round((x[i,0]-h[i,0]) /Lx)
            yu = (x[i,1]-h[i,1]) - Ly*np.round((x[i,1]-h[i,1]) /Ly)
            r_sim[i] = np.sqrt( (xu)**2 + (yu)**2   )
                
                
        #To get bare exciton binding energies,you only need en/N below
        
        #You can also comment out all of the charge-phonon interaction pieces
        #from the MC acceptance functions to speed up the bare exciton calculations
        Ener1.append(en/N + Est_self(x) + Est_self(h) + Est_cross(x,h))
        
        r_m1.append(r_sim)

#        print("Step: ", step, flush=True)           
        if(step >= 2000000):
            np.savetxt(Ener1)
            
            
print("Acceptance Ratio: ", count/steps)
            
#Print out mean energy. A plot of the mean energy as a function of MC sweeps (running average plot)
#can help determine  when equilibration occurs
print(r"<E> (meV): ",mean(Ener1[1000:])*1000*27.211399 )





