import numpy
import random
import matplotlib.pyplot as pyplot

#CONSTANTS
alpha= 0.0072973525664 # Fine structure constant 
Gf= 0.000011663787 # Fermi coupling constant in GeV^-2
Mz= 91.1876 #Z boson mass in GeV/c^2
Yz= 2.4952 #Total decay width of Z boson in GeV/c^2
sinsqthetaW= 0.2223 #sin of the weak mixing angle squared 
pi=numpy.pi

def dsigma(s,costheta,Qf):
    
    dsig= ((pi*(alpha**2)*(Qf**2))/(2*s))*(1+(costheta**2))
    return dsig
    
    
    
def dsigmaZ(s,costheta,Qf,Tfcubed):
    '''Return d(sigma)/d(costheta) for a e+e- annihilation into two massless fermions, f at given CoM energy, s.
    Accounts for Z channel. Formula from pg 55 of QCD and Collider Physics, R. K. Ellis'''
    k= ((2**0.5)*Gf*(Mz**2))/(16*pi*alpha)
    x1= k*(s*(s-(Mz**2)))/((s-(Mz**2))**2 +(Yz*Mz)**2)
    x2= (k**2)*(s**2)/((s-(Mz**2))**2 +(Yz*Mz)**2)
    Ve= -0.5- 2*Qf*sinsqthetaW
    Vf= -Tfcubed- 2*Qf*sinsqthetaW
    Ae= -0.5
    Af= Tfcubed
    
    dsig= ((pi*(alpha**2))/(2*s))*((1+(costheta**2))*((Qf**2)-(2*Qf*Ve*Af*x1)+ 
    ((Ae**2)+(Ve**2))*((Af**2)+(Vf**2))*x2)) + costheta*((-4)*Qf*Ae*Af*x1 + 8*Ae*Ve*Af*Vf*x2)
    
    return dsig

def MCintegration(N,s,Z):
    '''Returns MC integral of dsigmaZ with CoM energy s and N MC integration points used'''
    integral=0
     
    if Z == 1:
        for i in range(N):
            r=2*random.random()-1 #generating random value of costheta [-1,1] 
            integral= integral+ 2*dsigmaZ(s,r,-1,-0.5)  #Evaluating dsigmaZ at random value of costheta
        ans= integral/N
    else:
        for i in range(N):
            r=2*random.random()-1 #generating random value of costheta [-1,1] 
            integral= integral+ 2*dsigma(s,r,-1)  #Evaluating dsigmaZ at random value of costheta
        ans= integral/N    
    return ans    

crosssection= numpy.zeros(100000)
crosssectionZ= numpy.zeros(100000)
svalues= numpy.arange(1,10001,0.1)
crosssectionZ= MCintegration(10000,svalues,1) #Evaluating MC integral for a range of CoM energy values
crosssection= MCintegration(10000,svalues,0) 

pyplot.figure()
g1= pyplot.subplot(211)
#pyplot.plot(svalues**0.5, crosssection, color= 'blue')
g1.plot(svalues**0.5, crosssectionZ, color= 'red')

g2= pyplot.subplot(212)
g2.plot(svalues[0:5000]**0.5, crosssectionZ[0:5000], color= 'red')

pyplot.xlabel('CoM energy $GeV/c^2$')
pyplot.ylabel(r'$\sigma  (e^-e^+\to\mu^-\mu^+)$', size=15)
pyplot.show()