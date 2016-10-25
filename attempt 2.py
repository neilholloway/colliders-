import numpy
import random
import matplotlib.pyplot as pyplot
import time

start_time = time.time()

#CONSTANTS
alpha= 0.0072973525664 # Fine structure constant 
pi=numpy.pi
Mw=80.4 
Mz= 91.1876 #Z boson mass in GeV/c^2
Yz= 2.4952 #Total decay width of Z boson in GeV/c^2
sinsqthetaW= 0.2223 #sin of the weak mixing angle squared 
#Gf= 0.000011663787 # Fermi coupling constant in GeV^-2
Gf= 4*pi*alpha*numpy.sqrt(2.)/(8*sinsqthetaW*Mw**2)

#RAMBO CODE
def momentagen(outgoing,w):
    '''Generates an array of N 4 momenta, for N massless outgoing particles according to MC algorithm given in
    A NEW MONTE CARLO TREATMENT OF MULTIPARTICLE PHASE SPACE AT HIGH ENERGIES, R. KLEISS and W.J. STIRLING'''
    
    momenta = numpy.zeros((outgoing,4))
    for n in range(outgoing):
        r1= random.random()
        r2= random.random()
        r3= random.random()
        r4= random.random()
        
        #Generating the 4 momenta q as described in paper
        c= 2*r1 -1
        psi= 2*pi*r2
        E= -numpy.log(r3*r4)
        momenta[n,0]= E
        momenta[n,1]= E*numpy.cos(psi)*((1-c**2)**0.5)
        momenta[n,2]= E*numpy.sin(psi)*((1-c**2)**0.5)
        momenta[n,3]= E*c
        
    Qmom= numpy.zeros((1,3))
    Q4mom= numpy.sum(momenta,axis=0) #Sum of generated 4 momenta q, labelled Q in paper
    Q0= Q4mom[0] 
    Qmom= numpy.array([Q4mom[1],Q4mom[2],Q4mom[3]])
    
    #Transforming generated massless 4 momenta set q into set p with CoM energy s
    M= numpy.sqrt(Q4mom[0]**2-Q4mom[1]**2-Q4mom[2]**2-Q4mom[3]**2)
    b= -Qmom/M
    x= w/M
    y= Q0/M
    a=1/(1+y)
    
    transmom= numpy.zeros((outgoing,4))
    for j in range(outgoing):
        transmom[j,0]= x*(y*momenta[j,0]+numpy.dot(momenta[j,1:4],b))   
        transmom[j,1:4]= x*(momenta[j,1:4]+momenta[j,0]*b+a*numpy.dot(b,momenta[j,1:4])*b)
    
    P4mom= numpy.sum(transmom,axis=0)
    p1mom= numpy.array((0,0,w/2.))
    p3mom= transmom[0,1:4]
    
    costheta= numpy.dot(p1mom,p3mom)/(numpy.sqrt(numpy.dot(p1mom,p1mom)*numpy.dot(p3mom,p3mom)))

    return costheta

Te3 = -1./2.
Tf3 = -1./2.
Qe = 1.
Qf = 1.
      
def dsigmaZ(s,cosinetheta,Z_on=True):
    '''Return d(sigma)/d(costheta) for a e+e- annihilation into two massless fermions, f at given CoM energy, s.
    Accounts for Z channel. Formula from pg 55 of QCD and Collider Physics, R. K. Ellis'''
    k= (numpy.sqrt(2.)*Gf*Mz**2)/(16.*pi*alpha)
    x1= k*s*(s-Mz**2)/((s-Mz**2)**2 +(Yz*Mz)**2)
    x2= (k*s)**2/((s-Mz**2)**2 +(Yz*Mz)**2)
    Ve= Te3 - 2.*Qe*sinsqthetaW
    Vf= Tf3 - 2.*Qf*sinsqthetaW
    Ae= Te3
    Af= Tf3
    if not Z_on:
        Ve = Vf = Ae = Af = 0.
    
    term1 = (1+cosinetheta**2)*(Qf**2-2*Qf*Ve*Af*x1+(Ae**2+Ve**2)*(Af**2+Vf**2)*x2)
    term2 =  cosinetheta*(-4*Qf*Ae*Af*x1 + 8*Ae*Ve*Af*Vf*x2)
    dsig= ((pi*alpha**2)/(2*s))*(term1 + term2)
        
    return dsig

def MCintegration(N,s):
    '''Returns MC integral of dsigmaZ with CoM energy s and N MC integration points used'''
    integral=0
    for i in range(N):
        r= momentagen(2,numpy.sqrt(s))
        integral= integral+ 2*dsigmaZ(s,r)  #Evaluating dsigmaZ at random value of costheta
    ans= integral/N
        
    return ans    

points= 10000
rootsmin, rootsmax= 400, 10400
crosssection= numpy.zeros(points)
crosssectionZ= numpy.zeros(points)
svalues= numpy.arange(rootsmin,rootsmax,(rootsmax-rootsmin)/points)

for i in range(points):
    crosssectionZ[i]= MCintegration(100,svalues[i])*(3.894*10**8) #Evaluating MC integral for a range of CoM energy values


pyplot.figure()
g1= pyplot.subplot(111)
g1.semilogy(svalues**0.5, crosssectionZ, color= 'red')
pyplot.xlabel('CoM energy, $GeV/c^2$')
pyplot.ylabel(r'$\sigma  (e^-e^+\to\mu^-\mu^+),$ $pb$', size=15)
pyplot.show()

print "Run Time =", float((time.time() - start_time)/60), "minutes"