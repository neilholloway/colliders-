import numpy
import random
import matplotlib.pyplot as pyplot

pi= numpy.pi

def momentagen(N,w):
    '''Generates an array of N 4 momenta, for N massless outgoing particles according to MC algorithm given in
    A NEW MONTE CARLO TREATMENT OF MULTIPARTICLE PHASE SPACE AT HIGH ENERGIES, R. KLEISS and W.J. STIRLING'''
    
    momenta = numpy.zeros((N,4))
    for i in range(N):
        r1= random.random()
        r2= random.random()
        r3= random.random()
        r4= random.random()
        
        #Generating the 4 momenta q as described in paper
        c= 2*r1 -1
        psi= 2*pi*r2
        E= -numpy.log(r3*r4)
        momenta[i,0]= E
        momenta[i,1]= E*numpy.cos(psi)*((1-c**2)**0.5)
        momenta[i,2]= E*numpy.sin(psi)*((1-c**2)**0.5)
        momenta[i,3]= E*c
        
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
    
    transmom= numpy.zeros((N,4))
    for j in range(N):
        transmom[j,0]= x*(y*momenta[j,0]+numpy.dot(momenta[j,1:4],b))   
        transmom[j,1:4]= x*(momenta[j,1:4]+momenta[j,0]*b+a*numpy.dot(b,momenta[j,1:4])*b)
        #transmom[j,2]= x*(momenta[j,2]+a*(b[1]*momenta[j,2])*b[1])
        #transmom[j,3]= x*(momenta[j,3]+a*(b[2]*momenta[j,3])*b[2])
    
    P4mom= numpy.sum(transmom,axis=0)
    p1mom= numpy.array((0,0,w/2.))
    p3mom= transmom[0,1:4]
    
    costheta= numpy.dot(p1mom,p3mom)/(numpy.sqrt(numpy.dot(p1mom,p1mom)*numpy.dot(p3mom,p3mom)))
    #print P4mom
    return costheta
    
print momentagen(2,1000)
    