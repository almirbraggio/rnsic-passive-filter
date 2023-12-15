import scipy.integrate as integrate
from numpy import sqrt, sin, cos, pi, inf, arccos
from engineering_notation import EngNumber as eng

def f_eu(umain,omegat):
    return (umain*sin(omegat))

def f_ev(umain,omegat):
    return (umain*sin(omegat-((2/3)*pi)))

def f_ew(umain,omegat):
    return (umain*sin(omegat+((2/3)*pi)))

def f_imain(umain, pout):
    return (pout*(2/3)/umain)

def f_l1(umain=0, imain=0, freq=0, uout=0, omegat1=0):
    num_int1 = integrate.quad(lambda omegat: (f_ew(umain,omegat)-f_ev(umain,omegat)-uout), \
    (omegat1-pi/3), (pi/3))
    den_p1 = imain*2*pi*freq
    den_sin1 = sin(omegat1 + pi/3)
    den_sin2 = sin(omegat1)
    den_total = den_p1*(0-den_sin1-den_sin2+(sqrt(3)/2))
    return (num_int1[0]/den_total)

def f_c1(imain=0, freq=0, uout=0, omegat1=0, l1=0):
    num_int1 = integrate.quad(lambda omegat: (imain*(1-cos(omegat))), (omegat1-pi/3), (pi/3))
    den_p1 = l1*2*pi*freq*imain*(sin(omegat1)-sqrt(3)/2)
    den_p2 = (uout/3)*(0-omegat1+(2/3)*pi)
    den_int1 = integrate.quad(lambda omegat: (f_ev(umain, omegat)), (omegat1-pi/3), (pi/3))
    den_total = 6*2*pi*freq*(den_p1-den_p2-den_int1[0])
    return (num_int1[0]/den_total)

def f_l2(umain=0, imain=0, freq=0, uout=0, omegat1=0):
    num_int1 = integrate.quad(lambda omegat: (f_ev(umain,omegat)-f_ew(umain,omegat)+uout), \
    (pi/3), (omegat1))
    
    mult_num_int2 = (umain*((1/2)-cos(omegat1-((2/3)*pi))))+((uout/3)*(omegat1-(pi/3)))
    num_int2 = integrate.quad(lambda omegat: ((mult_num_int2)*((3*(1+cos(omegat+((2/3)*pi))))/(sin(omegat1)+sin(omegat+((2/3)*pi))-(sqrt(3)/2)))), \
    (pi/3), (omegat1), epsabs=inf)
    num_total = num_int1[0]-num_int2[0]
    
    den_p1 = (imain*2*pi*freq)
    den_sin1 = sin(omegat1-(2/3)*pi)-sin(omegat1+(2/3)*pi)+(sqrt(3)/2)
    den_p1_total = den_p1*den_sin1

    mult_den_int1 = (sin(omegat1-((2/3)*pi))+(sqrt(3)/2))
    den_int1 = integrate.quad(lambda omegat: ((mult_den_int1)*((3*imain*2*pi*freq*(1+cos(omegat+(2/3)*pi)))/(sin(omegat1)+sin(omegat+(2/3)*pi)-(sqrt(3)/2)))), \
    (pi/3), (omegat1), epsabs=inf)
    den_total = den_p1_total-den_int1[0]
    return (num_total/den_total)

def f_c2(imain=0, freq=0, uout=0, omegat1=0, l1=0):
    num_int1 = integrate.quad(lambda omegat: (imain*(1+cos(omegat+((2/3)*pi)))), \
    (pi/3), (omegat1))

    den_p1 = (3*2*pi*freq)
    den_p2 = (imain*l1*2*pi*freq*sin(omegat1-((2/3)*pi)))
    den_int1_mult = (2/3)*uout
    den_int1 = integrate.quad(lambda omegat: (den_int1_mult+f_ev(umain,omegat)), \
    (pi/3), (omegat1))
    den_total = den_p1*(den_int1[0]+den_p2)
    return (num_int1[0]/den_total)

def f_avg(l1c1=0, l2c2=0, omegat1=0):
    num_p1 = l1c1*(((2/3)*pi)-omegat1)
    num_p2 = l2c2*(omegat1-pi/3)
    return ((num_p1+num_p2)/(pi/3))

umain = 57.15       # amplitude da tensão de fase na entrada do retificador
uout = 180.0        # amplitude da tensão de saída do retificador
pout = 5320.0       # potencia de saída em watt
freq = 60.0         # freq em hertz

imain = f_imain(umain, pout)
#omegat1 = 0.518*pi # omega*t1 calculado
omegat1 = arccos(pi*umain/uout - 1)
print("omegat1 = " + str(omegat1/pi) + "*pi")

l1 = f_l1(umain, imain, freq, uout, omegat1)
c1 = f_c1(imain, freq, uout, omegat1, l1)

l2 = f_l2(umain, imain, freq, uout, omegat1)
c2 = f_c2(imain, freq, uout, omegat1, l1)

print (\
    "L1 = " + str(eng(l1)) + "H \t" + \
    "L2 = " + str(eng(l2)) + "H \t" + \
    "AVG = "+ str(eng(f_avg(l1, l2, omegat1))) + "H")
print (\
    "C1 = " + str(eng(c1)) + "F \t" + \
    "C2 = " + str(eng(c2)) + "F \t" + \
    "AVG = "+ str(eng(f_avg(c1, c2, omegat1))) + "F")
