import scipy.integrate as integrate
from numpy import sqrt, sin, cos, inf, arccos, real, imag, abs
from engineering_notation import EngNumber as eng

pi = 3.1415#9265

def f_eu(umain,omegat):
    return (umain*sin(omegat))

def f_ev(umain,omegat):
    return (umain*sin(omegat-(2.0*pi/3.0)))

def f_ew(umain,omegat):
    return (umain*sin(omegat+(2.0*pi/3.0)))

def f_imain(umain, pout):
    return (pout*(2.0/3.0)/umain)

def f_l1(umain=0.0, imain=0.0, freq=0.0, uout=0.0, omegat1=0.0):
    num_int1 = integrate.quad(lambda omegat: (f_ew(umain,omegat)-f_ev(umain,omegat)-uout), \
    (omegat1-pi/3.0), (pi/3.0))
    den_p1 = imain*2.0*pi*freq
    den_sin1 = sin(omegat1 + pi/3.0)
    den_sin2 = sin(omegat1)
    den_total = den_p1*(0.0-den_sin1-den_sin2+(sqrt(3.0)/2.0))
    return (num_int1[0]/den_total)

def f_c1(imain=0.0, freq=0.0, uout=0.0, omegat1=0.0, l1=0.0):
    num_int1 = integrate.quad(lambda omegat: (imain*(1.0-cos(omegat))), (omegat1-pi/3.0), (pi/3.0))
    den_p1 = l1*2.0*pi*freq*imain*(sin(omegat1)-sqrt(3.0)/2.0)
    den_p2 = (uout/3.0)*(0-omegat1+(2.0/3.0)*pi)
    den_int1 = integrate.quad(lambda omegat: (f_ev(umain, omegat)), (omegat1-pi/3.0), (pi/3.0))
    den_total = 6.0*2.0*pi*freq*(den_p1-den_p2-den_int1[0])
    return (num_int1[0]/den_total)

def f_l2(umain=0.0, imain=0.0, freq=0.0, uout=0.0, omegat1=0.0):
    num_int1 = integrate.quad(lambda omegat: (f_ev(umain,omegat)-f_ew(umain,omegat)+uout), \
    (pi/3.0), (omegat1))
    
    mult_num_int2 = (umain*((1.0/2.0)-cos(omegat1-((2.0/3.0)*pi))))+((uout/3.0)*(omegat1-(pi/3.0)))
    num_int2 = integrate.quad(lambda omegat: ( (mult_num_int2) * ( (3.0*(1.0+cos(omegat+(2.0*pi/3.0)))) / (sin(omegat1)+sin(omegat+(2.0*pi/3.0))-(sqrt(3.0)/2.0)) ) ), \
    (pi/3.0), (omegat1), epsabs=inf)
    num_total = num_int1[0]-num_int2[0]
    
    den_p1 = (imain*2.0*pi*freq)
    den_sin1 = sin(omegat1-(2.0*pi/3.0)) - sin(omegat1+(2.0*pi/3.0)) + (sqrt(3.0)/2.0)
    den_p1_total = den_p1*den_sin1

    mult_den_int1 = sin(omegat1-((2.0/3.0)*pi)) + (sqrt(3.0)/2.0)
    den_int1 = integrate.quad(lambda omegat: ( (mult_den_int1) * ( (3.0*imain*2.0*pi*freq*(1.0+cos(omegat+(2.0*pi/3.0)))) / (sin(omegat1)+sin(omegat+(2.0*pi/3.0))-(sqrt(3.0)/2.0)) ) ), \
    (pi/3.0), (omegat1), epsabs=inf)
    den_total = den_p1_total-den_int1[0]
       
    return abs(num_total/den_total)

def f_c2(imain=0.0, freq=0.0, uout=0.0, omegat1=0.0, l1=0.0):
    num_int1 = integrate.quad(lambda omegat: (imain*(1.0+cos(omegat+(2.0*pi/3.0)))), \
    (pi/3.0), (omegat1))

    den_p1 = (3.0*2.0*pi*freq)
    den_p2 = (imain*l1*2.0*pi*freq*sin(omegat1-(2.0*pi/3.0)))
    den_int1_mult = (2.0*uout/3.0)
    den_int1 = integrate.quad(lambda omegat: (den_int1_mult+f_ev(umain,omegat)), \
    (pi/3.0), (omegat1))
    den_total = den_p1*(den_int1[0]+den_p2)
    return abs(num_int1[0]/den_total)

def f_opt_l_c(x1=0.0, x2=0.0, omegat1=0.0):
    ret = x1*((2.0*pi/3.0)-(omegat1)) + x2*((omegat1)-(pi/3.0))
    return abs(ret/(pi/3.0))

umain = 180.0#*sqrt(2)       # amplitude da tensão de fase na entrada do retificador
uout = 350.0        # amplitude da tensão de saída do retificador
pout = 10000.0       # potencia de saída em watt
freq = 60.0         # freq em hertz

print (\
    "Upk = " + str(eng(umain)) + "V (" + str(eng(freq)) + " Hz)\t" + \
    "Uout = " + str(eng(uout)) + "V \t" + \
    "Pout=Pin = "+ str(eng(pout)) + "W \t" )

imain = f_imain(umain, pout)

omegat1_arr = [\
    # RNSIC-1
    arccos((pi*umain/uout) - 1.0), \
    # RNSIC-2
    arccos((pi*umain)/(uout*sqrt(3.0))) + (pi/6.0), \
    # RNSIC-3
    arccos((pi*umain)/(uout*sqrt(3.0)) - 1.0) - (pi/6.0)]

for omegat1 in omegat1_arr:
    try:
        print("\nomegat1 = " + str(eng(omegat1/pi)) + "*pi ")

        l1 = f_l1(umain, imain, freq, uout, omegat1)
        c1 = f_c1(imain, freq, uout, omegat1, l1)

        l2 = f_l2(umain, imain, freq, uout, omegat1)
        c2 = f_c2(imain, freq, uout, omegat1, l1)

        print (\
            "L1 = " + str(eng(l1)) + "H \t" + \
            "L2 = " + str(eng(l2)) + "H \t" + \
            "AVG => L(opt) = " + str(eng(f_opt_l_c(l1, l2, omegat1))) + "H")
        print (\
            "C1 = " + str(eng(c1)) + "F \t" + \
            "C2 = " + str(eng(c2)) + "F \t" + \
            "AVG => C(opt) = " + str(eng(f_opt_l_c(c1, c2, omegat1))) + "F")
    except:
        print("\n>>> invalid condition")
