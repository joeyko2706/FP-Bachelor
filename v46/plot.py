import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
from uncertainties import ufloat

def linregress(x, y):
    N = len(y)
    Delta = N * np.sum(x**2) - (np.sum(x))**2
    A = (N * np.sum(x * y) - np.sum(x) * np.sum(y)) / Delta
    B = (np.sum(x**2) * np.sum(y) - np.sum(x) * np.sum(x * y)) / Delta
    sigma_y = np.sqrt(np.sum((y - A * x - B)**2) / (N - 2))
    A_error = sigma_y * np.sqrt(N / Delta)
    B_error = sigma_y * np.sqrt(np.sum(x**2) / Delta)
    # print(A, A_error, B, B_error)

    return A*x + B

##BESTIMMUNG DER FLUSSDICHTE
a, B = np.genfromtxt("data/Flussdichte.txt", unpack = True)
xfluss = [-10, 10]
yfluss = [411, 411]

plt.figure()
plt.plot(a, B, 'x', label='Messwerte der Flussdichte')
plt.plot(xfluss, yfluss, label='Maximaler Wert')
plt.ylabel('B $[mT]$')
plt.xlabel('a $[mm]$')
plt.grid()
plt.legend(loc='best')
# plt.show()
plt.savefig('build/plot-Fluss.pdf')


g1, m1, lam = np.genfromtxt("data/GaAsRein.txt", unpack = True)
g2, m2 = np.genfromtxt("data/GaAsDot1.txt", unpack = True)
g3, m3 = np.genfromtxt("data/GaAsDot2.txt", unpack = True)

lamb = lam[:9] 
n = 1
while n < 10:
    lamb[n-1] = lam[2*n-1]
    n += 1

w1 = g1 + m1 / 60 # Winkel Dezimal in Grad
w2 = g2 + m2 / 60
w3 = g3 + m3 / 60

rw1 = w1  * 2 * np.pi / 360
rw2 = w2  * 2 * np.pi / 360
rw3 = w3  * 2 * np.pi / 360

wd1 = w1[:9] *2 #Arrays erstellen
wd2 = w2[:9] *2
wd3 = w3[:9] *2

n = 1
while n < 10:
    wd1[n - 1] = 0.5 * np.abs( rw1[2*n-2] - rw1[2*n-1])
    wd2[n - 1] = 0.5 * np.abs( rw2[2*n-2] - rw2[2*n-1])
    wd3[n - 1] = 0.5 * np.abs( rw3[2*n-2] - rw3[2*n-1])
    n += 1

d1 = 5.1*10**(-3)
d2 = 1.36*10**(-3)
d3 = 1.296*10**(-3)

nwd1 = wd1 / d1 #normieren
nwd2 = wd2 / d2
nwd3 = wd3 / d3
# print(nwd1, nwd2, nwd3)

#Tabellen:
# n = 1
# while n < 10:
#     print( lamb[n-1], ' & ', np.round(w1[2*n-2], 2), ' & ', np.round(w1[2*n-1], 2), ' & ', np.round(nwd1[n-1], 2), ' \')
#     n += 1
# n = 1
# while n < 10:
#     print( lamb[n-1], ' & ', np.round(w2[2*n-2], 2), ' & ', np.round(w2[2*n-1], 2), ' & ', np.round(nwd2[n-1], 2), ' \')
#     n += 1
# n = 1
# while n < 10:
#     print( lamb[n-1], ' & ', np.round(w3[2*n-2], 2), ' & ', np.round(w3[2*n-1], 2), ' & ', np.round(nwd3[n-1], 2), ' \')
#     n += 1


plt.figure()
plt.plot(lamb**2, nwd1, 'x', label='Erste Probe')
plt.plot(lamb**2, nwd2, 'x', label='Zweite Probe')
plt.plot(lamb**2, nwd3, 'x', label='Dritte Probe')
plt.xlabel(r'$\lambda^2 [\mu m^{2}]$')
plt.ylabel(r'$\theta [\frac{rad}{\si{\metre}}]$')      
plt.legend(loc='best')
plt.grid()
plt.savefig('build/plot-Werte.pdf')
# plt.show()
plt.close()

diff1 = nwd1 * 0
diff2 = nwd1 * 0

n = 0
while n < 9:
    diff1[n] = np.abs(nwd1[n] - nwd2[n])
    diff2[n] = np.abs(nwd1[n] - nwd3[n])
    n += 1


indexes = [6, 7]
rdiff1 = np.delete(diff1, indexes)
rdiff2 = np.delete(diff2, indexes)
rlamb = np.delete(lamb, indexes)



plt.figure()
plt.plot(lamb**2, diff1, 'x', label='Zweite Probe')
plt.plot(rlamb**2, linregress(rlamb**2, rdiff1), label='1.lineare Regression')
plt.plot(lamb**2, diff2, 'x', label='Dritte Probe')
plt.plot(rlamb**2, linregress(rlamb**2, rdiff2), label='2.lineare Regression')
plt.xlabel(r'$\lambda^2 [\mu m^{2}]$')
plt.ylabel(r'$\theta [\frac{rad}{\si{\metre}}]$')
plt.legend(loc='best')
plt.grid()
plt.savefig('build/plot-linReg.pdf')
# plt.show()

def get_m(a, N):
    return ((const.e**(3)*N*B)/(8*const.pi**(2)*const.c**(3)*const.epsilon_0*n*a))**(0.5)

B = 0.411
N1 = (1.2*10**(12))/(10**(-2))**3
N2 = (2.8*10**(12))/(10**(-2))**3
N3 = (1.2*10**(12))
N4 = (2.8*10**(12))
n = 3.3543
# a1 = ufloat(7160420549567.41, 5948184243297.86)
# a2 = ufloat(3508831908764.56, 6558733253011.16)

a1 = ufloat(15.6, 4.56)
a2 = ufloat(13.6, 3.83)

m1 = get_m(a1, N3)
m2 = get_m(a2, N4)

# print(m1, m2)

m_1 = m1/const.electron_mass
m_2 = m2/const.electron_mass

# print(m_1, m_2)







##BESTIMMUNG DER ROTATION: g Grad, s Sekunde, w Winkel(Dezimal), mw Winkelmittelwerte, bw Winkel im Bogenmaß, wd Winkeldifferenzen, lw Differenzen mit undotierter Probe
g1, s1 = np.genfromtxt("data/GaAshochrein.txt", unpack = True)
w1 = g1 + s1/60
mw1 = [0] * 18
g2, s2 = np.genfromtxt("data/GaAsDotiert1.txt", unpack = True)
w2 = g2 + s2/60
mw2 = [0] * 18
g3, s3 = np.genfromtxt("data/GaAsDotiert2.txt", unpack = True)
w3 = g3 + s3/60
mw3 = [0] * 18

n = 1
while(n < 19):
    mw1[n-1] = (w1[3*n-3] + w1[3*n-2] + w1[3*n-1]) / 3
    mw2[n-1] = (w2[3*n-3] + w2[3*n-2] + w2[3*n-1]) / 3
    mw3[n-1] = (w3[3*n-3] + w3[3*n-2] + w3[3*n-1]) / 3
    n += 1

bw1 = [0] * 18
bw2 = [0] * 18
bw3 = [0] * 18

n = 0
while(n < 18):
    bw1[n] = mw1[n] * 2 * 3.1416 / 360
    bw2[n] = mw1[n] * 2 * 3.1416 / 360
    bw3[n] = mw1[n] * 2 * 3.1416 / 360
    n += 1

wd1 = [0] * 9
wd2 = [0] * 9
wd3 = [0] * 9

n = 1
while(n < 10):
    wd1[n-1] = (0.5 * ( bw1[2*n-2] - bw1[2*n-1] ) / 5110) * 10000
    wd2[n-1] = (0.5 * ( bw2[2*n-2] - bw2[2*n-1] ) / 1360) * 10000
    wd3[n-1] = (0.5 * ( bw3[2*n-2] - bw3[2*n-1] ) / 1296) * 10000
    n += 1

lam = [1.06, 1.29, 1.45, 1.72, 1.96, 2.156, 2.34, 2.51, 2.65]
lam2 = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

n = 0
while(n < 9):
    lam2[n] = lam[n] ** 2
    n += 1

plt.figure()
plt.plot(lam2, wd1, 'x', label='Erste Probe')
plt.plot(lam2, wd2, 'x', label='Zweite Probe')
plt.plot(lam2, wd3, 'x', label='Dritte Probe')
plt.xlabel(r'$\lambda^2 [\mu m^{2}]$')
plt.ylabel(r'$\theta [\frac{rad \cdot 10^{-5}}{\mu m}]$')
plt.legend(loc='best')
plt.grid()
plt.savefig('build/plot-Werte.pdf')


lw1 = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
lw2 = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
n = 0
while(n < 9):
    lw1[n] = wd2[n] - wd1[n]
    lw2[n] = wd3[n] - wd1[n]
    n += 1

plt.figure()
plt.plot(lam2, lw1, 'x', label='Zweite Probe')
plt.plot(lam2, linregress(lam2, lw1), label='1.lineare Regression')
plt.plot(lam2, lw2, 'x', label='Dritte Probe')
plt.plot(lam2, linregress(lam2, lw2), label='2.lineare Regression')
plt.xlabel(r'$\lambda^2 [\mu m^{2}]$')
plt.ylabel(r'$\theta [\frac{rad \cdot 10^{-5}}{\mu m}]$')
plt.legend(loc='best')
plt.grid()
plt.savefig('build/plot-linReg.pdf')