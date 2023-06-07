import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import uncertainties as unp
from scipy import stats
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)


tem = pd.read_csv('data/TEM-Moden.txt', delimiter=';', header=0 )
pol = pd.read_csv('data/Polarisation.txt', delimiter=';', header=0 )
tem.columns = ['Position', 'I_Mode1', 'I_Mode0'] 
pol.columns = ['Winkel', 'Intensitaet']


def schoenerPlot():
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)


def model1(x, I0, x0, w):
    return I0 * np.exp(-(x-x0)**2/(2*w**2))


def model2(x, I1, x0, x1, w):
    return I1 * (8*(x-x0)**2/w**2) * np.exp(-(x-x1)**2/(2*w**2))


def poli(x, I, phi0):
    return (I * (np.sin(np.radians(x)-phi0)**2))


# Plotte/Fitte die nullte Mode
x = np.linspace(-6,15,1000)
parameters1, pcov1 = curve_fit(model1, tem['Position'] , tem['I_Mode0'], sigma=None)
I0 = unp.uarray(parameters1[0],pcov1[0,0])
x0_1 =unp.uarray(parameters1[1],pcov1[1,1])
w_1 =unp.uarray(parameters1[2],pcov1[2,2])

# print('Parameter der Grundmode: ', '\n', 'I = ', I0, '\n', 'x0_1 = ', x0_1, '\n', 'w_1 = ', w_1, '\n')

plt.plot(x,model1(x, *parameters1), label="Fit", color='red')
plt.plot(tem['Position'], tem['I_Mode0'], 'x', label='Mode')
# plt.xlabel(r'$d \,/\,\unit{\milli\metre}$')
# plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
plt.xlim(-6,15)
plt.xticks(np.linspace(-6,15,7),np.linspace(-6,15,7))
schoenerPlot()
plt.savefig("build/plot.pdf")
plt.clf()


# Plotte/Fitte die erste Mode
tem['Position'] -=5

x = np.linspace(-15,10,10000)
parameters2, pcov2 = curve_fit(model2, tem['Position'], tem['I_Mode1'], sigma=None)
I1= unp.uarray(parameters2[0], pcov2[0,0])
x0_2 =unp.uarray(parameters2[1],pcov2[1,1])
x1 =unp.uarray(parameters2[2],pcov2[2,2])
w_2 =unp.uarray(parameters2[2],pcov2[2,2])
# print(r'Parameter der Mode_01: ', '\n', 'I1 = ', I1, '\n', 'x0_2 = ', x0_2, '\n', 'x1 = ', x1, '\n', 'w_2 = ', w_2, '\n')


plt.plot(x,model2(x, *parameters2), label='Fit', color='red')
plt.plot(tem['Position'], tem['I_Mode1'], 'x', label='Mode')
# plt.xlabel(r'$d \,/\,\unit{\milli\metre}$')
# plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
plt.xticks([-15,-10,-5,0,5,10], [-10,-5,0,5,10,15])
plt.xlim(-15,10)
schoenerPlot()
plt.savefig("build/plot1.pdf")
plt.clf()

# Polarisation
x = np.linspace(0,360,1000)
xticks = np.linspace(0,360,5)
parameters3, pcov3 = curve_fit(poli, pol['Winkel'], pol['Intensitaet'], sigma=None)
I_pol =unp.uarray(parameters3[0],pcov3[0,0])
phi_0 =unp.uarray(np.degrees(parameters3[1]),np.degrees(pcov3[1,1]))
# print(r'Parameter der Polarisation: ', '\n', 'I_pol = ', I_pol, '\n', 'phi_0 = ', phi_0)


plt.plot(pol['Winkel'], pol['Intensitaet'], 'x', label='Werte')
plt.plot(x, poli(x, *parameters3), 'r', label='Fit')
plt.xticks(xticks, ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
plt.xlim(0,360)
# plt.xlabel(r'$\varphi \,/\,\unit{\radian}$')
# plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
schoenerPlot()
plt.savefig("build/plot2.pdf")
plt.clf()




# ----------------------------
####### Berechnungen #########


def wellenlaenge(b, n, d, L):
    return ( b/n * np.sin(np.arctan(d/(L))) )


gitter1 = np.array([2.4, 5.0, 7.7, 10.5, 13.1, 16.0, 19.0, 22.0])
gitter2 = np.array([3.2, 6.4, 9.9, 13.2, 16.9, 20.8, 24.9, 29.7])
gitter3 = np.array([8.5, 25.3])
gitter4 = 25.2

# berechne die Wellenlängen des ersten Gitters und trage sie im entsprechenden array ein
b = 1/80 * 10**(-3)
L = 63.4 * 10**(-2)
wellenlaenge1 = np.ones(len(gitter1))
n = 1
for i in gitter1:
    wellenlaenge1[n-2] = wellenlaenge(b,n, i*10**(-2), L)
    n += 1


# berechne die Wellenlängen des zweiten Gitters und trage sie im entsprechenden array ein
b = 1/100 * 10**(-3)
wellenlaenge2 = np.ones(len(gitter2))
n = 1
for i in gitter2:
    wellenlaenge2[n-2] = wellenlaenge(b,n, i*10**(-2), L)
    n += 1

# berechne die Wellenlängen des dritten Gitters und trage sie im entsprechenden array ein
b = 1/600 * 10**(-3)
wellenlaenge3 = np.ones(len(gitter3))
L = 47.75 *10*(-2)
n = 1
for i in gitter3:
    wellenlaenge3[n-2] = wellenlaenge(b,n, i*10**(-2), L)
    n += 1

# berechne die Wellenlängen des vierten Gitters und trage sie im entsprechenden array ein
b = 1/1200 * 10**(-3)
wellenlaenge4 = wellenlaenge(b,1, gitter4*10**(-2), L)

alleMeans = np.concatenate((wellenlaenge1, wellenlaenge2))

### Resionatorlänge bestimmen

c = 2.99792458*10**8
def laenge(f):
    return (c/(2*f*10**(6)))



print('Wellenlaenge 1:', wellenlaenge1, '\n', 'Mittelwert1:', np.mean(wellenlaenge1), '\n')
print('Wellenlaenge 2:', wellenlaenge2, '\n', 'Mittelwert2:', np.mean(wellenlaenge2), '\n')
print('Wellenlaenge 3:', wellenlaenge3, '\n', 'Mittelwert3:', np.mean(wellenlaenge3), '\n')
print('Wellenlaenge 4:', wellenlaenge4, '\n')
print('Mittelwert aller Verteilungen: ', np.mean(alleMeans))
print('Resonatorlängen:', '\n', laenge(300), '\n', laenge(200), '\n', laenge(150), '\n', laenge(125))
print('Abweichung der Resonatorlänge:', '\n', (laenge(300)/50)*100, '\n', (laenge(200)/75)*100, '\n', (laenge(150)/100)*100, '\n', (laenge(125)/125)*100)
# print(tem['I_Mode0'])