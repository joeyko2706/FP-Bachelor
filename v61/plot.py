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


tem = pd.read_csv('data/TEM-Moden.txt', delimiter=';', header=1)
pol = pd.read_csv('data/Polarisation.txt', delimiter=';', header=1)
stabi = pd.read_csv('data/stabilitaetsbedingung.txt', delimiter=';', header=1)
tem.columns = ['Position', 'I_Mode1', 'I_Mode0'] 
pol.columns = ['Winkel', 'Intensitaet']
stabi.columns = ['Abstand', 'Intensitaet']


def schoenerPlot():
    plt.grid()
    plt.legend(loc='best')
    plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)


def model1(x, I0, x0, w):
    return I0 * np.exp(-(x-x0)**2/(2*w**2))


def model2(x, I1, x0, x1, w):
    return I1 * (8*(x-x0)**2/w**2) * np.exp(-(x-x1)**2/(2*w**2))


def poli(x, I):
    return (I * (np.sin(np.radians(x))**2))


# Plotte/Fitte die nullte Mode
x = np.linspace(-10,15,1000)
parameters1, pcov1 = curve_fit(model1, tem['Position'] , tem['I_Mode0'], sigma=None)
I0 = unp.uarray(parameters1[0],pcov1[0,0])
x0_1 =unp.uarray(parameters1[1],pcov1[1,1])
w_1 =unp.uarray(parameters1[2],pcov1[2,2])

plt.plot(x,model1(x, *parameters1), label="Fit", color='red')
plt.plot(tem['Position'], tem['I_Mode0'], 'x', label='Mode')
plt.xlabel(r'$d \,/\,\unit{\milli\metre}$')
plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
plt.xlim(-10,15)
schoenerPlot()
plt.savefig("build/plot.pdf")
plt.clf()


# Plotte/Fitte die erste Mode
tem['Position'] -=5
x = np.linspace(-15,10,1000)
parameters2, pcov2 = curve_fit(model2, tem['Position'], tem['I_Mode1'], sigma=None)
I1= unp.uarray(parameters2[0], pcov2[0,0])
x0_2 =unp.uarray(parameters2[1],pcov2[1,1])
x1 =unp.uarray(parameters2[2],pcov2[2,2])
w_2 =unp.uarray(parameters2[2],pcov2[2,2])

plt.plot(x,model2(x, *parameters2), label='Fit', color='red')
plt.plot(tem['Position'], tem['I_Mode1'], 'x', label='Mode')
plt.xlabel(r'$d \,/\,\unit{\milli\metre}$')
plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
plt.xlim(-15,10)
schoenerPlot()
plt.savefig("build/plot1.pdf")
plt.clf()

# Polarisation
x = np.linspace(0,360,1000)
parameters3, pcov3 = curve_fit(poli, pol['Winkel'], pol['Intensitaet'], sigma=None)
xticks = np.linspace(0,360,5)

plt.plot(pol['Winkel'], pol['Intensitaet'], 'x', label='Werte')
plt.plot(x, poli(x, *parameters3), 'r', label='Fit')
plt.xticks(xticks, ['0', r'$\frac{\pi}{4}$', r'$\frac{\pi}{2}$', r'$\frac{3\pi}{2}$', r'$2\pi$'])
plt.xlabel(r'$\varphi \,/\,\unit{\radian}$')
plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
plt.xlim(0,360)
schoenerPlot()
plt.savefig("build/plot2.pdf")
plt.clf()
