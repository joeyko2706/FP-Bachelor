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


parameters1, pcov1 = curve_fit(model1, tem['Position'] , tem['I_Mode0'], sigma=None)
I0 = unp.uarray(parameters1[0],pcov1[0,0])
x0_1 =unp.uarray(parameters1[1],pcov1[1,1])
w_1 =unp.uarray(parameters1[2],pcov1[2,2])
x = np.linspace(-10,15,1000)

plt.scatter(tem['Position'], tem['I_Mode0'], marker='x', label='Mode')
plt.plot(x,model1(x, *parameters1), label="Fit", color='red')
# plt.xlabel(r'$d \,/\,\unit{\milli\metre}$')
# plt.ylabel(r'$I \,/\,\unit{\micro\watt}$')
schoenerPlot()
plt.savefig("build/plot.pdf")
