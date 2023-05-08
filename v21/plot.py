import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as const
import scipy.optimize as op


# Daten und Konstanten importieren / definieren
x = np.linspace(0.1,1,1000)

data1 = pd.read_csv('data/data1.txt', sep=" ", header=None)
data1.columns = ['f', 'B_sweep', 'B_horizontal']

data2 = pd.read_csv('data/data2.txt', sep=" ", header=None)
data2.columns = ['f', 'B_sweep', 'B_horizontal']

mu_0 = const.mu_0
R_sweep = 0.1639        # Radius der sweep Spule / m
N_sweep = 11            
R_horizontal = 0.1579   #Radius der horizontalen Helmholtz Spule /m
N_horizontal = 154
R_vertikal = 0.11735    #Radius der vertikalen Helmholtz Spule /m
N_vertiakl = 20

### Daten in die richtige Einheit bringen
def B(I, n, r): #Betrag der magnetischen Flussdichte in der Mitte einer Helmholtz-Spule mit Radius r, Windungszahl N und Stromstärke I
    return ((mu_0*8*n*I)/(np.sqrt(125)*r))


data1['B_horizontal'] *= 200 * 10**(-3) # von 200 Milliampere auf Ampere umrechnen
data2['B_horizontal'] *= 200 * 10**(-3)


I_1 = np.array(data1['B_sweep'] + data1['B_horizontal'])    #Strom des gesamten B-Feldes (nur der Strom, nicht Flussstärke)
I_2 = np.array(data2['B_sweep'] + data2['B_horizontal'])

B_1 = B(data1['B_sweep'], N_sweep, R_sweep) + B(data1['B_horizontal'], N_horizontal, R_horizontal)    #Berechne die Magnetfeldstärke aus den beiden Spulen
B_2 = B(data2['B_sweep'], N_sweep, R_sweep) + B(data2['B_horizontal'], N_horizontal, R_horizontal)


def f(x,m,b):   #Funktion für die lineare Ausgleichsgerade
    return (m*x+b)


# Daten plotten
# Daten des ersten Isotopes
params, pcov = op.curve_fit(f, data1['f'], B_1)
plt.plot(data1['f'], B_1, 'bx',label='Isotop 1')
plt.plot(x, f(x, *params), "b--", label="lineare Regression")

plt.xlim(0.1,1) #Plot schoener machen
plt.xlabel(r'Frequenz / $\si{\kilo\hertz}$')
plt.ylabel(r'magn. Flussdichte / $\si{\tesla}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Isotop1.pdf')
plt.close()


# Daten des zweiten Isotopes
params, pcov = op.curve_fit(f, data2['f'], B_2) 
plt.plot(data2['f'], B_2, 'rx',label='Isotop 2')
plt.plot(x, f(x, *params), "r--", label="lineare Regression")

plt.xlim(0.1,1) #Plot schoener machen
plt.xlabel(r'Frequenz / $\si{\kilo\hertz}$')
plt.ylabel(r'magn. Flussdichte / $\si{\tesla}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Isotop2.pdf')