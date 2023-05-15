import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants as const
import scipy.optimize as op
from uncertainties import ufloat
from sympy import symbols, solve
import uncertainties.unumpy as unp

"""

In den Daten gab es ein paar Ausreißer, die ich nach der ersten Korrektur rausgenommen habe.
Ich habe sie in den Daten einfach rausgenommen und unten einmal definiert.

"""
# Daten und Konstanten importieren / definieren
x = np.linspace(0.1*10**6,1*10**6,1000)

data1 = pd.read_csv('data/data1.txt', sep=" ", header=None)
data1.columns = ['f', 'B_sweep', 'B_horizontal']

data2 = pd.read_csv('data/data2.txt', sep=" ", header=None)
data2.columns = ['f', 'B_sweep', 'B_horizontal']

mu_0 = const.mu_0
R_sweep = 16.39 *10**(-2)        # Radius der sweep Spule / m
N_sweep = 11            
R_horizontal = 15.79 *10**(-2)  #Radius der horizontalen Helmholtz Spule /m
N_horizontal = 154
R_vertikal = 11.735 *10**(-2)    #Radius der vertikalen Helmholtz Spule /m
N_vertiakl = 20
h = const.Planck
mu_b = const.value('Bohr magneton')


### Daten in die richtige Einheit bringen

#Betrag der magnetischen Flussdichte in der Mitte einer Helmholtz-Spule mit Radius r, Windungszahl N und Stromstärke I
def B(I, n, r): 
    return ((mu_0*8*n*I)/(np.sqrt(125)*r))


data1['f'] *= 10**6
data2['f'] *= 10**6
data1['B_horizontal'] *=  10**(-3) # von 200 Milliampere auf Ampere umrechnen
data2['B_horizontal'] *=  10**(-3)

#Berechne die Magnetfeldstärke aus den beiden Spulen 
B_1 = B(data1['B_sweep'], N_sweep, R_sweep) + B(data1['B_horizontal'], N_horizontal, R_horizontal) 
B_2 = B(data2['B_sweep'], N_sweep, R_sweep) + B(data2['B_horizontal'], N_horizontal, R_horizontal) 
B_vertikal = B(235, N_vertiakl, R_vertikal)


# Ausreißer definieren
def ausreisser(sweep, horizontal):
    horizontal *= 10**(-3)
    return (B(sweep, N_sweep, R_sweep) + B(horizontal, N_horizontal, R_horizontal))


Isof1 = 0.8 * 10 **6                             # Ausreißer des ersten Isotopes -> Frequenz
Isof2 = np.array([0.5 * 10 **6, 0.8 * 10 **6])  # des zweiten

IsoB1 = ausreisser(0.909, 154.7)
IsoB2 = np.array([ausreisser(0.065, 77.7), ausreisser(0.476, 126.6)])


# print('Ausreißer1: ', Isof1, IsoB1)
# print('Ausreißer2: ', Isof2, IsoB2)
# print('Richtige Werte1: ', data1.head())
# print('Richtige Werte2: ', data2.head())


def f(x,m,b):   #Funktion für die lineare Ausgleichsgerade
    return (m*x+b)


# Daten plotten     # Daten des ersten Isotopes
paramsB1, pcovB1 = op.curve_fit(f, data1['f'], B_1)
errorsB1 = np.sqrt(np.diag(pcovB1))
mB1 = ufloat(paramsB1[0], errorsB1[0])
bB1 = ufloat(paramsB1[1], errorsB1[1])
plt.plot(data1['f'], B_1, 'bx',label='Messdaten')
plt.plot(x, f(x, *paramsB1), "r--", label="lineare Regression")
plt.plot(Isof1, IsoB1, 'gx')

plt.xlim(0.1*10**6,1*10**6) #Plot schoener machen
plt.xticks(np.linspace(0.1*10**6,1*10**6,10),[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.yticks(np.linspace(0.00004, 0.00018, 8), [40, 60, 80, 100, 120, 140, 160, 180])
plt.xlabel(r'Frequenz $f$ / $\si{\mega\hertz}$')
plt.ylabel(r'magn. Flussdichte $B$ / $\si{\micro\tesla}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Isotop1.pdf')
plt.close()


# Daten des zweiten Isotopes
paramsB2, pcovB2 = op.curve_fit(f, data2['f'], B_2) 
errorsB2 = np.sqrt(np.diag(pcovB2))
mB2 = ufloat(paramsB2[0], errorsB2[0])
bB2 = ufloat(paramsB2[1], errorsB2[1])
plt.plot(data2['f'], B_2, 'bx',label='Messdaten')
plt.plot(x, f(x, *paramsB2), "r--", label="lineare Regression")
plt.plot(Isof2, IsoB2, 'gx')

plt.xlim(0.1*10**6,1*10**6) #Plot schoener machen
plt.xticks(np.linspace(0.1*10**6,1*10**6,10),[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
plt.yticks(np.linspace(0.00005, 0.000225,8), [50, 75, 100, 125, 150, 175, 200, 225])
plt.xlabel(r'Frequenz $f$ / $\si{\mega\hertz}$')
plt.ylabel(r'magn. Flussdichte $B$ / $\si{\micro\tesla}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/Isotop2.pdf')


#Lande-Faktoren bestimmen über die Steigung
g1 = h/(mu_b*mB1)
g2 = h/(mu_b*mB2)
g = g1/g2

# Kernspin bestimmen

J = 0.5
S = 0.5
L = 0

gJ = 1+(J+J*(J+1)+S*(S+1)-L*(L+1))/(2*J*(J+1))

I1 = gJ / (4 * g1) - 1 + unp.sqrt((gJ / (4 * g1) - 1)**2+ 3 * gJ / (4 * g1) - 3 / 4)
I2 = gJ / (4 * g2) - 1 + unp.sqrt((gJ / (4 * g2) - 1)**2+ 3 * gJ / (4 * g2) - 3 / 4)

# Isotopenverhältnis ausrechnen:
Oszillos1 = 0.2
Oszillos2 = 0.5
Isotopenverhältnis = Oszillos1/Oszillos2

# quadratischen Zeeman-Effekt ausrechnen

U1 = g1*mu_b*np.max(B_1)+g1**2*mu_b**2*np.max(B_1)**2*(1-2*2)/(4.53e-24) # Zeeman-Aufspaltung
U2 = g2*mu_b*np.max(B_2)+g2**2*mu_b**2*np.max(B_2)**2*(1-2*3)/(2.01e-24)

"""

print('--------------')
print('Ausgerechnete Werte: ', '\n')
print('Erdmagnetfeld: ', B_vertikal)
print('Steigung1 = ', '\t', mB1)
print('Achsenabschnitt 1 =', '\t', bB1)
print('Steigung2 = ', '\t', mB2)
print('Achsenabschnitt 2 =', '\t', bB2)
print('LandeFaktor1 = ', '\t', f'{g1:.5f}')
print('LandeFaktor2 = ', '\t', f'{g2:.5f}')
print("Verhältnis 1/2: ", g)
print('g_j = ', '\t', gJ)
print('Kernspin_1= ', '\t', f'{I1:.5f}')
print('Kernspin_2= ', '\t', f'{I2:.5f}')
print('Abweichung Kernspin1: ', f'{abs(1-I1/1.5):.5f}')
print('Abweichung Kernspin1: ', f'{abs(1-I2/2.5):.5f}')
print('Isotopenverhältnis I_1/I_2 = ', Isotopenverhältnis, '\n')
print('--------------')
print('Quadratische Zeeman-Aufspaltung')
print('Maximales BFeld1: ', np.round(np.max(B_1)*10**6, 4))
print('Maximales BFeld2: ', np.round(np.max(B_2)*10**6, 4))
print('Quadratische Zeeman-Aufspaltung 1 in eV: ', U1/const.e)
print('Quadratische Zeeman-Aufspaltung 2 in eV: ', U2/const.e)
print('--------------', '\n')


#Daten für Latex ausgeben
# print('-------------')
# for i in range(10):   
#     print(data1.iloc[i, 0], ' & ', np.round(data1.iloc[i, 1], 4), ' & ', np.round(data1.iloc[i, 2],3), ' k')
    # print(data2.iloc[i, 0], ' & ', np.round(data2.iloc[i, 1], 4), ' & ', np.round(data2.iloc[i, 2],3), ' k')

"""