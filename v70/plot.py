import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as op
from uncertainties import ufloat
from scipy.stats import sem


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



####################### Leckrate Turbo
t, p1, p2, p3, p4 = np.genfromtxt("data/turboleck.txt", unpack = True)

p1 *= 1000
p2 *= 1000
p3 *= 1000
p4 *= 1000

fp1 = p1 * 0.3
fp2 = p2 * 0.3
fp3 = p3 * 0.3
fp4 = p4 * 0.3

# #Rechnung:
tv = 33
v = ufloat(tv, tv * 0.1)
rp1 = ufloat(4.94, 1.48)
m1 = ufloat(3.15, 0.14)
s1 = v * m1 / rp1
rp2 = ufloat(6.95, 2.09)
m2 = ufloat(5.66, 0.18)
s2 = v * m2 / rp2
rp3 = ufloat(1.02, 0.31)
m3 = ufloat(1.05, 0.29)
s3 = v * m3 / rp3
rp4 = ufloat(2.07, 0.62)
m4 = ufloat(2.38, 0.23)
s4 = v * m4 / rp4
# print(s1, s2, s3, s4)

plt.errorbar(t, p1, yerr=fp1, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p1), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar} \cdot 10^{-3}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotturboleck1.pdf')
plt.close() ###
plt.errorbar(t, p2, yerr=fp2, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p2), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar} \cdot 10^{-3}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotturboleck2.pdf')
plt.close() ###
plt.errorbar(t, p3, yerr=fp3, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p3), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar} \cdot 10^{-3}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotturboleck3.pdf')
plt.close() ###
plt.errorbar(t, p4, yerr=fp4, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p4), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar} \cdot 10^{-3}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotturboleck4.pdf')
plt.close() ###
# ##########################################################################################################

# ####################### Evakuierung Turbo

t, p1, p2, p3, pm, sa = np.genfromtxt("data/turboevak.txt", unpack = True)

n = 0
while n < 33:
    pm[n] = (p1[n] + p2[n] + p3[n]) / 3
    sa[n] = np.sqrt((1/(2*3)) * ((p1[n] - pm[n])**2 + (p2[n] - pm[n])**2 + (p3[n] - pm[n])**2) )
    # print(pm[n], sa[n])
    n += 1

pe = 4.2e-6
fpe = 1.26e-6
p0 = 5.01e-3
fp0 = 1.50e-3
ln = np.log((pm - pe) / (p0 - pe))
fln = ln *0.1
x = 13
t1 = t[:x]
ln1 = ln[:x]
t2 = t[x:32]
ln2 = ln[x:32]

sa = pm * 0.3

n = 0
while n < 33:
    fehler = np.sqrt((sa[n] / (pm[n] - pe))**2 + (fp0 / (p0 - pe))**2 + ((pm[n] - p0) * fpe / ((pe - p0) * (pe - pm[n])))**2)
    fln[n] = np.round(fehler, 2)
    # print(t[n], ln[n], fln[n])
    n += 1

# # # Tabelle:
# # n = 0
# # while n < 33:
# #     print(t[n], ' & $', np.round(ln[n], 2), ' \pm \, ', fln[n], '$ \ ')
# #     n += 1


#Rechnung:
tv = 33
v = ufloat(tv, tv * 0.1)
m1 = ufloat(0.302, 0.012)
s1 = (v * m1) 
m2 = ufloat(0.018, 0.001)
s2 = (v * m2)
# print(s1, s2)

plt.errorbar(t[:32], ln[:32], yerr=fln[:32], fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
# plt.plot(t[:32], ln[:32], 'b.', label = 'Messwerte')
plt.plot(t1, linregress(t1, ln1), "r-", label='lineare Regression1')
plt.plot(t2, linregress(t2, ln2), "g-", label='lineare Regression2')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'  $\ln\left(\frac{p(t) - p_E}{p_0 - p_E}\right)$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.show()
plt.savefig('build/plotturboevak.pdf')
plt.close() ##
##########################################################################################################

####################### Leckrate Dreh 1

t, p1, p2, p3, pm = np.genfromtxt("data/drehleck1.txt", unpack = True)

n = 0
while n < 21:
    pm[n] = (p1[n] + p2[n] + p3[n]) / 3
    n += 1

fp1 = p1 * 0.1
fp2 = p2 * 0.1
fp3 = p3 * 0.1
fpm = pm * 0.1

n = 0
while n < 21:
    fpm[n] = np.sqrt((1/(2*3)) * ((p1[n] - pm[n])**2 + (p2[n] - pm[n])**2 + (p3[n] - pm[n])**2) )
    # print(pm[n], fpm[n])
    n += 1
# print(np.round(fpm, 2))

# #Tabelle:
# n = 0
# while n < 21:
#     print(t[n], ' & $', p1[n], ' \pm \, ', np.round(fp1[n], 2), '$ & $', p2[n], ' \pm \, ', np.round(fp2[n], 2), '$ & $', p3[n], ' \pm \, ', np.round(fp3[n], 2), '$ & $',
#           np.round(pm[n], 2), ' \pm \, ', np.round(fpm[n], 2), '$ \ ')
#     n += 1

#Rechnung:
tv = 34
v = ufloat(tv, tv * 0.1)
rp = ufloat(0.50, 0.15)
m1 = ufloat(0.012, 0.001)
s1 = v * m1 / rp
# print(s1)

plt.errorbar(t, pm, yerr=fpm, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, pm), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotdrehleck1.pdf')
# plt.show()
plt.close() ###
# ##########################################################################################################

# ####################### Leckrate Dreh 2

t, p1, p2, p3 = np.genfromtxt("data/drehleck2.txt", unpack = True)


fp1 = 4
fp2 = 4
fp3 = 4

# # # #Tabelle:
# # # n = 0
# # # while n < 21:
# # #     print(t[n], ' & $', p1[n], ' \pm \, ', np.round(fp1[n], 2), '$ & $', p2[n], ' \pm \, ', np.round(fp2[n], 2), '$ & $', p3[n], ' \pm \, ', np.round(fp3[n], 2), '$ \ ')
# # #     n += 1

#Rechnung:
tv = 34
v = ufloat(tv, tv * 0.1)
rp1 = ufloat(9.9, 4.0)
m1 = ufloat(0.289, 0.015)
s1 = v * m1 / rp1
rp2 = ufloat(50.3, 4.0)
m2 = ufloat(1.804, 0.006)
s2 = v * m2 / rp2
rp3 = ufloat(99.7, 4.0)
m3 = ufloat(3.309, 0.025)
s3 = v * m3 / rp3
# # print(s1, s2, s3)

plt.errorbar(t, p1, yerr=fp1, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p1), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.show()
plt.savefig('build/plotdrehleck2.pdf')
plt.close() ###
plt.errorbar(t, p2, yerr=fp2, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p2), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotdrehleck3.pdf')
plt.close() ###
plt.errorbar(t, p3, yerr=fp3, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1.5, capsize = 2, capthick = 1.5)
plt.plot(t, linregress(t, p3), "r-", label='lineare Regression')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'Druck $p$ / $\si{\milli\bar}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plotdrehleck4.pdf')
plt.close() ###
##########################################################################################################

####################### Evakuierung Dreh

t, p = np.genfromtxt("data/drehevak.txt", unpack = True)
fp = p * 0.1
fp[:15] = 4
pe = 1.6e-2
fpe = 0.16e-2
p0 = p[0]

ln = np.log((p - pe) / (p[0] - pe))
fln = ln * 0.1

n = 0
while n <61:
    fehler = np.sqrt((fp[n] / (p[n] - pe))**2 + (fp[0] / (p0 - pe))**2 + ((p[n] - p0) * fpe / ((pe - p0) * (pe - p[n])))**2)
    fln[n] = fehler
    # print(t[n], fln[n])
    n += 1

# #Tabelle:
# n = 0
# while n < 31:
#     print(t[n], ' & $', p[n], ' \pm \, ', np.round(fp[n], 2), '$ & ', t[n + 30], ' & $', p[n + 30], ' \pm \, ', np.round(fp[n + 30], 2), '$ \ ')
#     n += 1

# #Tabelle:
# n = 0
# while n < 31:
#     print(t[n], ' & $', np.round(ln[n], 2), ' \pm \, ', np.round(fln[n], 3), '$ & ', t[n + 30], ' & $', np.round(ln[n + 30], 2), ' \pm \, ', np.round(fln[n + 30], 3), '$ \ ')
#     n += 1

#Rechnung:
tv = 34
v = ufloat(tv, tv * 0.1)
m1 = ufloat(0.0296, 0.0004)
s1 = (v * m1) 
m2 = ufloat(0.0058, 0.0002)
s2 = (v * m2)
# print(s1, s2)

x = 23

plt.errorbar(t, ln, yerr=fln, fmt = '.', label = 'Messwerte', ecolor='black', elinewidth = 1, capsize = 1, capthick = 1)
plt.plot(t, ln, 'b.', label = 'Messwerte')
plt.plot(t[:x], linregress(t[:x], ln[:x]), "r-", label='lineare Regression1')
plt.plot(t[x:], linregress(t[x:], ln[x:]), "g-", label='lineare Regression2')
plt.xlabel(r'Zeit $t$ / $\si{\second}$')
plt.ylabel(r'  $\ln\left(\frac{p(t) - p_E}{p_0 - p_E}\right)$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
# plt.show()
plt.savefig('build/plotdrehevak.pdf')
plt.close() 



##################################################################################################################

s1 = p2
fs1 = p3
pa = t
fpa = t * 0.1

s1[0] = 10
s1[1] = 21
s1[2] = 27
s1[3] = 34
s1[4] = 38
s1[5] = 0.59

fs1[0] = 1.1
fs1[1] = 7
fs1[2] = 9
fs1[3] = 14
fs1[4] = 13
fs1[5] = 0.07

pa[0] = 2.5e-3
pa[1] = 5e-5
pa[2] = 7e-5
pa[3] = 1e-4
pa[4] = 2e-4
pa[5] = 9e-6

fpa = pa * 0.3
fpa[0] = 2.4e-3

plt.errorbar(pa[:6], s1[:6], yerr=fs1[:6], xerr=fpa[:6], fmt = '.', label = 'Saugvermögen', ecolor='black', elinewidth = 1, capsize = 1, capthick = 1)
plt.xlabel(r'Druck $p$ / $\si{\milli\bar}$')
plt.ylabel(r' Saugvermögen $S(p)$ / $\si{\litre\per\second}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/saug1.pdf')
# plt.show()
plt.close() 

s2 = p2
fs2 = p3
pb = t

s2[0] = 1.01
s2[1] = 0.197
s2[2] = 0.82
s2[3] = 0.99
s2[4] = 1.2
s2[5] = 1.1

fs2[0] = 0.11
fs2[1] = 0.021
fs2[2] = 0.09
fs2[3] = 0.32
fs2[4] = 0.4
fs2[5] = 0.4

pb[0] = 501
pb[1] = 1
pb[2] = 0.5
pb[3] = 10
pb[4] = 50
pb[5] = 100

fpb = pb * 0.1
fpb[0] = 499
fpb[2] = 4
fpb[3] = 4
fpb[4] = 4
fpb[5] = 4

plt.errorbar(pb[:6], s2[:6], yerr=fs2[:6], xerr=fpb[:6], fmt = '.', label = 'Saugvermögen', ecolor='black', elinewidth = 1, capsize = 1, capthick = 1)
plt.xlabel(r'Druck $p$ / $\si{\milli\bar}$')
plt.ylabel(r' Saugvermögen $S(p)$ / $\si{\litre\per\second}$')
plt.legend(loc='best')
plt.tight_layout(pad=0, h_pad=1.08, w_pad=1.08)
plt.savefig('build/saug2.pdf')
# plt.show()
