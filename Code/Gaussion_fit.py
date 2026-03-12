import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad

radius = 5
bin_size = 10
wbin = 0.4
pc = [800000, 0.3]
binsz = 0.1
x_ra = 100
y_dec = 100

psf_data = result2 - result1

psf1 = []
error1 = []
omega = []
theta = []

for i in range(int(radius / wbin)):
    counts1 = []
    for x in range(psf_data.shape[0]):
        for y in range(psf_data.shape[1]):
            r2 = (x - x_ra)**2 + (y - y_dec)**2
            r_max = ((i + 1) * wbin / binsz)**2
            r_min = (i * wbin / binsz)**2
            if r2 < r_max and r2 >= r_min:
                counts1.append(psf_data[x, y])
    
    dtheta = (i + 1) * wbin
    theta.append(dtheta)
    domega = 2 * np.pi * (np.cos(i * wbin * np.pi / 180) - np.cos((i + 1) * wbin * np.pi / 180))
    omega.append(domega)
    psf1.append(sum(counts1) / domega)
    error1.append(np.sqrt(abs(sum(counts1))) / domega)

def gaussian(x, a, sigma):
    return a / (2 * np.pi * sigma) * np.exp(-0.5 * (x)**2 / sigma**2)

xx = np.linspace(0, bin_size, 1000)
popt1, pcov1 = curve_fit(gaussian, theta, psf1, p0=pc)

def f(x):
    return popt1[0] / (2 * np.pi * popt1[1]) * np.exp(-0.5 * x**2 / popt1[1]**2)

results1, error2 = quad(f, 0, radius)

x = np.arange(0, radius, 0.01)
y = []
ra = []
for i in range(len(x)):
    yy, _ = quad(f, 0, 0.01 * i)
    if yy / results1 >= 0.99998:
        y.append(i * 0.01)
        ra.append(yy / results1)

results, error = quad(f, 0, y[0])
print('ratio=', ra[0], ', theta=', y[0], ', excess=', results)

count1 = []
count_on = []
count_off = []
for i in range(result1.shape[0]):
    for j in range(result1.shape[1]):
        if (i - 100)**2 + (j - 100)**2 <= (y[0] * 10)**2:
            counts1 = result2[i, j] - result1[i, j]
            counts_on = result2[i, j]
            counts_off = result1[i, j]
            count1.append(counts1)
            count_on.append(counts_on)
            count_off.append(counts_off)

print('ex=', sum(count_on))
alpha = 1
lamda = sum(count_on) * np.log((1 + alpha) / alpha * sum(count_on) / (sum(count_on) + sum(count_off))) + \
        sum(count_off) * np.log((1 + alpha) * sum(count_off) / (sum(count_on) + sum(count_off)))
S1 = np.sqrt(2 * lamda)
S2 = sum(count1) / np.sqrt(sum(count_off))
print('Sig=', S1)
print('sig=', S2)

plt.figure(figsize=(6, 4), dpi=100)

plt.errorbar(theta, psf1, fmt='b+', yerr=error1, label='Excess 1 - 100 GeV', markersize=10, alpha=1)
plt.plot(xx, gaussian(xx, *popt1), 'c-', label=f'Gaussian $\sigma$ = {popt1[1]:.2f}')

plt.text(0.15, 0.9, 'c', horizontalalignment='center', transform=plt.gca().transAxes, fontweight='bold', fontsize=15)
plt.text(0.75, 0.7, 'The Stacked GRBs', horizontalalignment='center', transform=plt.gca().transAxes, fontsize=12, color='k')

plt.xlim(0, 5)
plt.xlabel(r'$\Delta\Theta$ [degree]', fontweight='bold', fontsize=12)
plt.ylabel('dN/d$\u03A9$ [sr$^{-1}$]', fontweight='bold', fontsize=12)
plt.tick_params(axis='both', labelsize=10)
plt.legend(fontsize=10)

plt.show()