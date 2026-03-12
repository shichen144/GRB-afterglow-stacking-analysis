import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def joint(x, a, b, c, d, e, f):
    return a * np.power(x, b) * np.power((x/c)**5 + 1, -d/5) * np.power((x/e)**5 + 1, -f/5)

lcxt = [10**((i*0.1-1)) for i in range(70)]
lcx = np.array(lcxt)
lcx1 = np.arange(0.1, 2e5, 0.1)

pc = [0, 0.1, 1, 0.1, 100, 0.1]
popt, pcov = curve_fit(joint, ttime[:11], LC_flux[0][0:11], 
                       sigma=LC_flux_err[0][:11], absolute_sigma=True, maxfev=100000)
perr = np.sqrt(np.diag(pcov))

print("Parameters:", popt)
print("Errors:", perr)

LC3 = joint(lcx, *popt)
LC2 = joint(lcx1, *popt)

plt.figure(figsize=(6, 4), dpi=100)
plt.plot(lcx1, LC2, color='grey', alpha=0.8, linestyle='--')
plt.plot(lcx, LC3, 'k-')
plt.errorbar(ttime[:12], LC_flux[0], yerr=LC_flux_err[0], 
             fmt='b+', alpha=0.8, label='100MeV - 1GeV', capsize=2)

plt.xlabel('Time [s]', fontweight='bold', fontsize=12)
plt.ylabel(r'Flux [erg cm$^{-2}$ s$^{-1}$]', fontweight='bold', fontsize=12)

alpha1_err = np.sqrt(perr[1]**2 + perr[5]**2 - 2*pcov[1,5])
alpha2_err = np.sqrt(perr[1]**2 + perr[3]**2 + perr[5]**2 - 
                     2*pcov[1,3] - 2*pcov[1,5] + 2*pcov[3,5])

plt.text(0.25, 0.4, r'$\alpha_1$ = ' + f'{-popt[1]+popt[5]:.2f} $\pm$ {alpha1_err:.2f}', 
         horizontalalignment='center', transform=plt.gca().transAxes, fontsize=12)
plt.text(0.25, 0.3, r'$\alpha_2$ = ' + f'{-popt[1]+popt[3]+popt[5]:.2f} $\pm$ {alpha2_err:.2f}', 
         horizontalalignment='center', transform=plt.gca().transAxes, fontsize=12)
plt.text(0.32, 0.2, r'$t_{break}$ = ' + f'{popt[2]:.2f} $\pm$ {perr[2]:.2f} s', 
         horizontalalignment='center', transform=plt.gca().transAxes, fontsize=12)

plt.xlim(0.1, 100000)
plt.ylim(2e-10, 2e-4)
plt.tick_params(axis='both', labelsize=12)
plt.xscale("log", base=10)
plt.yscale("log", base=10)
plt.legend(fontsize=12)
plt.show()
