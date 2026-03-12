import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def power_law_e(x, a, b, c):
    return a * (x/100)**(-b) * np.exp(-x/c)

xe1 = np.arange(100, 100000, 100)
pc = [1e1, 1, 10000]
popt, pcov = curve_fit(power_law_e, Ex[:6], flux[:6], p0=pc, sigma=error[:6])
perr = np.sqrt(np.diag(pcov))

ye = power_law_e(xe1, *popt)
ye_up = power_law_e(xe1, *(popt + perr))
ye_low = power_law_e(xe1, *(popt - perr))

print("Parameters:", popt)
print("Errors:", perr)

plt.figure(figsize=(6, 4), dpi=100)
plt.plot(xe1, power_law_e(xe1, *popt), 'b--', label='PLEC fit')
plt.fill_between(xe1, ye_low, ye_up, color='gray', alpha=0.5)
plt.errorbar(Ex[:6], flux[:6], fmt='bo', xerr=[Exerrl[:6], Exerru[:6]], 
             yerr=error[:6], label='Data')
plt.xlabel('E [MeV]', fontsize=12)
plt.ylabel(r'dN/dE [erg$^{-1}$ cm$^{-2}$ s$^{-1}$]', fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.legend(loc='upper right')
plt.show()