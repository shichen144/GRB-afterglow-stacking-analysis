from scipy import stats, optimize
import numpy as np
import matplotlib.pyplot as plt

#read all Background TS values
GRBTS = open('tsvalues220.txt', 'r')
ABts = []
lat_lines = GRBTS.readlines()
for line in lat_lines:
    a = line.strip().split(',')
    if float(a[1]) > 0:
        ABts.append(float(a[1]))
GRBTS.close()

logbins = np.arange(10, 70, 4.7)
counts, bin_edges = np.histogram(ABts, bins=logbins)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
bin_widths = (bin_edges[1:] - bin_edges[:-1])
print(sum(counts))

counts1, bin_edges1 = np.histogram(ABts, bins=np.arange(10, 100, 1))
bin_centers1 = (bin_edges1[:-1] + bin_edges1[1:]) / 2

def chi2(params):
    nu, scale = params
    expected1 = []
    for i in range(len(bin_centers)):
        cdf_low = stats.chi2.cdf(bin_edges[i], nu) * scale
        cdf_high = stats.chi2.cdf(bin_edges[1 + i], nu) * scale
        prob = cdf_high - cdf_low
        expected1.append(149 * prob)
    residuals = counts - expected1
    expected_safe = np.maximum(expected1, 0.1)
    return np.sum(residuals**2 / expected_safe)

def ts_cdf(x, nu, scale):
    if nu <= 0 or scale <= 0:
        return np.inf
    expected = []
    for i in range(len(bin_centers1)):
        cdf_low = stats.chi2.cdf(bin_edges1[i], nu) * scale
        cdf_high = stats.chi2.cdf(bin_edges1[1 + i], nu) * scale
        prob = cdf_high - cdf_low
        expected.append(149 * prob * 3.7)
    return expected

initial_p = [2, 0.5]
bounds = [(0.1, 40), (0.1, 10)]
result = optimize.minimize(chi2, initial_p, bounds=bounds)
nu_fit, scale_fit = result.x
print(nu_fit, scale_fit)

xt = np.arange(10, 70, 1)

plt.figure(figsize=(6, 4), dpi=100)
plt.hist(ABts, bins=logbins, color='k', histtype='step', label='Background')
plt.plot(bin_centers1, ts_cdf(xt, nu_fit, scale_fit), label=f'$\chi_v^2$  v={nu_fit:.2f}')
plt.text(0.75, 0.5, r'220 LAT detected GRBs', horizontalalignment='center', 
         style='italic', transform=plt.gca().transAxes, fontsize=12)
plt.text(0.05, 0.9, 'a', horizontalalignment='center', 
         transform=plt.gca().transAxes, fontweight='bold', fontsize=15)

plt.xlabel(r'TS value', fontweight='bold', fontsize=12)
plt.ylabel('Numbers', fontweight='bold', fontsize=12)
plt.xlim(10, 100)
plt.legend()
plt.show()
