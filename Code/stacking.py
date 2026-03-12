import numpy as np
import os
from astropy.io import fits
from scipy.optimize import curve_fit
from scipy.integrate import quad
from scipy.stats import poisson

# ==================== 1. Parameter Settings ====================
map_size = 200
binsz = 0.1  # deg/pixel
path1 = '/path/to/background/files/'
path2 = '/path/to/signal/files/'

# PSF fitting parameters
radius = 5  # fitting radius (deg)
wbin = 0.25  # bin width (deg)
pc = [8000, 0.1, 0.1, 0.1]  # initial fit parameters [N0, a, sigma1, sigma2]
x_ra, y_dec = 100, 100  # source center (pixel coordinates)

# Energy bin and time parameters
Ebin = [10**2.5-10**2, 10**3-10**2.5, 10**3.5-10**3, 
        10**4-10**3.5, 10**4.5-10**4, 10**5-10**4.5]  # MeV
tbin = 10  # time bin width (s)
Omega = 2.4  # LAT field of view (sr)

# ==================== 2. Stack FITS Files ====================
def stack_fits_files(path, map_size):
    """Stack all FITS files in the given directory"""
    result = np.zeros((map_size, map_size))
    filelist = os.listdir(path)
    for filename in filelist:
        with fits.open(path + filename) as hdu:
            result += hdu[0].data
    return result

print("Stacking background files...")
result1 = stack_fits_files(path1, map_size)  # background map
print("Stacking signal files...")
result2 = stack_fits_files(path2, map_size)  # signal map

# ==================== 3. Compute Angular Distribution (PSF) ====================
psf_data = result2 - result1  # excess count map

psf1, error1, theta = [], [], []
n_bins = int(radius / wbin)

for i in range(n_bins):
    counts1 = []
    r_min = (i * wbin / binsz) ** 2
    r_max = ((i + 1) * wbin / binsz) ** 2
    
    for x in range(psf_data.shape[0]):
        for y in range(psf_data.shape[1]):
            r2 = (x - x_ra)**2 + (y - y_dec)**2
            if r_min <= r2 < r_max:
                counts1.append(psf_data[x, y])
    
    dtheta = (i + 1) * wbin
    theta.append(dtheta)
    domega = 2 * np.pi * (np.cos(i * wbin * np.pi / 180) - 
                          np.cos((i + 1) * wbin * np.pi / 180))
    psf1.append(sum(counts1) / domega)
    error1.append(np.sqrt(abs(sum(counts1))) / domega)

# ==================== 4. Fit PSF ====================
def gaussian(x, N0, a, sigma1, sigma2):
    """Double Gaussian PSF model"""
    return N0 * (a / (2 * np.pi * sigma1) * np.exp(-0.5 * (x)**2 / sigma1**2) + 
                 (1 - a) / (2 * np.pi * sigma2) * np.exp(-0.5 * (x)**2 / sigma2**2))

popt1, _ = curve_fit(gaussian, theta, psf1, p0=pc)

def psf_integral(x):
    """Angular integral of PSF (for containment fraction calculation)"""
    return 2 * np.pi * (1 - np.cos(x * np.pi / 180)) * popt1[0] * (
        popt1[1] / (2 * np.pi * popt1[2]) * np.exp(-0.5 * x**2 / popt1[2]**2) +
        (1 - popt1[1]) / (2 * np.pi * popt1[3]) * np.exp(-0.5 * x**2 / popt1[3]**2)
    )

# Compute radius containing 90% of counts
total_int, _ = quad(psf_integral, 0, radius)
theta_90 = None
for r in np.arange(0, radius, 0.01):
    partial_int, _ = quad(psf_integral, 0, r)
    if partial_int / total_int >= 0.9:
        theta_90 = r
        break

print(f'Radius containing 90% of counts: {theta_90:.2f} deg')

# ==================== 5. Extract Source Region Counts ====================
source_counts = []
source_radius_pix = theta_90 / binsz

for x in range(result1.shape[0]):
    for y in range(result1.shape[1]):
        if (x - x_ra)**2 + (y - y_dec)**2 <= source_radius_pix**2:
            source_counts.append(result2[x, y] - result1[x, y])

excess_counts = sum(source_counts)
print(f'Excess counts in source region: {excess_counts:.1f}')

# ==================== 6. Compute Flux ====================
# Note: The following arrays need to be defined based on your actual data:
# theta_c(Parameters related to the effective are), excess, counts, theta_err, etc.

# Example structure (uncomment and modify with your actual data):
# for i in range(6):
#     flux = (theta_c[i] * ex[i] / cou[i]) / Ebin[i] / tbin / 1.6 * 1e6
#     error = np.sqrt(theta_err[i]**2 + (0.05 * theta_c[i] * ex[i] / cou[i])**2) / Ebin[i] / tbin / 1.6 * 1e6

# ==================== 7. Significance Evaluation ====================
def poisson_cdf(k, lam):
    """Poisson cumulative distribution function"""
    return sum(poisson.pmf(np.arange(k+1), lam))

# Significance calculation would go here using your actual data

print("Processing complete!")
