import numpy as np
from numpy import arange

smooth_r = 0.4
sigma_w = 0.3
binsz = 0.1

def calculate_significance_map(data_map, bkg_map, smooth_r, binsz):
    r_pix = smooth_r / binsz
    half_width = int(r_pix)
    
    nx, ny = data_map.shape
    sum_on = np.zeros_like(data_map, dtype=float)
    sum_off = np.zeros_like(data_map, dtype=float)
    sig = np.zeros_like(data_map, dtype=float)
    Sig = np.zeros_like(data_map, dtype=float)
    
    for x in range(nx):
        for y in range(ny):
            smooth_counts = []
            bkg_counts = []
            
            for i in range(-half_width, half_width + 1):
                for j in range(-half_width, half_width + 1):
                    if (x + i >= 0 and x + i < nx and 
                        y + j >= 0 and y + j < ny and 
                        i**2 + j**2 <= r_pix**2):
                        smooth_counts.append(data_map[x + i][y + j])
                        bkg_counts.append(bkg_map[x + i][y + j])
            
            sum_on[x][y] = sum(smooth_counts)
            sum_off[x][y] = sum(bkg_counts)
            
            sig[x][y] = (sum_on[x][y] - sum_off[x][y]) / np.sqrt(sum_off[x][y])
            
            if sum(smooth_counts) < 0 or sum(bkg_counts) < 0:
                lamda = (sum(smooth_counts) * np.log(2 * sum(smooth_counts) / (sum(smooth_counts) + sum(bkg_counts))) + 
                        sum(bkg_counts) * np.log(2 * sum(bkg_counts) / (sum(smooth_counts) + sum(bkg_counts))))
                if sum(smooth_counts) >= sum(bkg_counts):
                    Sig[x][y] = np.sqrt(2 * lamda)
                else:
                    Sig[x][y] = -np.sqrt(2 * lamda)
            else:
                Sig[x][y] = (sum_on[x][y] - sum_off[x][y]) / np.sqrt(sum_off[x][y])
    
    return sig, Sig

# Example usage:
# sig_gaussian, sig_lima = calculate_significance_map(result2, result1, smooth_r, binsz)