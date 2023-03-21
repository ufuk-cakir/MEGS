import numpy as np

from swiftsimio.visualisation.projection import scatter as scatter2D
from swiftsimio.visualisation.volume_render import scatter as scatter3D

import matplotlib.pyplot as plt
import os 
from tqdm import trange, tqdm

def image2D(coordinates, R_half, weights, smoothing_length, plot_factor = 10, res = 64):
    
    plot_range = plot_factor*R_half
    
    x = coordinates[:,0].copy()
    y = coordinates[:,1].copy()
    
    m =  weights
    
    h = smoothing_length.copy()
    
    #Transform Particles s.t -factor*r_halfmassrad < x <factor*r_halfmassrad -> 0 < x <1
    x = x/(2*plot_range) +1/2  
    y = y/(2*plot_range) +1/2

    h = h/(2*plot_range)
    
    SPH_hist = scatter2D(x=x, y = y,h = h, m = m ,res= res)
        
    return(SPH_hist)






def clip_image(data, lower = 0.1, upper = 1.):
    """Clip image to [lower,upper] quantile.
    """    
    hist = data.copy()
    L,U = np.quantile(hist,[lower,upper])
    hist = np.clip(hist, L, U)
    return(hist)




def norm(x, takelog = True, plusone = True, clip = True, lower = 0.1, upper = 1.):
    x = np.nan_to_num(x)
    x = x+1 if plusone else x
    mask = np.where(x!=0)
    
    x[mask] = np.log10(x[mask]) if takelog else x[mask]
    
    x[mask] = clip_image(x[mask], lower = lower, upper = upper) if clip else x[mask]
    
    
    x[mask] -= x[mask].min()
    x[mask]/=x[mask].max()
    return(x)