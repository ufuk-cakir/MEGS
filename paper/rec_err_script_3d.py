import numpy as np
import matplotlib.pyplot as plt
from tqdm import trange

import sys
sys.path.insert(0, '../src')

from megs.data import image, DataLoader, Galaxy






IMG_SHAPE = (64,64,64)

norm = image.norm # Normalization function
lower = 0.25
upper = 1.0
_norm_function_kwargs = {"Masses": {"takelog": True, "plusone": True,"lower": lower, "upper": upper},
                     "GFM_Metallicity": {"takelog": True, "plusone": True, "lower": 0.25, "upper": upper},
                     "GFM_StellarFormationTime": {"takelog": True, "plusone": True, "lower": 0.25, "upper": upper},
                    
}



data = DataLoader("/export/home/ucakir/MEGS/MEGS/src/megs/data/galaxy_data.hdf5", m_min = 8)

scores = np.load('scores.npy')
eigen3d = np.load("eigengalaxies.npy")
means = np.load("means.npy").reshape(-1)
eigengalaxies = eigen3d.reshape(512,-1)



IMG_SHAPE = (3,64,64,64)
def reconstruct(index,ncomp = 215):
    score = scores[index][:ncomp]
    eig = eigengalaxies[:ncomp]
    reconstructed = np.dot(score, eig) + means
    reconstructed = reconstructed.reshape(IMG_SHAPE)
    return reconstructed

def get_original(index):
    metal = data.get_image("stars", "GFM_Metallicity", index, dim = 3).flatten()
    age = data.get_image("stars", "GFM_StellarFormationTime", index, dim = 3).flatten()
    mass = data.get_image("stars", "Masses", index, dim = 3).flatten()
    metal = norm(metal,**_norm_function_kwargs["GFM_Metallicity"])
    age = norm(age,**_norm_function_kwargs["GFM_StellarFormationTime"])
    mass = norm(mass,**_norm_function_kwargs["Masses"])
    
    return np.concatenate([metal,age,mass])
def sse(reconstructed, original):
    diff = np.sum((original - reconstructed.reshape(-1)) ** 2)
    sse = diff/np.sum(original)
    return sse



num_gal = data.get_attribute("mass").shape[0]
from tqdm import trange
import os
def rec_error_loop(ncomp = 215):
    os.makedir("rec_err_3d",exist_ok=True)
    num_gal = data.get_attribute("mass").shape[0]
    RE = []
    for i in trange(num_gal):
        s = sse(reconstructed=reconstruct(i,ncomp=ncomp),original=get_original(i))
        RE.append(s)
    np.save(f'rec_err_3d/reconstruction_err_{ncomp}_3d.npy',RE)
    return RE


if __name__ == "__main__":
    RE = []
    for i in trange(num_gal):
        s = sse(reconstructed=reconstruct(i,ncomp=215),original=get_original(i))
        RE.append(s)

