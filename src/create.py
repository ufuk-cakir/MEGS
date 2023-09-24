
from megs.data import Gamma 
import numpy as np 
import h5py

halo_ids = np.load("/export/home/ucakir/MEGS/MEGS/src/gamma_ids_outlier_removed.npy")
outlier_ids = np.load("/export/home/ucakir/MEGS/MEGS/src/outlier_ids.npy")

old_datapath = "/export/home/ucakir/MEGS/MEGS/src/megs/data/galaxy_data.hdf5"
new_gamma_path = "/export/data/ucakir/final_gamma_data"
import os
os.makedirs(new_gamma_path, exist_ok=True)

old_data = Gamma(old_datapath)
ids_in_old = old_data.get_attribute("halo_id")
gamma_ids = np.load("/export/home/ucakir/MEGS/MEGS/src/megs/gamma_ids.npy")
print(gamma_ids.shape)
mask = np.isin(ids_in_old, gamma_ids)
# filter out the outliers
mask = np.logical_and(mask, ~np.isin(ids_in_old, outlier_ids))

old_data.mask = mask
old_data.get_attribute("halo_id").shape

new_halo_ids = old_data.get_attribute("halo_id")


from tqdm import tqdm
FIELDS = ["GFM_Metallicity", "GFM_StellarFormationTime", "Masses"]
attributes = ["mass", "halo_id"]

particle = "stars"
from megs.data.generate import _create_data_structure

# Create the data structure
_create_data_structure(
    n_galaxies=len(new_halo_ids),
    image_res = 64,
    galaxy_parameters=attributes,
    particle_types = ["stars"],
    fields = FIELDS,
    path = new_gamma_path,
)

# OPen new file
# overvrite the old file


with h5py.File(f"{new_gamma_path}/galaxy_data.hdf5","a") as f:
    for parameter in f["Galaxies/Attributes"].keys():
        print(parameter)
        f["Galaxies/Attributes"][parameter][:] = old_data.get_attribute(parameter)
    
    for d in tqdm(f["Galaxies"]["Particles"][particle]["Images"].keys()):
        for field in tqdm(f["Galaxies"]["Particles"][particle]["Images"][
                        d
                    ].keys()):
                        print(field)
                        # Get the image
                        dim = int(d[-1]) 
                        print(dim)
                        f["Galaxies"]["Particles"][particle]["Images"][d][field][:] = old_data.get_image(particle, field, dim = dim)
                        
                        
    