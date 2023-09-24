from megs.data.image import norm
import h5py
from tqdm import tqdm

lower = 0.25
upper = 1.0
norm_function_args = {"Masses": {"takelog": True, "plusone": True,"lower": lower, "upper": upper},
                     "GFM_Metallicity": {"takelog": True, "plusone": True, "lower": 0.25, "upper": upper},
                     "GFM_StellarFormationTime": {"takelog": True, "plusone": True, "lower": 0.25, "upper": upper},
}

datapath = "/export/data/ucakir/final_gamma_data/galaxy_data.hdf5"

particle = "stars"

# Open the HDF5 file in append mode
with h5py.File(datapath, "a") as f:
    for d in tqdm(f["Galaxies"]["Particles"][particle]["Images"].keys()):
        for field in tqdm(f["Galaxies"]["Particles"][particle]["Images"][d].keys()):
            print(field)
            # Get the image
            dim = int(d[-1])
            print(dim)


            for index in tqdm(range(f["Galaxies"]["Particles"][particle]["Images"][d][field].shape[0])):

                image =    f["Galaxies"]["Particles"][particle]["Images"][d][field][index]
                
                # Get the max and min values
                max_value = image.max()
                min_value = image.min()
                
                # Store the max and min values
                f["Galaxies/ImageAttributes"][field]["max"][index] = max_value
                f["Galaxies/ImageAttributes"][field]["min"][index] = min_value
                
                
                # Normalize the image
                normalized_image = norm(image, **norm_function_args[field])

                # Save the normalized image back to the HDF5 file
                f["Galaxies"]["Particles"][particle]["Images"][d][field][index] = normalized_image
