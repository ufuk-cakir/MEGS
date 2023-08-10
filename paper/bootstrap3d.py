



import sys
sys.path.insert(0, '../src')

from megs.model.mPCA import mPCA
from megs.data import image, DataLoader, Galaxy



import numpy as np
import joblib


from tqdm import trange,tqdm


particle_type = "stars"
_dim = "dim3"
dim = 3
_IMG_SHAPE = (64,64,64) if dim == 3 else (64,64)
print("Loading data...")





data = DataLoader("/export/home/ucakir/MEGS/MEGS/src/megs/data/galaxy_data.hdf5", m_min = 8)
print("Create datamatrix...")
datamatrix = np.zeros((data.get_attribute("mass").shape[0], 3 * np.prod(_IMG_SHAPE))) 


_norm_function = image.norm # Normalization function
lower = 0.25
upper = 1.0
_norm_function_kwargs = {"Masses": {"takelog": True, "plusone": True,"lower": lower, "upper": upper},
                    "GFM_Metallicity": {"takelog": True, "plusone": True, "lower": 0.25, "upper": upper},
                    "GFM_StellarFormationTime": {"takelog": True, "plusone": True, "lower": 0.25, "upper": upper},
                    
}

_IMG_ORDER = data._image_fields[particle_type][_dim]
_num_gal = data.get_attribute("mass").shape[0]

n_jobs = -1  # Use all cores
for i, field in tqdm(enumerate(_IMG_ORDER)):
    # Get the image of the specified particle type and field
    image = data.get_image(particle_type, field, dim = dim)
    norm_params = _norm_function_kwargs[field] if field in _norm_function_kwargs.keys() else {}
    # Normalize the images in parallel
    images_normalized = joblib.Parallel(n_jobs=n_jobs)(joblib.delayed(_norm_function)(img, **norm_params) for img in image)
    # Flatten the images
    images_flattened = np.array([img.flatten() for img in images_normalized])
    # Insert the flattened images into the correct place in the datamatrix
    
    datamatrix[:, i*np.prod(_IMG_SHAPE):(i+1)*np.prod(_IMG_SHAPE)] = images_flattened
                
            
print("Datamatrix created.")


# Create 1000 PCA fits with random subsamples of the data
from tqdm import trange
from sklearn.decomposition import PCA
data = datamatrix
pca = PCA(n_components = 215)
explained_variance_ratio=[]
for i in trange(400):
    random_samples = np.random.choice(data.shape[0], int(data.shape[0]*0.75), replace=False)
    subsample = data[random_samples]
    pca.fit(subsample)
    explained_variance_ratio.append(pca.explained_variance_ratio_)
evr = np.array(explained_variance_ratio)
np.save("EVR_boostrap_3d_215_eigengalaxies.npy",evr)