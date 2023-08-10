
import sys
sys.path.insert(0, '../src')

from megs.model.mPCA import mPCA
from megs.data import image, DataLoader, Galaxy



import numpy as np
import joblib


from tqdm import trange,tqdm
import pickle

scores = np.load('scores.npy')
eigen3d = np.load("eigengalaxies.npy")
means = np.load("means.npy").reshape(-1)
eigengalaxies = eigen3d.reshape(512,-1)




particle_type = "stars"
_dim = "dim2"
dim = 2
_IMG_SHAPE = (64,64,64) if dim == 3 else (64,64)
print("Loading data...")

PICKLE_DATA_PATH = "2dmorphmodel400.pkl"

if PICKLE_DATA_PATH is not None:
    with open(PICKLE_DATA_PATH, "rb") as f:
        model = pickle.load(f)
        datamatrix = model.datamatrix
        scores = model.get_scores()
        means = model.get_means().reshape(-1)
        eigengalaxies = model.get_eigengalaxies().reshape(scores.shape[1], -1)
        _num_gal = scores.shape[0]

else:

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


def reconstruct(index,ncomp = 215):
    score = scores[index][:ncomp]
    eig = eigengalaxies[:ncomp]
    reconstructed = np.dot(score, eig) + means
    return reconstructed

def sse(reconstructed, original):
    diff = np.sum((original - reconstructed) ** 2)
    sse = diff/np.sum(original)
    return sse



import os 
os.makedirs(f"rec_err_{dim}d",exist_ok=True)
for ncomp in range(1,215,10):
    rec_err = []
    for i in trange(_num_gal):
        rec = reconstruct(i,ncomp)
        orig = datamatrix[i]
        rec_err.append(sse(rec,orig))
    np.save(f"rec_err_{dim}d/rec_err_{ncomp}.npy",rec_err)