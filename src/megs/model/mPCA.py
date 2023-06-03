import h5py
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

from ..data import DataLoader

class mPCA():
    '''Class to perform a PCA on the data.
    
    This class is a wrapper around the sklearn PCA class. It allows to perform a PCA on the data loaded by the DataLoader class.
    
    Parameters:
    -----------
    data : DataLoader
        The data to perform the PCA on. The data must be an instance of DataLoader.
    particle_type : str, optional
        The particle type on which to apply the PCA. If None, the first particle type in the dataset is used.
    norm_function : function, optional
        The function to use to normalize the data before applying the PCA. If None, no normalization is applied.
    norm_function_kwargs : dict, optional
        The keyword arguments to pass to the normalization function.
    mask : np.ndarray, optional
        The mask to apply to the data before applying the PCA. If None, no mask is applied.
    dim : int, optional
        The dimension of the data. If 2, the data is assumed to be 2D. If 3, the data is assumed to be 3D.

    Examples:
    ---------
    >>> data = DataLoader("data.hdf5")
    >>> pca = mPCA(data, particle_type = "stars",dim=2)
    >>> pca.fit()
    >>> eigengalaxies = pca.get_eigengalaxies()
    '''
    
    def _identity(self, x, **kwargs):
        '''Identity function
        
        Used as default if no normalization function is specified.
        '''
        return x
    
    def __init__(self,data, particle_type = None, norm_function = None, norm_function_kwargs=None, mask = None, dim=2):
        # Check if data is instance of DataLoader
        if not isinstance(data, DataLoader):
            raise ValueError("data must be an instance of DataLoader")
        
        # Set the particle type on which to apply the PCA (TODO maybe add option to apply PCA on multiple particle types)
        if particle_type is None:
           # Check if data has multiple particle types
            if len(data._particles_keys) > 1:
                raise ValueError("PCA can only be applied to a single particle type. Please select a single particle type when initializing the PCA object.")
            else:
                # Set the particle type to the only particle type in the dataset
                self.particle_type = data._particles_keys[0]
        else:
            #Check if the particle type is valid
            if particle_type not in data._particles_keys:
                raise ValueError(f"Particle type {particle_type} not found. Valid particle types are: {data._particles_keys}")
            else:
                self.particle_type = particle_type
                
        # Check if mask is set
        if mask is None:
            self.mask = np.ones(data.get_attribute("mass").shape[0]).astype(bool)
        else: 
            self.mask = mask
        
        self.data = data
        self._dim = "dim2" if dim == 2 else "dim3"
        self._norm_function = norm_function if norm_function is not None else self._identity # Set the normalization function to the identity function if no normalization function is specified
        self._norm_function_kwargs = norm_function_kwargs
        self._IMG_ORDER = self.data._image_fields[self.particle_type][self._dim]
        self._IMG_SHAPE = self.data.get_image(self.particle_type, self.data._image_fields[self.particle_type][self._dim][0], index = 0, dim = dim).shape
        # Initialize the datamatrix to of shape (n_galaxies, 0)
        
        self.datamatrix = np.empty((data.get_attribute("mass").shape[0], 0)) 
        # Create datamatrix
        self._create_datamatrix(self._dim)
        
        #Get the image shape of one galaxy image
        

            
    def _create_datamatrix(self, dim):
        # Get the fields images of the specified particle type
        
        print("Creating datamatrix with the following fields:")
        print("===============================================")
        print("Particle type: ", self.particle_type)
        print("Fields: ", self._IMG_ORDER)
        print("Dimension: ", dim)
        #print("norm_function_kwargs: ", self._norm_function_kwargs)
        # Say that for the fields that are not specified in the norm_function_kwargs, the default arguments are used
        print("Default arguments are used for the fields that are not specified in the norm_function_kwargs")
        print("===============================================")
        
    
    
        for field in self._IMG_ORDER:
            # Get the image of the specified particle type and field
            image = self.data.get_image(self.particle_type, field, dim = dim)
            norm_params = self._norm_function_kwargs[field] if field in self._norm_function_kwargs.keys() else {}
            image = np.array([self._norm_function(img, **norm_params).flatten() for img in image])
            
            # Normalize the image
           
           # if self._norm_function is not None:
                # Check if the norm function has arguments for the specified field
           #     if field in self._norm_function_kwargs.keys():
           #         image = self._norm_function(image, **self._norm_function_kwargs[field])
           #     else:
                    # If no arguments are specified, use the default arguments
           #         image = self._norm_function(image)
            # Add the image to the datamatrix
            
            self.datamatrix = np.concatenate((self.datamatrix, image), axis=1)
         
        print("Created datamatrix with shape: ", self.datamatrix.shape)
        
        
        
    def fit(self, n_components=None, show_results = True,**kwargs):
        '''Fit the PCA to the datamatrix
        
        This function fits the PCA to the datamatrix. The scores, eigengalaxies and inverse_transformed_datamatrix are calculated and stored as attributes of the PCA object.
        
        Parameters:
        -----------
        n_components : int, optional
            Number of components to keep. If None, all components are kept.
        **kwargs : dict, optional
            Keyword arguments to pass to the sklearn PCA object.
        '''
        self.pca = PCA(n_components=n_components, **kwargs)
        self.scores = self.pca.fit_transform(self.datamatrix)
        self.eigengalaxies = self.pca.components_.reshape(self.pca.components_.shape[0],len(self.data._image_fields[self.particle_type][self._dim]), *self._IMG_SHAPE)
        self.inverse_transformed_datamatrix = self.pca.inverse_transform(self.scores)
        if show_results:
            self.show_results()
            
    def show_results(self):
        '''Show the results of the PCA'''
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        ax[0].plot(self.pca.explained_variance_ratio_)
        ax[0].set_xlabel("Component")
        ax[0].set_ylabel("Explained variance ratio")
        ax[1].plot(np.cumsum(self.pca.explained_variance_ratio_))
        ax[1].set_xlabel("Component")
        ax[1].set_ylabel("Cumulative explained variance ratio")
        plt.show()
        
        
        field_length = len(self.data._image_fields[self.particle_type][self._dim]) 
        
        #Loop over different image fields
        for index,field in enumerate(self.data._image_fields[self.particle_type][self._dim]):

            

            rows = int(np.ceil(np.sqrt(self.pca.n_components)))
            cols = int(np.ceil(self.pca.n_components/rows))
            fig, ax = plt.subplots(rows, cols, figsize=(cols*3, rows*3))
            for i in range(rows):
                for j in range(cols):
                    if i*cols+j < self.pca.n_components:
                        ax[i, j].imshow(self.eigengalaxies[i*cols+j][index])
                        ax[i, j].set_title(f"Component {i*cols+j}")
                        ax[i, j].axis("off")
            fig.suptitle(f"Eigengalaxies:{field}")
            plt.show()
            
            
        
    
        # Plot the mean galaxy
        
        fig, ax = plt.subplots(1,field_length, figsize = (field_length*3, 3))
        for index,field in enumerate(self.data._image_fields[self.particle_type][self._dim]):
            ax[index].imshow(self.pca.mean_.reshape(field_length, *self._IMG_SHAPE)[index])
            ax[index].set_title(f"{field}")
            ax[index].axis("off")
        fig.suptitle("Mean Galaxy")
        plt.show() 
                
                
        #Calculate residue of random galaxy
        self.compare()
            
    def compare(self, index = None):
        '''Compare a galaxy with its reconstruction
        
        This function compares a galaxy with its reconstruction. If no index is specified, a random galaxy is chosen.
        
        Parameters:
        -----------
        index : int, optional
            Index of the galaxy to compare. If None, a random galaxy is chosen.
        '''
        field_length = len(self.data._image_fields[self.particle_type][self._dim])
        #Calculate residue of random galaxy
        if index is None:
            randomind = np.random.randint(0, self.datamatrix.shape[0])
        else:
            randomind = index
        inverse_images = self.inverse_transformed_datamatrix[randomind].reshape(field_length, *self._IMG_SHAPE)
        for index, field in enumerate(self.data._image_fields[self.particle_type][self._dim]):
            fig, ax = plt.subplots(1, 3, figsize=(6, 3))
            
            original = self.datamatrix[randomind].reshape(field_length, *self._IMG_SHAPE)[index]
            
            ax[0].imshow(original)
            ax[0].set_title(f"Original galaxy")
            ax[0].axis("off")
            ax[1].imshow(inverse_images[index])
            ax[1].set_title(f"Reconstructed galaxy")
            ax[1].axis("off")
            
            # Calculate residue
            residue = np.abs(original - inverse_images[index])
            #Divide with original where original != 0
            residue[original != 0] = residue[original != 0] / original[original != 0]
            print(residue.max(), residue.min())
            ax[2].imshow(residue, vmin=0, vmax=1)
            ax[2].set_title(f"Residue")
            ax[2].axis("off")
            
            fig.suptitle(f"Comparison of {field} field")
            plt.show()
            
    def get_eigengalaxies(self):
        '''Return the eigengalaxies'''
        return self.eigengalaxies