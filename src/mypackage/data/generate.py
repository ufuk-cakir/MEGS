'''This file contains the functions to generate the galaxy data and save it to the HDF5 file.
'''

import h5py
import os
import numpy as np
from tqdm import tqdm

from galaxy import Galaxy


def _create_data_structure(n_galaxies, image_res,galaxy_parameters, particle_types, fields, path="./"):
    '''Creates the HDF5 data structure for the galaxy data
    #TODO: Add the ability to save the images as a 3D array
    #TODO: Add the ability to set different fields for different particle types
    Parameters
    ----------
    n_galaxies: int
        Number of galaxies to be saved
    image_res: tuple
        Image resolution (x,y). #TODO: Change to be flexible (x,y,z) for 3D images
    galaxy_parameters: list
        List of galaxy parameters to be saved. These need to be attributes of the Galaxy class and are saved in the "Galaxies/Attributes" group.
    particle_types: list
        List of particle types to be saved e.g. ["stars", "gas"]. 
    fields: dict
        Dictionary of fields to be saved, where the key is the field name and the value is a boolean indicating if the field is mass weighted or not. 
        The fields are the same for all particle types and are saved in the "Images" group. The value is used to calculate the image later. #TODO maybe change this
    path: str
        Path to the HDF5 file
    
    
    Example:
    --------
    >>> create_data_structure(1000, (64,64), galaxy_paramters=["mass","halo_id"], particle_types =["stars"], {"Masses":False, "GFM_Metallicity":False, "GFM_StellarFormationTime":True})
    
    This will create a "galaxy_data.hdf5" HDF5 file with the following structure and save it to the current directory:
    ---------------------------------------------------------
    Galaxies
        Attributes
            mass: (1000,)
            halo_id: (1000,)
        Particles
            stars
                Images
                    Masses: (1000,64,64)
                    GFM_Metallicity: (1000,64,64)
                    GFM_StellarFormationTime: (1000,64,64)
    ---------------------------------------------------------
    '''
    # Open the HDF5 file in "w" mode to create a new file
    with h5py.File(os.path.join(path, "galaxy_data.hdf5"), "w") as f:
        # Create the Galaxies group
        galaxies_group = f.create_group("Galaxies")
        galaxy_attributes = galaxies_group.create_group("Attributes")
        # Create the datasets for the galaxy parameters
        for parameter in galaxy_parameters:
            galaxy_attributes.create_dataset(parameter, shape=(n_galaxies,), maxshape=(None,))    
        
        particles_group = galaxies_group.create_group("Particles")
        # Create the Particle Types group
        for particle_type in particle_types:
            particle_type_group = particles_group.create_group(particle_type)
            
            # Create the Images group
            images_group = particle_type_group.create_group("Images")
            
            # Create the datasets for the images
            for field in fields:
                images_group.create_dataset(field, shape=(n_galaxies, *image_res), maxshape=(None, None,None))
                



def _calculate_images(simulation,halo_ids,fields,plot_factor, image_res, path="./",norm_params = {},**kwargs):
    '''Calculates the images for the galaxies and saves them to the HDF5 file. Needs to be run after _create_data_structure() method.
    
    Parameters
    ----------
    simulation: str
        Simulation name (e.g. IllustrisTNG). Used to initialise the Galaxy class. The Galaxy class for a specific simulation should be defined in the simulations.py file. #TODO: Change name
    halo_ids: list
        List of halo IDs to calculate the images for. The halo_ids are used to load the galaxy data from the simulation.
    fields: dict
        Dictionary of fields to be saved, where the key is the field name. The values of the dictionary are passed to the get_image() method of the Galaxy class.
        Example: fields = {"Masses":{"mass_weighted":False,
                                     "normed":True}, 
                  "GFM_Metallicity":{"mass_weighted":True,
                                     "normed":True}, 
                  "GFM_StellarFormationTime":{"mass_weighted":True,
                                              "normed":False}
                }
    plot_factor: float
        Factor for the image range. The image range is calculated as halfmass_radius*plot_factor and the image is centred on the galaxy centre.
        For the halfmass_radius only the particle type specified in the particle_type argument are used.
    path: str
        Path to the HDF5 file to save the data to. This file should be created using create_data_structure() method.
    **kwargs: dict
        Keyword arguments passed to the Galaxy class. Halo ID and particle type are overwritten in the loop.
        e.g. {"base_path":basePath,"halo_id":0,"particle_type": "stars", "snapshot":99} for IllustrisTNG

    '''
    # Check if the HDF5 file exists, which should be created using create_data_structure() method.
    if not os.path.exists(os.path.join(path, "galaxy_data.hdf5")):
        raise FileNotFoundError(f"{os.path.join(path, 'galaxy_data.hdf5')} does not exist. This should have been created using the _create_data_structure() method.")
    
    n_galaxies = len(halo_ids)
    # Open the HDF5 file in "append" mode
    with h5py.File(os.path.join(path,"galaxy_data.hdf5"), "a") as f:
        
        # Check if the "index_position" attribute exists
        if "index_position" in f.attrs:
            index_position = f.attrs["index_position"]
            
            # Ask the user if they want to continue from the last index position
            if input(f"Continue from index position {index_position}? (y/n): ") == "y":
                # Check if the index position is valid
                if index_position > n_galaxies:
                    raise ValueError(f"Index position {index_position} is greater than the number of galaxies {n_galaxies}")        
            else:
                # Reset the index position since the user does not want to continue from the last index position 
                index_position = 0
        else:
            #Attribute does not exist so set the index position to 0 to start the loop from the beginning
            index_position = 0
            
            
        f.attrs["index_position"] = index_position
        # Loop through the galaxies
        for index, haloid in tqdm(enumerate(halo_ids[index_position:])):
            # Create the galaxy object
            kwargs["halo_id"] = haloid
            
            g = Galaxy(simulation= simulation, **kwargs) 
            
            # Get the galaxy parameters
            for parameter in f["Galaxies/Attributes"].keys():
                if hasattr(g,parameter):
                    f["Galaxies/Attributes"][parameter][index] = getattr(g,parameter)
                else: 
                    raise ValueError(f"Galaxy class does not have the attribute {parameter}")
                
            # Get the particle data
            for particle_type in f["Galaxies"]["Particles"].keys():
                # Get the particle data
                for field in f["Galaxies"]["Particles"][particle_type]["Images"].keys():
                    # Get the image
                    #TODO: Currently img_res is a list of 2 elements. Need to change this to a single value
                    image = g.get_image(field=field, plotfactor= plot_factor, res =image_res[0],**fields[field])
                    f["Galaxies"]["Particles"][particle_type]["Images"][field][index] = image
            
            # Update the index position
            f.attrs["index_position"] += 1
            
        # Show the user that the images have been calculated
        print("Images calculated and saved to HDF5 file: ", os.path.join(path,"galaxy_data.hdf5"))
        
        
      
def generate_data(simulation, halo_ids, fields, plot_factor, image_res, galaxy_parameters, particle_types, path="./", **kwargs):
    '''Generates the data for the galaxies and saves it to an HDF5 file.
    
    This method creates the HDF5 file data structure and saves the galaxy parameters and images to the file. 
    
    Parameters
    ----------
    simulation: str
        Simulation name (e.g. IllustrisTNG). Used to initialise the Galaxy class. The Galaxy class for a specific simulation should be defined in the simulations.py file. #TODO: Change name
    halo_ids: list
        List of halo IDs to calculate the images for. The halo_ids are used to load the galaxy data from the simulation.
    fields: dict
        Dictionary of fields to be saved, where the key is the field name. The values of the dictionary are passed to the get_image() method of the Galaxy class.
        Example: fields = {"Masses":{"mass_weighted":False,
                                        "normed":True}, 
                            "GFM_Metallicity":{"mass_weighted":True,
                                                "normed":True},
                            "GFM_StellarFormationTime":{"mass_weighted":True,
                                                        "normed":False}
                            }
       For more information on the fields see the get_image() method of the Galaxy class.                         
    plot_factor: float
        Factor for the image range. The image range is calculated as halfmass_radius*plot_factor and the image is centred on the galaxy centre.
        For the halfmass_radius only the particle type specified in the particle_type argument are used.
    path: str
        Path to the HDF5 file to save the data to. This file should be created using create_data_structure() method.
    **kwargs: dict
        Keyword arguments passed to the Galaxy class. Halo ID and particle type are overwritten in the loop.
        e.g. {"base_path":basePath,"halo_id":0,"particle_type": "stars", "snapshot":99} for IllustrisTNG

    Example
    -------
    # Set up the parameters to call the generate_data() method
    >>> simulation = "IllustrisTNG" # Simulation name, used to initialise the Galaxy class
    >>> halo_ids = [0,1,2,3,4,5,6,7,8,9] # List of halo IDs to calculate the images for
    >>> particle_types = ["stars", "gas"]
    >>> img_res = [64,64]
    >>> plot_factor = 10
    >>> path = "./" # Path where the HDF5 file will be saved
    
    # Define the fields to be saved
    >>> fields = {"Masses":{"mass_weighted":False,
                            "normed":True},
                "GFM_Metallicity":{"mass_weighted":True,
                                    "normed":True},
                "GFM_StellarFormationTime":{"mass_weighted":True,
                                            "normed":False}
                }
                
    >>> galaxy_parameters = ["mass", "halo_id"] # List of galaxy parameters to be saved
    >>> generate_data(simulation, halo_ids, fields, plot_factor, img_res, galaxy_parameters, particle_types, path)
    
    
    
    '''
    n_galaxies = len(halo_ids)
    # Create the HDF5 file and the data structure
    _create_data_structure(n_galaxies=n_galaxies, image_res=image_res, galaxy_parameters=galaxy_parameters, particle_types=particle_types, fields=fields, path=path)
    
    # Calculate the images
    _calculate_images(simulation, halo_ids, fields, plot_factor, image_res, path, **kwargs)
        
        
    

    
    
    
# Function to load the data of a specific galaxy
def load_single_galaxy_data(path, galaxy_index):
    '''Loads the data of a specific galaxy from the HDF5 file created using generate_data() method.
    
    Parameters
    ----------
    path: str
        Path to the HDF5 file
    galaxy_index: int
        Index of the galaxy to load the data for. The index is the position of the galaxy in the HDF5 file.
    
    Returns
    -------
    galaxy_data: dict
        Dictionary of the galaxy data. The keys are the galaxy parameters and particle types.
        e.g. {"mass": 1e12, "halo_id": 1234, "stars": {"Masses":image, "GFM_Metallicity":image, "GFM_StellarFormationTime":image}}
    '''
    #Check if the HDF5 file exists
    if not os.path.exists(os.path.join(path,"galaxy_data.hdf5")):
        raise FileNotFoundError(f"{path} does not exist.")
    
    
    # Open the HDF5 file in "r" mode
    with h5py.File(os.path.join(path,"galaxy_data.hdf5"), "r") as f:
        # Get the galaxy parameters
        galaxy_data = {parameter:f["Galaxies/Attributes"][parameter][galaxy_index] for parameter in f["Galaxies/Attributes"].keys()}
        
        # Get the particle data
        for particle_type in f["Galaxies"]["Particles"].keys():
            # Get the particle data
            galaxy_data[particle_type] = {field:f["Galaxies"]["Particles"][particle_type]["Images"][field][galaxy_index] for field in f["Galaxies"]["Particles"][particle_type]["Images"].keys()}
            
    return galaxy_data