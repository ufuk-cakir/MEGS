'''This module contains the class for specific simulations. Currently only IllustrisTNG is implemented.
You can add your own simulation by creating a new class based on the following template:

class YourSimulation():
    def __init__(self, halo_id, particle_type, data_path):
        self.data_path = data_path #Path to the data.
        self.halo_id= halo_id #Halo ID of a subhalo
        self.particle_type = particle_type #Particle type of the subhalo
        
        #Add any other attributes you need. 
        # These can later be accessed by the galaxy object using the getattr function, and can be saved in the HDF5 file.
        
        self._load_data()
    def _load_data(self):
        #Load the data from the snapshot and subhalo catalog.
        #The data needs to be stored in the following attributes:
        
        self.particle_coordinates #Coordinates of the particles,used for face-on rotation and image rendering.
        self.particle_masses #Masses of the particles, used for face-on rotation and image rendering.
        self.center #Center of the subhalo. Used for centering the galaxy
        self.halfmassrad #Half mass radius of the subhalo. Used for the image rendering
        self.hsml #Smoothing length of the subhalo. Used for the image rendering.
        
        
    def get_field(self, field):
        #Get the field from the snapshot.
        #The field should be returned as a numpy array and converted to physical units.
        #This field is then later used to create the weighted image.

        return field

Note that the class name should be the same as the simulation name, since this is used to create the galaxy object.
'''



import numpy as np
from astropy.cosmology import Planck15 as cosmo
import illustris_python as il 


def scale_to_physical_units(x, field):
    '''get rid of the Illustris units.'''

    if field == 'Masses':
        return x * 1e10 / (cosmo.H(0).value)

    elif field == 'Coordinates':
        return x / (cosmo.H(0).value)

    elif field == 'SubfindHsml':
        return x / (cosmo.H(0).value)

    elif field == 'SubfindDensity':
        return x * 1e10 * (cosmo.H(0).value) * (cosmo.H(0).value)

    elif field == 'GFM_StellarFormationTime':
        #Calculates Age of Stars
        return (cosmo.age(0).value-cosmo.age(1 / x - 1).value)*1e9 #Gyr
    elif field =="GFM_Metallicity":
        return(x/0.0127) #Solar Metallicity
    else:
        print("No unit conversion for Field {}. Return without changes.".format(field))
        return x


class IllustrisTNG():
    '''Class for the IllustrisTNG simulation.'''
    
    def __init__(self, halo_id, particle_type,base_path, snapshot):
        self.base_path = base_path
        self.halo_id= halo_id
        self.particle_type = particle_type
        self.snapshot = snapshot
       
        self._load_data() 
    
    def _load_data(self):
        '''Load the data from the snapshot and subhalo catalog.'''
        self.subhalo = il.groupcat.loadSingle(self.base_path, self.snapshot, subhaloID=self.halo_id)
        self.center = self.subhalo['SubhaloPos']
        self.mass = scale_to_physical_units(self.subhalo['SubhaloMassType'][il.util.partTypeNum(self.particle_type)], 'Masses')
        self.halfmassrad = self.subhalo['SubhaloHalfmassRadType'][il.util.partTypeNum(self.particle_type)]
        self.halfmassrad_DM= self.subhalo['SubhaloHalfmassRadType'][il.util.partTypeNum('DM')]
        self.particles = il.snapshot.loadSubhalo(self.base_path, self.snapshot, self.halo_id, self.particle_type)
    
        if self.particle_type == "stars": 
            # Get only real stars, not wind particles
            self.real_star_mask = np.where(self.particles["GFM_StellarFormationTime"]>0)[0]
            self.hsml = self.particles["StellarHsml"][self.real_star_mask]
        else:
            self.real_star_mask = np.ones(len(self.particles["Coordinates"]), dtype=bool)
            #Is this correct? Is this the smoothing length used for visualization?
            self.hsml = self.particles["SubfindHsml"]
        
        self.particle_coordinates = self.particles["Coordinates"][self.real_star_mask]
        self.particle_masses = scale_to_physical_units(self.particles["Masses"][self.real_star_mask], 'Masses')
    def get_field(self, field):
        '''Load a field from the particle data. Used for the image generation.
        The field is returned as a numpy array and converted to physical units.

        Parameters
        ----------
        field : str
            Name of the field to load. The field should be stored in the snapshot.  

        Returns
        -------
        numpy.array
            The field converted to physical units. 
        

        Examples
        --------
        >>> galaxy = Galaxy("IllustrisTNG", halo_id=0, particle_type="stars")
        >>> galaxy.get_field("GFM_StellarFormationTime")

        '''
        return scale_to_physical_units(self.particles[field][self.real_star_mask], field)
