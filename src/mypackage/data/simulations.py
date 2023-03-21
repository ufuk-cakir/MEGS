
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
            self.hsml = self.particles["StellarHsml"]
        else:
            self.real_star_mask = np.ones(len(self.particles["Coordinates"]), dtype=bool)
            #Is this correct? Is this the smoothing length used for visualization?
            self.hsml = self.particles["SubfindHsml"]
        
        self.particle_coordinates = self.particles["Coordinates"][self.real_star_mask]
        self.particle_masses = scale_to_physical_units(self.particles["Masses"][self.real_star_mask], 'Masses')
    def get_field(self, field):
        '''Load a field from the particle data. Used for the image generation.'''
        return scale_to_physical_units(self.particles[field][self.real_star_mask], field)
