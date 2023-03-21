

from simulations import * #Import all Galaxy Classes of the simulations 

import numpy as np
import sys

from rotation import face_on_rotation, horizontal_rotation
from image_modules import image2D, norm



def str_to_class(classname):
    '''Converts a string to a class.'''
    if hasattr(sys.modules[__name__], classname):
        return getattr(sys.modules[__name__], classname)
    else:
        raise ValueError("Simulation {} not supported.".format(classname))
        


class Galaxy():
    '''Base class for all galaxies.
    
    Need following attributes:
     - get_field(): returns a field from the snapshot. Used to create the images.
     - coordinates: Coordinates of the particles in the galaxy.
    
    '''
    def __init__(self,simulation="IllustrisTNG", **kwargs):
        self.simulation = simulation
        # Gemeral Galaxy Properties??
        self.rotated_flag = False
        
        #Load the defined Galaxy Class
        self.galaxy_object = str_to_class(self.simulation)(**kwargs)
        self.coordinates = self._rotate_galaxy()
        
        
        #Set Atributes for the image rendering
        self.plot_factor =10
        self.res = 64
        if hasattr(self.galaxy_object, "hsml"):
            self.smoothing_length = self.galaxy_object.hsml
        else:
            #Calculate smoothing length maybe later (TODO)
            raise AttributeError("Galaxy object does not have a smoothing length.")
        
        
         
    
    def _rotate_galaxy(self, _plotfactor=10):
        '''Rotate the galaxy to face-on and horizontal orientation.'''
        if self.rotated_flag:
            return self.coordinates
        
        face_on_rotated_coords = face_on_rotation(coordinates=self.particle_coordinates,particle_masses=self.particle_masses, 
                                                        rHalf=self.halfmassrad, subhalo_pos=self.center)
        #maybe horizontal rotataion is not working properly
        horizontal_rotated_coords = horizontal_rotation(coordinates=face_on_rotated_coords, halfmassrad=self.halfmassrad, plotfactor=_plotfactor)        
        self.rotated_flag = True
        return horizontal_rotated_coords

    def __getattr__(self, name):
        #Delegate all other attributes to the galaxy object
        return getattr(self.galaxy_object, name)


    #-----------------Image Rendering-----------------#

    def _render_image_2D(self, field):
        '''Image Render Module for 2D images. Can be changed. Uses the image2D function from image_modules.py as default.'''
        img = image2D(coordinates=self.coordinates, R_half=self.halfmassrad, weights=field,
                        smoothing_length=self.smoothing_length, plot_factor=self.plot_factor, res=self.res)
        return img
    
                      
    def get_image(self, field, mass_weighted = True,normed= False,**kwargs):
        '''Get the image of a given field.'''
        if mass_weighted:
            #First create the mass image
            masses = self.get_field("Masses")
            mass_img = self._render_image_2D(masses)
            if field == "Masses":
                #If the field is mass, return the mass image
                image = mass_img
            else:
                #Create the mass weighted field image.
                weights = self.get_field(field)*masses #mass weighted weights
                weights_img= self._render_image_2D(weights) 
                image = weights_img/mass_img
                
        else:
            #Not mass weighted. Return the field image
            weights = self.get_field(field)
            image = self._render_image_2D(weights)
        
        if normed:
            image = norm(image, **kwargs)

        return image
    

    





