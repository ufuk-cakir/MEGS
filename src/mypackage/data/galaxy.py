
''' 
Galaxy Class
------------
This module contains the Galaxy class as the main class of the package. It is used to load a galaxy from a simulation and to render images of the galaxy.
One can specify the simulation the galaxy is from and the halo id of the galaxy. The Galaxy class will then load the galaxy data from the simulation based 
on the Class defined in the simulatoins module and rotate the galaxy to face-on and horizontal orientation.
The get_image method can then be used to render images of different fields of the galaxy, and chose if they should be mass weighted or not.

Example
-------
>>> galaxy = Galaxy("IllustrisTNG", halo_id=0, particle_type="stars") # Load a galaxy from the IllustrisTNG simulation. Note that the Class IllustrisTNG needs to be defined in the simulations.py file
>>> galaxy.get_image("GFM_StellarFormationTime", mass_weighted=True) # Render an image of the age of the stars in the galaxy
>>> galaxy.mass # Get the mass of the galaxy
'''

from simulations import * #Import all Galaxy Classes of the simulations 

import numpy as np
import sys


from image_modules import image2D, norm, face_on_rotation, horizontal_rotation


def str_to_class(classname):
    '''Converts a string to a class.'''
    #Check if the class is defined in the simulations.py file
    if hasattr(sys.modules[__name__], classname):
        #Return the class
        return getattr(sys.modules[__name__], classname)
    else:
        raise ValueError("Simulation {} not supported.".format(classname))
        


class Galaxy():
    '''Base class for all galaxies.
    
    Parameters
    ----------
    simulation : str, optional
        The simulation the galaxy is from. Curenntly only "IllustrisTNG" is supported. You can add your own simulation by adding a new class to the simulations.py file. 
    
    **kwargs : dict
        Additional arguments for the galaxy class. This is used to initialize the galaxy class of a specific simulation. See the documentation of the galaxy class defined in simulations.py for more information.
    
    '''
    def __init__(self,simulation="IllustrisTNG", **kwargs):
        self.simulation = simulation
        # Gemeral Galaxy Properties??
        self.rotated_flag = False
        
        #Load the defined Galaxy Class
        self.galaxy_object = str_to_class(self.simulation)(**kwargs)
        self.coordinates = self._rotate_galaxy()
        
        
        #Set default Atributes for the image rendering
        self.plot_factor =10
        self.res = 64
        if hasattr(self.galaxy_object, "hsml"):
            self.smoothing_length = self.galaxy_object.hsml
        else:
            #Calculate smoothing length maybe later (TODO)
            raise AttributeError("Galaxy object does not have a smoothing length.")
        
        
         
    
    def _rotate_galaxy(self, _plotfactor=10):
        '''Rotate the galaxy to face-on and horizontal orientation.
        
        This function is called when the galaxy object is initialized. It is not necessary to call it again. First the galaxy is rotated face-on and then horizontal. 
        The rotation is done using the rotation.py module.
        
        Parameters
        ----------
        _plotfactor : int, optional
            Factor used in the horizontal_rotation method defined in the rotation.py module to scale the image. The default is 10. For more information see the documentation of the horizontal_rotation method.
        
        Returns
        -------
        numpy.array
            The rotated coordinates of the galaxy.
        '''
        if self.rotated_flag:
            return self.coordinates
        
        face_on_rotated_coords = face_on_rotation(coordinates=self.particle_coordinates,particle_masses=self.particle_masses, 
                                                        rHalf=self.halfmassrad, subhalo_pos=self.center)
        #maybe horizontal rotataion is not working properly
        horizontal_rotated_coords = horizontal_rotation(coordinates=face_on_rotated_coords, halfmassrad=self.halfmassrad, plotfactor=_plotfactor)        
        self.rotated_flag = True
        return horizontal_rotated_coords

    def __getattr__(self, name):
        '''Delegate all other attributes to the galaxy object. This is used to access the attributes of the simulation galaxy class defined in the simulations.py file, 
        if they are not defined in the general Galaxy class.
        '''
        return getattr(self.galaxy_object, name)


    #-----------------Image Rendering-----------------#

    def _render_image_2D(self, field):
        '''Image Render Module for 2D images.
        This function is called by the get_image function. It renders the image using the image2D function from the image_modules.py file.
        You can change the image rendering by implementing your own image rendering function here.
        
        For more information of the default render method see the documentation of the image2D function.
        '''
        img = image2D(coordinates=self.coordinates, R_half=self.halfmassrad, weights=field,
                        smoothing_length=self.smoothing_length, plot_factor=self.plot_factor, res=self.res)
        return img
    
                      
    def get_image(self, field, mass_weighted = True,normed= False, res = None, plotfactor = None,**kwargs):
        '''Get the image of a given field.

        This function renders the image of a given field, which can be any field that is available in the snapshot and can be accessed by the get_field function of the galaxy object.
        The images can be mass weighted and normalized.
        
        It uses the _render_image_2D function to render the image. You can change the image rendering by implementing your own image rendering function there.
        
        Parameters
        ----------
        field : str
            The field to be rendered. Can be any field that is available in the snapshot. Used to call the get_field function of the galaxy object.
        mass_weighted : bool, optional
            If True, the image is mass weighted. The default is True. 
        normed : bool, optional
            If True, the image is normalized. The default is False.
        res : int, optional
            The resolution of the image. The default is None. If None, the resolution is set to the default value of the Galaxy class, which is 64.
        plotfactor : int, optional
            The plotfactor used to scale the image. The default is None. If None, the plotfactor is set to the default value of the Galaxy class, which is 10.
        **kwargs : dict
            Additional arguments for the normalization function. This is only used if normed is True. For more information see the documentation of the norm function defined in the image_modules.py file.
        
        Returns
        -------
        numpy.array
            The image of the given field.
        
        Examples
        --------
        >>> galaxy = Galaxy(halo_id=0, particle_type="stars", base_path="data", snapshot=99)
        >>> image = galaxy.get_image("Masses", mass_weighted=False, normed=True)
        >>> plt.imshow(image)
        '''
        #Set the resolution and plotfactor
        if res is not None:
            self.res = res
        if plotfactor is not None:
            self.plot_factor = plotfactor
        
        
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
    

    





