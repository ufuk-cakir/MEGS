from sklearn.decomposition import PCA
import numpy as np





#------------- Coordinate Rotation -------------
def radial_distance(coords,center):
    d = coords-center
    r_i = d**2
    r= np.sqrt(np.sum(r_i, axis = 1))
    return(r)

def moment_of_intertia_tensor(coordinates, particle_masses, rHalf, subhalo_pos): 
    """ Calculate the moment of inertia tensor (3x3 matrix) for a subhalo-scope particle set."""

    rad = radial_distance(coords = coordinates, center = subhalo_pos) 
    wGas = np.where((rad <= 2.0*rHalf))[0] 
    masses = particle_masses[wGas] 
    xyz = coordinates[wGas,:]   

    xyz = np.squeeze(xyz)


    if xyz.ndim == 1:
        xyz = np.reshape( xyz, (1,3) )

    for i in range(3):
        xyz[:,i] -= subhalo_pos[i]

    I = np.zeros( (3,3), dtype='float32' )

    I[0,0] = np.sum( masses * (xyz[:,1]*xyz[:,1] + xyz[:,2]*xyz[:,2]) )
    I[1,1] = np.sum( masses * (xyz[:,0]*xyz[:,0] + xyz[:,2]*xyz[:,2]) )
    I[2,2] = np.sum( masses * (xyz[:,0]*xyz[:,0] + xyz[:,1]*xyz[:,1]) )
    I[0,1] = -1 * np.sum( masses * (xyz[:,0]*xyz[:,1]) )
    I[0,2] = -1 * np.sum( masses * (xyz[:,0]*xyz[:,2]) )
    I[1,2] = -1 * np.sum( masses * (xyz[:,1]*xyz[:,2]) )
    I[1,0] = I[0,1]
    I[2,0] = I[0,2]
    I[2,1] = I[1,2]

    return I



def rotation_matrix(inertiaTensor, return_value = "face-on"):
    """ Calculate 3x3 rotation matrix by a diagonalization of the moment of inertia tensor.
    Note the resultant rotation matrices are hard-coded for projection with axes=[0,1] e.g. along z. """

    # get eigen values and normalized right eigenvectors
    eigen_values, rotation_matrix = np.linalg.eig(inertiaTensor)

    # sort ascending the eigen values
    sort_inds = np.argsort(eigen_values)
    eigen_values = eigen_values[sort_inds]

    # permute the eigenvectors into this order, which is the rotation matrix which orients the
    # principal axes to the cartesian x,y,z axes, such that if axes=[0,1] we have face-on
    new_matrix = np.matrix( (rotation_matrix[:,sort_inds[0]],
                             rotation_matrix[:,sort_inds[1]],
                             rotation_matrix[:,sort_inds[2]]) )

    phi = np.random.uniform(0, 2*np.pi)
    theta = np.pi / 2
    psi = 0

    A_00 =  np.cos(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.sin(psi)
    A_01 =  np.cos(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.sin(psi)
    A_02 =  np.sin(psi)*np.sin(theta)
    A_10 = -np.sin(psi)*np.cos(phi) - np.cos(theta)*np.sin(phi)*np.cos(psi)
    A_11 = -np.sin(psi)*np.sin(phi) + np.cos(theta)*np.cos(phi)*np.cos(psi)
    A_12 =  np.cos(psi)*np.sin(theta)
    A_20 =  np.sin(theta)*np.sin(phi)
    A_21 = -np.sin(theta)*np.cos(phi)
    A_22 =  np.cos(theta)

    random_edgeon_matrix = np.matrix( ((A_00, A_01, A_02), (A_10, A_11, A_12), (A_20, A_21, A_22)) )

    # prepare return with a few other useful versions of this rotation matrix
    r = {}
    r['face-on'] = new_matrix
    r['edge-on'] = np.matrix( ((1,0,0),(0,0,1),(0,-1,0)) ) * r['face-on'] # disk along x-hat
    r['edge-on-smallest'] = np.matrix( ((0,1,0),(0,0,1),(1,0,0)) ) * r['face-on']
    r['edge-on-y'] = np.matrix( ((0,0,1),(1,0,0),(0,-1,0)) ) * r['face-on'] # disk along y-hat
    r['edge-on-random'] = random_edgeon_matrix * r['face-on']
    r['phi'] = phi
    
    return r[return_value]


def face_on_rotation(rHalf, subhalo_pos,coordinates, particle_masses):
    """ Rotate the subhalo-scope particle set to face-on orientation."""
    I = moment_of_intertia_tensor(coordinates=coordinates, rHalf=rHalf,particle_masses=particle_masses, subhalo_pos=subhalo_pos)
    rot_matrix = rotation_matrix(inertiaTensor=I, return_value = "face-on")
    
    
    #Rotate Particles to face-on with the calculated Rotation Matrix
    pos = coordinates- subhalo_pos
    rot_pos= np.dot(rot_matrix, pos.T).T
    rotated_particles = np.asarray(rot_pos)
    return rotated_particles    






#------------- Vertical Rotation -------------

def calc_rotation_matrix(angle):
    '''Calculate rotation matrix around z-axis for a given angle(rad).'''
    rot_mat = np.array([[np.cos(angle),-np.sin(angle),0],
                           [np.sin(angle), np.cos(angle),0],
                           [0 , 0, 1]])
    return(rot_mat)

def get_horizontal_angle(img):
    fit=PCA(n_components=2).fit(np.argwhere(img>=np.quantile(img,.75)))
    return np.arctan2(*fit.components_[0])


def horizontal_rotation(coordinates, halfmassrad,plotfactor=10):
    #First Create Dummy hist
    hist, xedges, yedges = np.histogram2d(coordinates[:,0], coordinates[:,1], bins=(128,128), range =((-halfmassrad*plotfactor,halfmassrad*plotfactor),(-halfmassrad*plotfactor,halfmassrad*plotfactor)))
    angle = get_horizontal_angle(hist)
    
    #Rotate
    horizontal_rotation_matrix = calc_rotation_matrix(-angle)
    rotated_coordinates = np.dot(horizontal_rotation_matrix, coordinates.T).T
    return rotated_coordinates

    