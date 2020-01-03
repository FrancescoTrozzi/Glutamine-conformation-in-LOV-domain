import numpy as np
import mdtraj as md 

def radians_to_degrees(radians):
    return radians * (180/np.pi)


def gln_array(*args):
    if len(args) == 4:
        dihedral = np.array([[args[0], args[1], args[2], args[3]]])
        gln_dihedral_array = md.compute_dihedrals(traj, dihedral)
        return gln_dihedral_array

    elif len(args) == 3:
        angle = np.array([[args[0], args[1], args[2]]])
        gln_angle_array = md.compute_angles(traj, angle)
        return gln_angle_array
    
    else:
        print('The amount of arguments (number of atoms) is not compatible with the selected light state')


def exp_vs_burI(dihedral):
        for k in dihedral:
            if k > 0:
                return 1, 0
            else:
                return 0, 1


def burII_vs_burIII(angle):
    if angle < 25:
        return 1, 0
    else:
        return 0, 1
