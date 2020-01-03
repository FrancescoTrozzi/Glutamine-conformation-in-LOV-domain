import frequency_calculator

#indices must be numpy array of shape (any, 3) with rows being atom number of Donor, Hydrogen, Acceptor

class Hbonds(object):
    """
    Compute angle and distances of atoms involved in HBonding and then check the existance of H-bond according the BH criteria.

    Attributes
    ----------
        trajectory : numpy array
            frame of the MD trajectory
        indexDonor : integer
            0th index of the Donor atom in the Hydrogen Bond
        indexHydrogen : integer
            0th index of the Hydrogen atom in the Hydrogen Bond
        indexAcceptor : integer
            0th index of the Acceptor atom in the Hydrogen Bond
        indexPreAcceptor : integer
            0th index of the Acceptor atom in the Hydrogen Bond

    Methods
    -------
        array_generator : 
            creates numpy arrays of the atoms involved in the distances and angles used to define the hydrogen bond
        geometry_criteria :
            creates numpy arrays of the distances and angles to be checked for the existance of the hydrogen bond 
        criteria_checker :
            check if given angles and distances provided, the Hydrogen bond exists
    """

    def __init__(self, trajectory, indexDonor, indexHydrogen, indexAcceptor, indexPreAcceptor):
        """Initialization of all data necessary to compute excistence of hydrogen bonds"""
        self.trajectory = trajectory
        self.indexDonor = indexDonor
        self.indexHydrogen = indexHydrogen
        self.indexAcceptor = indexAcceptor
        self.indexPreAcceptor = indexPreAcceptor

    def array_generator(self):
        """Generate distances and angle arrays to check
        
        Parameters
        ----------
            MD trajectory and indices of atoms involved in the H-Bond (see class attributes)

        Returns
        -------
            DA_array : numpy array
                array composed of the indices of the Donor and Acceptor atoms
            HA_array : numpy array
                array composed of the indices of the Hydrogen and Acceptor atoms
            DHA_array : numpy array
                array composed of the indices of the Donor, Hydrogen and Acceptor atoms
            aAD_array : numpy array
                array composed of the indices of the pre-Acceptor, Acceptor and Donor atoms
            aAH_array : numpy array
                array composed of the indices of the pre-Acceptor, Acceptor and Hydrogen atoms
        """
        DA_array = np.array([[self.indexDonor, self.indexAcceptor]])
        HA_array = np.array([[self.indexHydrogen, self.indexAcceptor]])
        DHA_array = np.array([[self.indexDonor, self.indexHydrogen, self.indexAcceptor]])
        aAD_array = np.array([[self.indexPreAcceptor, self.indexAcceptor, self.indexDonor]])
        aAH_array = np.array([[self.indexPreAcceptor, self.indexAcceptor, self.indexHydrogen]])
        return DA_array, HA_array, DHA_array, aAD_array, aAH_array

    def geometry_criteria(self):
        """Generate arrays of distances and angles that define the H-bond
        
        Parameters
        ---------- 
            Numpy arrays composed of the 0th atomic indices of atoms involved in the distances and angles needed to compute the H-Bond  (see returns of array_generetor method)
        
        Returns
        ------- 
            DA_dist : numpy array
                array of the H-bond distance between the Donor and Acceptor in the given trajectory frame
            HA_dist : numpy array
                array of the H-bond distance between the Hydrogen and Acceptor in the given trajectory frame
            DHA_angle : numpy array
                array of the H-bond angle between the Donor, Hydrogen and Acceptor in the given trajectory frame
            aAD_angle : numpy array
                array of the H-bond angle between the pre-Acceptor, Acceptor and Donor in the given trajectory frame
            aAH_angle : numpy array 
                array of the H-bond angle between the pre-Acceptor, Acceptor and Hydrogen in the given trajectory frame 
        """
        DA_array, HA_array, DHA_array, aAD_array, aAH_array = self.array_generator()
        DA_dist = md.compute_distances(self.trajectory, DA_array)*10
        HA_dist = md.compute_distances(self.trajectory, HA_array)*10
        DHA_angle = radians_to_degrees(md.compute_angles(self.trajectory, DHA_array))
        aAD_angle = radians_to_degrees(md.compute_angles(self.trajectory, aAD_array))
        aAH_angle = radians_to_degrees(md.compute_angles(self.trajectory, aAH_array))
        return DA_dist, HA_dist, DHA_angle, aAD_angle, aAH_angle

    def criteria_checker(self):
        """Check if the hydrogen bond can exists based on the angles and distances arrays.
        
        Parameters
        ---------- 
            Numpy array containing the values of distances and angle needed to compute the H-bond (see returns from geometry_criteria method)

        Returns
        -------
            True is H-bond exists. 
        """
        if (self.indexDonor, self.indexHydrogen, self.indexAcceptor, self.indexPreAcceptor):
            DA_dist, HA_dist, DHA_angle, aAD_angle, aAH_angle = self.geometry_criteria()
        elif ("DA_dist", "HA_dist", "DHA_angle", "aAD_angle", "aAH_angle" in args):
            DA_dist = np.load(DA_dist)
            HA_dist = np.load(HA_dist)
            DHA_angle = np.load(DHA_angle)
            aAD_angle = np.load(aAD_angle)
            aAH_angle = np.load(aAH_angle)
        else:
            print("The arguments given are not enough to check the presence of hydrogen bond according to BH criteria")
        
        paramenters = np.concatenate((DA_dist, HA_dist, DHA_angle),  axis=1)
        paramenters = np.concatenate((paramenters, aAD_angle, aAH_angle), axis=1)

        for i in range(len(paramenters)):
            if paramenters[i,0] < 3.9 and paramenters[i,1] < 2.5 and paramenters[i,2] > 90 and paramenters[i,3] > 90 and paramenters[i,4] > 90:
                return True