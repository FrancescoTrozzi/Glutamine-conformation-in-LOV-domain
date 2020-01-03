import frequency_calculator as fc
import oop_hb_detector as hbd
import argparse
import sys

parser = argparse.ArgumentParser(description= 'This script check first the presence of specific HB between Gln and FMN and then check which functioning state Gln is based on dihedrals or angles')
parser.add_argument("-s", "--state", required = True, help="either dark or light (all lowercase)")
parser.add_argument("-t", "--trajectory", required = True, help="structure name of the dcd file")
parser.add_argument("-i", "--indices", required = True, help="numpy array of indices. In the order: Donor, Hydrogen, Acceptor and PreAcceptor in the Hydrogen Bond")
parser.add_argument("-gi", "--gln_indices", required = True, help="numpy array of Glutamine indices.")
args = parser.parse_args()


#Load data
traj = md.load(args.trajectory+'.dcd', top=args.trajectory+'.pdb')
indices = np.load(args.indices)
gln_indices = np.load(args.gln_indices)


def main():
    if args.state == 'dark': 
        dihedrals_radians_A = fc.gln_array(gln_indices[0,0], gln_indices[0,1], gln_indices[0,2])
        dihedrals_degrees_A = fc.radians_to_degrees(angles_radians)
        dihedrals_radians_B = fc.gln_array(gln_indices[1,0], gln_indices[1,1], gln_indices[1,2])
        dihedrals_degrees_B = fc.radians_to_degrees(angles_radians)

        # for each frame of the trajectory check the existance of the H-bond
        for i in range(len(traj)):
            dark_I_tot = 0
            dark_II_tot = 0
            coor = traj[i]
            bonding_frames_n = 0
            active_frames_n = 0
            h_A = hbd.Hbonds(coor, indexDonor=indices[0,0], indexHydrogen=indices[0,1],
                           indexAcceptor=indices[0,2], indexPreAcceptor=indices[0,3])
            h_B = hbd.Hbonds(coor, indexDonor=indices[1,0], indexHydrogen=indices[1,1],
                           indexAcceptor=indices[1,2], indexPreAcceptor=indices[1,3])

            # if H-Bond is present calculate the type of GLN conformation
            if h_A.criteria_checker():
                bonding_frames_n += 1
                dark_I, dark_II = fc.Exp_vs_burI(angles_degrees_A[i])
                dark_I_tot += dark_I
                dark_II_tot += dark_II
            
	    if h_B.criteria_checker():
                bonding_frames_n += 1
                dark_I, dark_II = fc.Exp_vs_burI(angles_degrees_B[i])
                dark_I_tot += dark_I
                dark_II_tot += dark_II

        dark_I_freq = (dark_I_tot/bonding_frames_n)*100
        dark_II_freq = (dark_II_tot/bonding_frames_n)*100
        print(dark_I_freq)
        print(dark_II_freq)

    elif args.state == 'light':
        angles_radians_A = fc.gln_array(gln_indices[0,0], gln_indices[0,1], gln_indices[0,2])
        angles_degrees_A = fc.radians_to_degrees(angles_radians)
        angles_radians_B = fc.gln_array(gln_indices[1,0], gln_indices[1,1], gln_indices[1,2])
        angles_degrees_B = fc.radians_to_degrees(angles_radians)

        # for each frame of the trajectory check the existance of the H-bond
        for i in range(len(traj)):
            light_I_tot = 0
            light_II_tot = 0
            coor = traj[i]
            bonding_frames_n = 0
            active_frames_n = 0
            h_A = hbd.Hbonds(coor, indexDonor=indices[0,0], indexHydrogen=indices[0,1],
                           indexAcceptor=indices[0,2], indexPreAcceptor=indices[0,3])
            h_B = hbd.Hbonds(coor, indexDonor=indices[1,0], indexHydrogen=indices[1,1],
                           indexAcceptor=indices[1,2], indexPreAcceptor=indices[1,3])

             # if H-Bond is present calculate the type of GLN conformation
            if h_A.criteria_checker():
                bonding_frames_n += 1
                light_I, light_II = fc.burII_vs_burIII(angles_degrees[i])
                light_I_tot += light_I
                light_II_tot += light_II

            if h_B.criteria_checker():
                bonding_frames_n += 1
                light_I, light_II = fc.burII_vs_burIII(angles_degrees[i])
                light_I_tot += light_I
                light_II_tot += light_II

        # print the frequencies in which the two light conformations are present
        light_I_freq = (light_I_tot/bonding_frames_n)*100
        light_II_freq = (light_II_tot/bonding_frames_n)*100
        print(light_I_freq)
        print(light_II_freq)


# Execute
main()
