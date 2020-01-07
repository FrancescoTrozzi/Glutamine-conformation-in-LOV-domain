import frequency_calculator as fc
import oop_hb_detector as hbd
import argparse
import sys
import mdtraj as md
import numpy as np

parser = argparse.ArgumentParser(description= 'This script check first the presence of specific HB between Gln and FMN and then check which functioning state Gln is based on dihedrals or angles')
parser.add_argument("-s", "--state", required = True, help="either dark or light (all lowercase)")
parser.add_argument("-t", "--trajectory", required = True, help="structure name of the dcd file")
parser.add_argument("-i", "--indices", required = True, help="numpy array of indices. In the order: Donor, Hydrogen, Acceptor and PreAcceptor in the Hydrogen Bond")
parser.add_argument("-gi", "--gln_indices", required = True, help="numpy array of Glutamine indices.")
args = parser.parse_args()


#Load data
traj1 = md.load(args.trajectory+'_01_50000_tot.dcd', top=args.trajectory+'.pdb')
traj2 = md.load(args.trajectory+'_02_50000_tot.dcd', top=args.trajectory+'.pdb')
traj3 = md.load(args.trajectory+'_03_50000_tot.dcd', top=args.trajectory+'.pdb')

#Sum the trajectories
traj = traj1 + traj2 + traj3

#Load numpy arrays for of indices for H-Bonds and Gln conformation
indices = np.load(args.indices)
gln_indices = np.load(args.gln_indices)


def main():
    if args.state == 'dark':
        dihedrals_radians_A = fc.gln_array(traj, gln_indices[0,0], gln_indices[0,1], gln_indices[0,2], gln_indices[0,3])
        dihedrals_degrees_A = fc.radians_to_degrees(dihedrals_radians_A)
        dihedrals_radians_B = fc.gln_array(traj, gln_indices[1,0], gln_indices[1,1], gln_indices[1,2], gln_indices[1,3])
        dihedrals_degrees_B = fc.radians_to_degrees(dihedrals_radians_B)

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
                dark_I, dark_II = fc.exp_vs_burI(dihedrals_degrees_A[i])
                dark_I_tot += dark_I
                dark_II_tot += dark_II

            if h_B.criteria_checker():
                bonding_frames_n += 1
                dark_I, dark_II = fc.exp_vs_burI(dihedrals_degrees_B[i])
                dark_I_tot += dark_I
                dark_II_tot += dark_II

        dark_I_freq = (dark_I_tot/bonding_frames_n)*100
        dark_II_freq = (dark_II_tot/bonding_frames_n)*100
        print(str(args.trajectory)+ ' Exposed: ' +str(dark_I_freq))
        print(str(args.trajectory)+ ' Semi-Buried: ' +str(dark_II_freq))

    elif args.state == 'light':
        angles_radians_A = fc.gln_array(traj, gln_indices[0,0], gln_indices[0,1], gln_indices[0,2])
        angles_degrees_A = fc.radians_to_degrees(angles_radians_A)
        angles_radians_B = fc.gln_array(traj, gln_indices[1,0], gln_indices[1,1], gln_indices[1,2])
        angles_degrees_B = fc.radians_to_degrees(angles_radians_B)

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
        print(str(args.trajectory)+ 'Buried I: ' +str(light_I_freq))
        print(str(args.trajectory)+ 'Buried II: ' +str(light_II_freq))

# Execute
main()
