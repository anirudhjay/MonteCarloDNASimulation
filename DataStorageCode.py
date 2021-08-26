""" This file contains codes to store all data attained from the simulations """

import HelperFunction1 as Hf1


def store_positions(polymers, filename, timestep):

    file = open(filename + '.xyz', 'a')
    if not hasattr(polymers, "polymers"):
        raise TypeError(f"input given is type {type(polymers)} but required type is PolymerClass")
    file.write('timestep = ' + str(timestep))
    file.write('\n')
    pol_num = 0
    for pol, beads in polymers.polymers.items():
        for bead in beads:
            identify, pos = str(bead.id), bead.pos
            file.write(str(pol) + '    ')
            for cord in pos:
                if cord >= 0:
                    file.write(f"{cord:.5f}    ")
                else:
                    file.write(f"{cord:.5f}    ")
            file.write('\n')
        pol_num += 1
    file.write('\n')
    file.close()
    return


def store_bead_char(polymers, filename, timestep):

    file = open(filename + '.txt', 'a')
    if not hasattr(polymers, "polymers"):
        raise TypeError(f"input given is type {type(polymers)} but required type is PolymerClass")
    file.write('timestep = ' + str(timestep) + '\n')
    for pol, beads in polymers.polymers.items():
        for bead in beads:
            file.write(f'{pol},{bead.id},{bead.bead_class},{bead.size}\n')
    file.write('\n')
    file.close()
    return


def store_bonds(polymers, filename, timestep):

    file = open(filename + '.txt', 'a')
    bonds = ['BackBone', 'Exclusion', 'BeadBead', 'Contain', 'BeadFix' ]
    if not hasattr(polymers, "polymers"):
        raise TypeError(f"input given is type {type(polymers)} but required type is PolymerClass")
    file.write('timestep = ' + str(timestep) + '\n')
    all_bonds = [polymers.spring_potentials, polymers.non_overlap_potentials, polymers.bead_bead_potentials1,
                      polymers.sphere_limit_potentials, polymers.bead_positioning_potentials]
    for bond_sets in all_bonds:
        file.write(bonds[all_bonds.index(bond_sets)] + '\n')
        for bond in bond_sets:
            file.write(str(bond) + '\n')
    file.write('\n')
    return


def store_pol_char(polymers, filename, timestep):
    file = open(filename + '.txt', 'a')
    if not hasattr(polymers, "polymers"):
        raise TypeError(f"input given is type {type(polymers)} but required type is PolymerClass")
    file.write('timestep = ' + str(timestep) + '\n')
    for pol, beads in polymers.polymers.items():
        file.write(f'{pol},{Hf1.cal_radius_of_gyration(beads)}\n')
    file.write('\n')
    file.close()
    return




