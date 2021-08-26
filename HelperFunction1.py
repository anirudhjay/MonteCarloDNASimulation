import numpy as np
import random as rdm
import math
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import OrderedDict


def distance(x, y):

    """
    Code to find the 3D distance between any two points
    :param x: a list with the 3D coordinates of a point x
    :param y: a list with the 3D coordinates of a point y
    :return: The distance 'r' between x and y
    """
    r = ((x[0] - y[0])**2 + (x[1] - y[1])**2 + (x[2] - y[2])**2)**0.5
    return r


def within_limit(new_position, limit):

    dist = distance([0, 0, 0], new_position)
    if dist >= limit:
        return False
    return True


def plt_sphere(list_center, ax, list_radius=tuple([0.5]), alpha=0.5, color='blue'):

    for c, r in zip(list_center, list_radius):

        # draw sphere
        u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:50j]
        x = r*np.cos(u)*np.sin(v)
        y = r*np.sin(u)*np.sin(v)
        z = r*np.cos(v)

        ax.plot_surface(x-c[0], y-c[1], z-c[2], color=color, alpha=alpha)


def overlap(polymers, new_monomer, monomer_length):
    """
    To check if the position for the 'new_monomer', in the 'polymer', is greater than 'monomer_length', so that it gets
    accepted.
    :param polymers: A list of all the bead position present in the polymer (in order)
    :param new_monomer: A new monomer position ascertained from the function 'select_position'
    :param monomer_length: the length of the monomer to be added (this is also the distance, lesser than which,
    the beads are rejected
    :return: Whether the new monomer is overlapping with any other bead in the polymer
    """
    # loop to check if the new_monomer is within 1*monomer_length distance with any other bead in the polymer
    if not isinstance(polymers, OrderedDict):
        raise TypeError(
            f'class of the polymer input should be "Polymers", but type of polymers given is {type(polymers)}'
        )
    for pol, beads in polymers.items():
        for bead in beads:
            # distance function to be used
            dist = distance(bead(), new_monomer)
            if dist < monomer_length:
                return True
    return False


def select_position(monomer, monomer_length=1):
    """
    A Class Polymer helper function that identifies a random position for the next monomer, given the previous
    monomer in the polymer.
    :param monomer: This is the position of the last bead added in the 'Polymer' class(this is a helper function)
    :param monomer_length: The length of the random monomer two be added to the growing Polymer
    :return: The new position of a random monomer fitted with the previous monomer of size 'monomer_length'
    """

    ml = monomer_length
    theta, phi = rdm.uniform(0, 2*math.pi), rdm.uniform(0, 2*math.pi)
    next_position = [ml * math.cos(theta) * math.cos(phi), ml * math.sin(theta) * math.cos(phi),
                     ml * math.sin(phi)]
    return [monomer[i] + next_position[i]for i in range(3)]


def capture_bonds(bead, bonds):

    all_bonds_of_bead = [bond for bond in bonds if bond.id[0] == bead or bond.id[1] == bead]
    return all_bonds_of_bead


def check_length(old_energy, new_energy):

    if old_energy >= new_energy:
        return True
    else:
        r = rdm.random()
        p = math.e**(old_energy - new_energy)
        if r <= p:
            return True
        else:
            return False


def capture_contacts(lj_rep_bonds, polymer_size, resolution=1.5):

    contact_matrix = np.full((polymer_size, polymer_size), 0.0)
    for bond in lj_rep_bonds:
        length = bond.current_length
        if length <= resolution:
            beads = bond.id
            contact_matrix[int(beads[0]) - 1][int(beads[1]) - 1], contact_matrix[int(beads[1]) - 1][int(beads[0]) - 1]\
                = 1, 1
    return contact_matrix


def show_contacts(matrix_file, title=''):
    """

    :param matrix_file:
    :param title:
    :return:
    """
    plt.figure()
    df = pd.read_excel(matrix_file)
    ax = plt.axes()
    sns.color_palette("YlOrBr", as_cmap=True)
    sns.heatmap(df , ax=ax)
    ax.set_title(title)
    return


def store(element, file_type='doc', file_name=''):

    if file_type == 'doc':
        file = open(file_name + '.txt', 'w')
        file.write(str(element))
        file.close()
        return
    elif file_type == 'excel':
        df = pd.DataFrame(element)
        df.to_excel(excel_writer=file_name + '.xlsx')
        return


def store_polymer(polymer, filename):

    file = open(filename + '.xyz', 'w')
    # if not hasattr(polymer, ""):
    #     raise TypeError(f"input given is type {type(polymer)} but required type is PolymerClass")
    file.write(str(len(polymer)) + '\n' + 'DNA Polymer' + '\n')
    print(polymer.polymer)
    for bead in polymer.polymer:
        pos = polymer.polymer[bead]()
        file.write('H    ')
        for cord in pos:
            if cord >= 0:
                file.write(f" {cord:.5f}    ")
            else:
                file.write(f"{cord:.5f}    ")
        file.write('\n')
    file.close()


def check_avg_bead_dist(polymer):

    fig=plt.figure()
    ax = plt.axes()
    capture = []
    for bead in polymer:
        for bead2 in polymer:
            if int(bead2) > int(bead):
                capture.append(distance(polymer[bead](), polymer[bead2]()))
    ax.hist(capture, bins=np.arange(0, 10, 0.5))
    return


def cal_centre_of_mass(beads):

    x, y, z = 0, 0, 0
    for bead in beads:
        x += bead()[0]
        y += bead()[1]
        z += bead()[2]
    x, y, z = x/len(beads), y/len(beads), z/len(beads)
    return [x, y, z]


def cal_radius_of_gyration(beads):

    com = cal_centre_of_mass(beads)
    return sum([distance(com, bead()) for bead in beads])/len(beads)


