import logging
import time

import BeadClass as Bc
import EnergyBondClass as Ebc
import numpy as np
from matplotlib import pyplot as plt
import random as rdm
import HelperFunction1 as Hf1
from collections import OrderedDict
import pandas as pd
import string
import DataStorageCode as Dsc

logging.basicConfig(level=logging.INFO)
alphabet_string = string.ascii_uppercase
extra_strings = ()
alphabet_list = list(alphabet_string)


class Polymers:

    def __init__(self,
                 polymers,
                 bead_size=1.0,
                 bead_id_list=False,
                 bead_specific_class_list=False,
                 spherical_constraint=([0, 0, 0], 0.6),
                 **kwargs):

        self.checklist = {'PolDev': False,
                          'SpringPot': False,
                          'ExcludePot': False,
                          'SpherePot': False,
                          'SphereBindPot': False,
                          'BreakFreq': False,
                          'TransFreq': False
                          }

        logging.info("Developing the Polymer... \n")

        global alphabet_list
        self.new_name_list = alphabet_list
        self.breaks_counter = 0
        if isinstance(polymers, dict):
            if isinstance(polymers[list(polymers.keys())[0]][0], int):
                for pol, seq in polymers.items():
                    polymers[pol] = [np.array([0, 0, 0])]*(seq[0])
        else:
            raise TypeError(
                f'polymers input should be a dictionary with keys representing the names of the polymers and values'
                f'representing the number of beads in the polymer.'
                f'If the positions are known, then just add the file to the polymer input'
            )

        if bead_specific_class_list is True:
            if not isinstance(bead_specific_class_list, dict):
                raise TypeError(
                    f'bead_specific_class_list should either be false or should be a dictionary with keys representing'
                    f'polymer names and values representing a list of bead classes for that polymer.'
                    f'The type provided is {type(bead_specific_class_list)}'
                )
        else:
            bead_specific_class_list = dict()
            for pol, seq in polymers.items():
                bead_specific_class_list[pol] = ['None']*len(seq)

        if bead_id_list is True:
            if not isinstance(bead_id_list, dict):
                raise TypeError(
                    f'bead_id_list should either be false or should be a dictionary with keys representing'
                    f'polymer names and values representing a list of bead_ids for that polymer.'
                    f'The type provided is {type(bead_id_list)}'
                )
        else:
            bead_id_list = dict()
            bead_counter = 1
            for pol, seq in polymers.items():
                bead_id_list[pol] = [i for i in range(bead_counter, bead_counter + len(seq), 1)]
                bead_counter += len(seq)

        self.polymers = OrderedDict()
        for pol, seq in polymers.items():
            self.polymers[pol] = [Bc.Bead(polymers[pol][i], bead_specific_class_list[pol][i],
                                          identity=bead_id_list[pol][i],
                                          diameter=bead_size) for i in range(len(seq))]
        self.gen_size = self.polymers[list(self.polymers.keys())[-1]][-1].id
        self.bead_size = bead_size
        self._spherical_radius, self._spherical_centre = (self.gen_size**spherical_constraint[1])*self.bead_size, \
                                                         np.array(spherical_constraint[0])

        self.spring_potentials = []
        self.non_overlap_potentials = []
        self.sphere_limit_potentials = []
        self.bead_positioning_potentials = []
        self.bead_bead_potentials1 = []
        self.break_freq, self.trans_freq = np.full((self.gen_size, self.gen_size), 0.003), np.full((self.gen_size,
                                                                                                     self.gen_size),
                                                                                                    0.99)
        self.Trans_history, self.Break_history = np.zeros((self.gen_size, self.gen_size)), np.zeros((self.gen_size,
                                                                                                    self.gen_size))
        self.all_polymers, self.brokenBeadsA, self.brokenBeadsB = [], [], []
        self.all_broken_lists = [self.brokenBeadsA, self.brokenBeadsB]

        connections = {
            'SpringPot': self.spring_potentials,
            'PolDev': self.polymers,
            'ExcludePot': self.non_overlap_potentials,
            'SpherePot': self.sphere_limit_potentials,
            'SphereBindPot': self.bead_positioning_potentials,
            'BreakFreq': self.break_freq,
            'TransFreq': self.trans_freq
        }

        logging.info("Energy matrices Initialized... \n")

        self._system_energy, self._polymer_size = [], []
        self._energy_average, self._pol_size_average = [], []

        for name, value in kwargs:
            if name in self.checklist:
                connections[str(name)] = value
                if str(name) == 'BreakFreq' or str(name) == 'TransFreq':
                    if not isinstance(value, np.ndarray):
                        raise TypeError(
                            f'Given argument for {name} is of type {type(value)}, type should be np.ndarray'
                        )
                    if value.size != (self.gen_size, self.gen_size):
                        raise ValueError(
                            f'The size of the array for {name} given is {value.size}. Size required is '
                            f'{(self.gen_size, self.gen_size)}.'
                        )
                elif str(name) == 'PolDev':
                    if not isinstance(value, (OrderedDict, dict)):
                        raise ValueError(
                            f'type of {name} should be an OrderedDict or dict, but was given as {type(value)}'
                        )
                if str(name) in ['SpringPot', 'ExcludePot', 'SpherePot', 'SphereBindPot']:
                    if not isinstance(value, (tuple, list)):
                        raise TypeError(
                            f'type of {name} should be a sequence like a list, but type given is {type(value)}.'
                        )
                logging.info(f'User defined {str(name)} array was updated... \n')
                self.checklist[str(name)] = True
            else:
                raise ValueError(
                    f'Given argument not valid. Valid extra arguments include the following {self.checklist.keys()}'
                )
        self.all_polymers = list(self.polymers.keys())

    def __len__(self):

        return self.gen_size

    def create_random_polymer(self):

        """
        uses the monomer class to add N monomers in a random
        position within the polymer and identifies the
        difference between the first and last monomer of the
        polymer

        :return: None, updates the self.polymer_pos
        """
        if self.checklist['PolDev'] is True:
            logging.info(f'polymer positions have already been given. Random walk is not going to be performed.')
            return
        logging.info("Polymer is being developed... \n")
        new_position = self._spherical_centre  # initialize position
        for pol, beads in self.polymers.items():  # Move through all of the beads
            for i in range(len(beads)):
                if beads[i].id == 1:
                    beads[i].pos = new_position  # If it is the first bead then allow it to stay in the centre of the
                    # sphere
                    continue
                while Hf1.overlap(self.polymers, new_position, self.bead_size) is True:  # Check if there is overlap
                    new_position = beads[i].setup_bead(beads[i-1])  # Change into a new position
                beads[i].pos = new_position  # Confirm position when overlap is not found.
        logging.info(f'Polymers have been given a random configurations')
        self.checklist['PolDev'] = True
        return

    def initialize_spring_energy(self, constant=30):

        if self.checklist['SpringPot'] is True:
            logging.info(f'spring potentials have already been added. No change done')
            return
        for pol, beads in self.polymers.items():  # loop over all polymers
            for i in range(len(beads) - 1):  # loop over all beads of the polymer in question
                bond = Ebc.SpringBond(beads[i], beads[i+1], constant, self.bead_size)  # Create a SpringBond class
                length = Hf1.distance(beads[i](), beads[i+1]())
                bond.make_change(length)  # Updates the energy of the bond after checking distance between the beads
                self.spring_potentials.append(bond)  # Updates the SpringBond
                # matrix
        logging.info(f'SpringEnergies have been updates')
        self.checklist['SpringPot'] = True
        return

    def initialize_exclusion_energy(self, constant=1):

        if self.checklist['ExcludePot'] is True:
            logging.info(f'Exclusion potentials have already been added. No change done')
            return
        all_beads = []
        for beads in self.polymers.values():  # loop over all total beads
            all_beads += beads
        for i in range(len(all_beads)):  # loop over all beads
            for j in range(i+1, len(all_beads)):
                bond = Ebc.RepLjBond(all_beads[i], all_beads[j], constant, self.bead_size)  # Create an RepLJBond
                # class
                length = Hf1.distance(all_beads[i](), all_beads[j]())
                bond.make_change(length)  # Updates the energy of the bond after checking distance
                # between the beads
                self.non_overlap_potentials.append(bond)  # Updates
                # the Exclusion matrix
        logging.info(f'RJ exclusion Energies have been updates')
        self.checklist['ExcludePot'] = True
        return

    def initialize_sphere_limit_energy(self, constant=100):

        if self.checklist['SpherePot'] is True:
            logging.info(f'Sphere limiting potentials have already been added. No change done')
            return
        all_beads = []
        for beads in self.polymers.values():  # loop over all total beads
            all_beads += beads
        for bead in all_beads:
            bond = Ebc.NucleusConfinement(bead, centre=Bc.Bead(self._spherical_centre, identity=0),
                                          energy_constant=constant, radius=self._spherical_radius)
            length = Hf1.distance(bead(), self._spherical_centre)
            bond.make_change(length)  # Updates the energy of the bond after checking distance
            self.sphere_limit_potentials.append(bond)
        logging.info(f'Bead type specific nuclear binding has been initialized... \n')
        self.checklist['SpherePot'] = True
        return

    def initialize_nucleus_binding_energy(self, constant=10, bead_type='Cen'):

        if self.checklist['SphereBindPot'] is True:
            logging.info(f'Specific nuclear binding potentials have already been added. No change done... \n')
            return
        all_beads = []
        for beads in self.polymers.values():  # loop over all total beads
            all_beads += beads
        for bead in all_beads:
            if bead.bead_class != bead_type:
                continue
            bond = Ebc.NucleusBinding(bead, centre=Bc.Bead(self._spherical_centre, identity=0),
                                      energy_constant=constant, radius=self._spherical_radius)
            length = Hf1.distance(bead(), self._spherical_centre)
            bond.make_change(length)  # Updates the energy of the bond after checking distance
            self.bead_positioning_potentials.append(bond)
        logging.info(f'Spherical Limit For the beads have been placed... \n')
        self.checklist['SphereBindPot'] = True
        return

    def attractive_forces(self, constant=1):

        all_beads = []
        for beads in self.polymers.values():  # loop over all total beads
            all_beads += beads
        for i in range(len(all_beads)):  # loop over all beads
            for j in range(i + 1, len(all_beads)):
                bond = Ebc.AttLjBond(all_beads[i], all_beads[j], constant, self.bead_size)  # Create an RepLJBond
                # class
                length = Hf1.distance(all_beads[i](), all_beads[j]())
                bond.make_change(length)  # Updates the energy of the bond after checking distance
                # between the beads
                self.bead_bead_potentials1.append(bond)  # Updates
                # the Exclusion matrix
        logging.info(f'RJ exclusion Energies have been updates')
        return

    def move(self, breaks=False):

        for dummy_i in range(self.gen_size):
            if breaks:
                self._bond_break()
            polymer = rdm.choice(list(self.polymers.keys()))
            bead = rdm.choice(self.polymers[polymer])  # choose a random bead
            old_pos = bead()
            new_pos = bead.try_move()  # try moving the bead to a position
            # capture all the bonds associated to the bead of interest
            all_bonds = self._capture_all_bonds(bead.id)
            old_energy = sum([bond.energy for bond in all_bonds])
            new_energy = sum([bond.try_length(Hf1.distance(new_pos, other.pos)) for bond in all_bonds for other
                              in bond.bond_beads if other != bead])
            # capture the new bond energy after trying the new move
            allow_length = Hf1.check_length(old_energy, new_energy)  # check if the move length is acceptable by prob.
            if allow_length:
                bead.pos = new_pos  # move the bead
                change_energies = [bond.make_change(Hf1.distance(new_pos, other.pos)) for bond in all_bonds for
                                   other in bond.bond_beads if other != bead]  # change all the energy values and
                # lengths of the bonds
            if breaks:
                self._bond_fix()
        # self._update_energy()
        # self._update_size()
        return

    def simulate(self, sim_number, breaks=False):

        logging.info(f'Simulation has begun with {sim_number} monte carlo steps \n')
        if not breaks:
            logging.info(f'Average waiting time is approximately {(10**(-5))*(self.gen_size**2.15)*(sim_number)}'
                         f' seconds.\n')
        for idx in range(sim_number):
            self.move(breaks)
        logging.info(f'Simulation is complete. \n')

        return

    def _bond_break(self):

        polymer = rdm.choice(list(self.polymers.keys()))
        selected_beads = list(enumerate(self.polymers[polymer]))
        bead_idx, bead = rdm.choice(selected_beads)
        if (bead_idx, bead) == selected_beads[-1]:
            return
        r = rdm.random()
        bead1, bead2 = bead, selected_beads[bead_idx + 1][1]

        if r < self.break_freq[bead1.id - 1, bead2.id - 1]:
            self.polymers['Brk' + str(self.breaks_counter)] = self.polymers[polymer][:bead_idx + 1]
            self.polymers[polymer] = self.polymers[polymer][bead_idx + 1:]
            [self.spring_potentials.remove(bond) for bond in self.spring_potentials if bond.id == (bead1.id, bead2.id)
             or bond.id == (bead2.id, bead1.id)]
            self.collect_beads(bead1, bead2)
            self.breaks_counter += 1

        self.Break_history[bead2.id - 1, bead1.id - 1] += 1
        self.Break_history[bead1.id - 1, bead2.id - 1] += 1
        return

    def collect_beads(self, bead1, bead2):

        self.brokenBeadsA.append(bead1)
        self.brokenBeadsB.append(bead2)

        return

    def _bond_fix(self):

        """
        Tries to simulate the translocation process between two double strand breaks in a polymer
        :return: None, makes changes to self.polymers, self.brokenBeadsA, self.brokenBeadB and self.spring_potentials
        """
        # ------------------------------------------------------------------------------------------------------------
        if not self.brokenBeadsA:
            return
        contacts = Hf1.capture_contacts(self.non_overlap_potentials, self.gen_size, 1.2)  # array, contact map of beads
        bead_idx = rdm.choice(range(len(self.brokenBeadsA)))  # integer, selects an index from the broken beads
        bead1, bead2 = self.brokenBeadsA[bead_idx], self.brokenBeadsB[bead_idx]  # Beads, captures the beads (it and its
        # previous partner)
        contacts1, contacts2 = contacts[int(bead1.id) - 1], contacts[int(bead2.id) - 1]  # 1D array, captures contact
        # row of the beads
        possible_1, possible_2 = np.where(contacts1 == 1), np.where(contacts2 == 1)  # 1D array, captures only
        # beads that have contacts with bead in question
        possible_1, possible_2 = [i+1 for i in possible_1[0]], [i+1 for i in possible_2[0]]
        poss_beads1, poss_beads2 = [], []
        for beadset in [possible_1, possible_2]:
            for p in beadset:
                for bead in self.brokenBeadsA:
                    if p == bead.id:
                        if beadset == possible_1:
                            poss_beads1.append(bead)
                        else:
                            poss_beads2.append(bead)
                for bead in self.brokenBeadsB:
                    if p == bead.id:
                        if beadset == possible_1:
                            poss_beads1.append(bead)
                        else:
                            poss_beads2.append(bead)
        # ------------------------------------------------------------------------------------------------------------

        possibilities = []

        # ------------------------------------------------------------------------------------------------------------

        # To Capture all possible pairs of contacts. If a and b beads were selected, then the code search for a set or
        # pair of broken beads previously together, but also such that a in contact with c and b in contact with d or
        # a in contact with d and b in contact with c. Here, a and b are the initial beads selected by the previous set
        # of codes.

        for bead in poss_beads1:
            if bead in self.brokenBeadsA:
                if self.brokenBeadsB[self.brokenBeadsA.index(bead)] in poss_beads2:
                    idx = self.brokenBeadsA.index(bead)
                    possibilities.append([(bead1, bead, self.brokenBeadsA, idx), (bead2, self.brokenBeadsB[idx],
                                                                                  self.brokenBeadsB, idx)])
            if bead in self.brokenBeadsB:
                if self.brokenBeadsA[self.brokenBeadsB.index(bead)] in poss_beads2:
                    idx = self.brokenBeadsB.index(bead)
                    possibilities.append([(bead1, bead, self.brokenBeadsB, idx), (bead2, self.brokenBeadsA[idx],
                                                                                  self.brokenBeadsA, idx)])

        # ------------------------------------------------------------------------------------------------------------

        # Main lines of code which help manipulate the self.polymers dictionary such that the translocation takes place
        # Here, one key gets deleted as two polymers are now combined to one longer polymer.
        if not possibilities:
            return
        selection = rdm.choice(possibilities)
        if selection[0][0] == selection[1][1] and selection[0][1] == selection[1][0]:
            selection.pop(1)
        for trans in selection:
            set1, set2 = [], []
            # For selecting the orientation of the beads to be bound, for appending purposes.
            for polymer, bead_set in self.polymers.items():
                if trans[0].id == bead_set[0].id:
                    orientation = [i for i in bead_set[::-1]]
                    set1 = (polymer, orientation)
                elif trans[0].id == bead_set[-1].id:
                    set1 = (polymer, bead_set)
                if trans[1].id == bead_set[0].id:
                    set2 = (polymer, bead_set)
                elif trans[1].id == bead_set[-1].id:
                    orientation = [i for i in bead_set[::-1]]
                    set2 = (polymer, orientation)
            self.polymers[set1[0]] = set1[1] + set2[1]
            del self.polymers[set2[0]]

            # Adds a new spring potential between the two beads to be bonded
            self.spring_potentials.append(Ebc.SpringBond(trans[0], trans[1]))

            # Removes the beads from the broken beads lists.
            [broken_list.pop(trans[3]) for broken_list in self.all_broken_lists]
            self.Trans_history[trans[0].id - 1, trans[1].id - 1] += 1
            self.Trans_history[trans[1].id - 1, trans[0].id - 1] += 1
        return

        # ------------------------------------------------------------------------------------------------------------

    def _update_energy(self):

        all_energies = self.spring_potentials + self.non_overlap_potentials + self.sphere_limit_potentials + \
                            self.bead_positioning_potentials + self.bead_bead_potentials1
        if len(self._system_energy) > 10000:
            self._system_energy, self._energy_average = self._system_energy[1:], self._energy_average[1:]
        self._system_energy.append(sum([bond() for bond in all_energies]))
        self._energy_average.append(sum(self._system_energy)/len(self._system_energy))
        return

    def _update_size(self):

        if len(self._polymer_size) > 10000:
            self._polymer_size, self._pol_size_average = self._polymer_size[1:], self._pol_size_average[1:]
        self._polymer_size.append(self.get_radius_of_gyration())
        self._pol_size_average.append(sum(self._polymer_size)/len(self._polymer_size))
        return

    def plot_energy(self, ax, other_naming=''):

        # fig = plt.figure()
        # ax3 = plt.axes()
        plt.title('Energy Graph')
        ax.plot(self._system_energy, '--', label='Instantaneous Energy Value ' + str(other_naming))
        ax.plot(self._energy_average, label='Average Energy Value ' + str(other_naming), linewidth=4)
        return

    def plot_polymer_size(self, ax, other_naming=''):

        # fig = plt.figure()
        plt.title('Radius of Gyration Graph')
        ax.plot(self._polymer_size, '--', label='Instantaneous RoG Value ' + str(other_naming))
        ax.plot(self._pol_size_average, label='Average RoG Value ' + str(other_naming), linewidth=4)
        return

    def class_update(self, beads, bead_type='Cen'):

        for pol, bead_list in self.polymers.items():
            for bead in bead_list:
                if bead.id in beads:
                    bead.bead_class = bead_type
        return

    #
    def other_bead_bead_attractions(self, bead_type='Cen', constant=10):

        for pol1, beads1 in self.polymers.items():
            for pol2, beads2 in self.polymers.items():
                for bead1 in beads1:
                    for bead2 in beads2:
                        if bead1.bead_class == bead_type and bead2.bead_class == bead_type:
                            if bead1.id > bead2.id:
                                bond = Ebc.SpringBond(bead1, bead2, constant, self.bead_size)
                                length = Hf1.distance(bead1(), bead2())
                                bond.make_change(length)
                                self.bead_bead_potentials1.append(bond)
        return

    def _capture_all_bonds(self, bead_name):

        all_bonds = Hf1.capture_bonds(bead_name, self.spring_potentials) + Hf1.capture_bonds(bead_name,
                                                                                             self.non_overlap_potentials
                                                                                             ) + Hf1.capture_bonds\
            (bead_name, self.sphere_limit_potentials) + Hf1.capture_bonds(bead_name, self.bead_positioning_potentials)\
        + Hf1.capture_bonds(bead_name, self.bead_bead_potentials1)

        return all_bonds

    def show_polymer(self, ax):

        for pol, beads in self.polymers.items():
            x, y, z = [], [], []
            for bead in beads:
                x.append(bead()[0])
                y.append(bead()[1])
                z.append(bead()[2])
            plt.plot(x, y, z)

        return

    def _get_centre_of_mass(self):

        x, y, z = 0, 0, 0
        for pol, beads in self.polymers.items():
            for bead in beads:
                x += bead()[0]
                y += bead()[1]
                z += bead()[2]
        x, y, z = x/self.__len__(), y/self.__len__(), z/self.__len__()
        return [x, y, z]

    def get_radius_of_gyration(self):

        com = self._get_centre_of_mass()
        return sum([Hf1.distance(com, bead()) for pol, beads in self.polymers.items() for bead in beads])/self.__len__()
    #
    # def __call__(self):
    #
    #     return [self.polymer]
    #
    # def __repr__(self):
    #
    #     string = ''
    #     polymer = self.polymer.items()
    #     for key, value in polymer:
    #         string = string + str(key) + '->' + str(value) + '\n'
    #     return str(string)
    #

    def print_nucleus(self, ax):

        Hf1.plt_sphere([self._spherical_centre], ax, [self._spherical_radius], 0.2, 'black')
        return
    #
    # def __class__(self: 'PolymerClass'):
    #     return "ClassPolymer"


polymer1 = Polymers({'chr1': [50],
                     'chr2': [50],
                     'chr3': [50],
                     'chr4': [50]
                     })
plt.figure()
ax = plt.axes(projection='3d')
polymer1.create_random_polymer()
polymer1.initialize_spring_energy()
polymer1.initialize_exclusion_energy()
polymer1.initialize_sphere_limit_energy()
# polymer1.class_update([40, 80, 120, 170])
polymer1.class_update([30, 31, 32, 33, 34, 35, 36, 75, 76, 77, 78, 79, 80, 81, 120, 121, 122, 123, 124, 125,
                       176, 177, 178, 179, 180, 181])
polymer1.other_bead_bead_attractions()
polymer1.initialize_nucleus_binding_energy()
print(polymer1.bead_bead_potentials1)
a = time.time()
timestep = 0
for snap in range(1000):
    polymer1.simulate(100)
    timestep += 1000
    print(timestep)
    Dsc.store_positions(polymer1, 'YeastGenomeExamplePos', timestep)
    Dsc.store_bead_char(polymer1, 'YeastGenomeExampleBeads', timestep)
    Dsc.store_bonds(polymer1, 'YeastGenomeExampleBonds', timestep)
    Dsc.store_pol_char(polymer1, 'YeastGenomeExampleRoG', timestep)
b = time.time()
print(b - a)
# polymer1.show_polymer(ax)
# Hf1.plt_sphere([polymer1._spherical_centre], ax, list_radius=[polymer1._spherical_radius], alpha=0.5, color='blue')
# plt.figure()
# ax2 = plt.axes()
# polymer1.plot_energy(ax2)
# plt.figure()
# ax3 = plt.axes()
# polymer1.plot_polymer_size(ax3)
# print(polymer1.Break_history, '\n', polymer1.Trans_history)
# print(polymer1.polymers)
# # Dsc.store_positions(polymer1, 'Test_File')
# # Hf1.store_polymer(polymer1, 'polymer1_test_file')
# # # print(polymer1.polymer)
# plt.show()










