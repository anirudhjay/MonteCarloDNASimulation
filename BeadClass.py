"""
BeadClass
---------------

This module contains the monomer class called 'BeadClass' that represents a monomer and its characteristics.
It has 4 total parameters of which only 1 is mandatory:

Position: Should be an np array (3, 0)
BeadClass: The DNA type that the particular bead consists. (Check the dictionary 'color', to add more DNA element types.
Diameter: Size of the Bead.
Identity: a bead marker, that is used when multi-polymer simulations are done.
"""

import HelperFunction1 as Hf1
import random as r
import numpy as np


color = {'None': 'blue',
         'Cen': 'black',
         'Tel': 'grey',
         'Het': 'brown',
         'Euc': 'purple'}


class Bead:
    
    """
    A class 'Bead' that captures the characteristics of coarse-grained beads

    """

    def __init__(self,
                 position,
                 bead_class='None',
                 diameter=1.0,
                 identity=None):

        """
        :param position: The 3D position of the bead (array)
        :param bead_class: The type of DNA represented by the bead
                example inputs: Cen, Tel, Het, Euc etc. (string)
        :param diameter: size of the bead (float64)
        :param identity: A unique integer that identifies the bead in a larger polymer (int64)
        """

        global color
        if isinstance(position, np.ndarray):
            if len(position) == 3:
                self.pos = position
            else:
                raise ValueError(
                    f'length of the array should be 3, but array length is {len(position)}'
                )
        else:
            raise TypeError(
                f'pos parameter should be an numpy array, but type of pos given is {type(position)}'
            )
        if type(diameter) != float:
            raise TypeError(
                f'diameter type should be float, but type provided is {type(diameter)}'
            )
        self.size = diameter
        if bead_class in color.keys():
            self.bead_class = bead_class
            self._color = color[self.bead_class]
        else:
            raise ValueError(
                f' class type is not present in the genome. \n'
                f' current classes include the following \n'
                f' {color.keys()} \n'
                f' Class given is {bead_class}\n'
                f' Please go to addClass to add new genome specific classes'
            )
        if identity is None or isinstance(identity, int):
            self.id = identity
        else:
            raise TypeError(
                f' identity should be of type int or None, but type provided is {type(identity)}'
            )

    def move(self, pos):

        """
        :param pos: Provides a new 3D position for the bead (array)
        """

        if isinstance(pos, np.ndarray):
            if len(pos) == 3:
                self.pos = pos
            else:
                raise ValueError(
                    f'length of the array should be 3, but array length is {len(pos)}'
                )
        else:
            raise TypeError(
                f'pos parameter should be an numpy array, but type of pos given is {type(pos)}'
            )

    def try_move(self):

        """
        Function to randomly move a bead in 3D space from it's original position (Uniform random function used)
        :return: The new position
        """

        # calculates a small change in all three axes randomly
        dr = [r.uniform(-self.size / 2, self.size / 2) for dummy_i in range(3)]
        # adds the change to the new position
        new_pos = [self.pos[i] + dr[i] for i in range(3)]
        return new_pos

    def setup_bead(self, previous_bead):

        """
        Initial function to setup polymers before simulations
        :param previous_bead: The precursor bead within the polymer (Bead)
        :return: A new acceptable position for the bead at hand, this is such that they don't initially coincide
        """

        if isinstance(previous_bead, Bead):
            return Hf1.select_position(previous_bead(), self.size)
        else:
            raise TypeError(
                f'Input type should be a bead class but type is {type(previous_bead)}'
            )

    def update_class(self, bead_class):

        """

        :param bead_class: A DNA type within the genome polymer (String, should be present in dict(color))
        :return: Updates the bead's class
        """

        if bead_class in color.keys():
            self.bead_class = bead_class
            self._color = color[self.bead_class]
        else:
            raise ValueError(
                f' class type is not present in the genome. \n'
                f' current classes include the following \n'
                f' {color.keys()} \n'
                f' Please go to addClass to add new genome specific classes'
            )
        return

    def __str__(self):

        return f'{self.pos}, {self.bead_class}, {self.size}, {self.id}'

    def __call__(self):

        return self.pos

    def __repr__(self):

        return f'Bead: {self.id}'
