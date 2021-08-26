"""
EnergyBondClass
-----------------------

This module contains 4 different kinds of interactions that can affect the beads in the yeast genome.

Bead - Bead interactions include: SpringBond and RepLjBond

The module also contains 2 bead - nucleus interactions that allow placement of beads at a specific location

Bead - nucleus interactions include: NuclearConfinement (To confine polymers within the nucleus)
and NuclearBinding (allows certain beads to be bound or remain close to the nuclear periphery (Centromeres
in the YeastGenome).

Additional Bonds can be added provided they have similar codes as compared to the classes below.
The Bonds should also be added separately in the simulation class in the YeastGenome module.

All Bonds have 4 parameters: interacting Beads 1 and 2, potential parameter and natural length scale of the bond.
For three or more bead interactions, make sure the NewBond.id method is a tuple of all beads involved.
"""

import HelperFunction1 as Hf1
import numpy as np
import BeadClass as Bc


class SpringBond:

    """
    The SpringBond Class
    """

    def __init__(self,
                 bead1,
                 bead2,
                 energy_constant=10,
                 bond_length=1.0):
        """

        :param bead1: First bead in the bond (Bead)
        :param bead2: Second bead in the bond (Bead)
        :param energy_constant: Energy Parameter for the SpringBond (float/int)
        :param bond_length: The natural bond length of the SpringBond (float)
        """

        if bond_length < 0 or type(bond_length) != float:
            raise ValueError('bond_length not possible')
        if isinstance(bead1, Bc.Bead) and isinstance(bead2, Bc.Bead):
            self._bead1, self._bead2 = bead1, bead2
            self.bond_beads = self._bead1, self._bead2
        else:
            raise TypeError(f'Bead identities should be of type Bead, but found to be {type(bead1)}, {type(bead2)}')
        if not (isinstance(energy_constant, float) or isinstance(energy_constant, int)):
            raise TypeError(f'energy constant should be type int or float, but type found was {type(energy_constant)}')
        self._ene_con = energy_constant
        self.id = (bead1.id, bead2.id)
        self.energy, self.type = 0, '__Spring__'
        if not isinstance(bond_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(bond_length)}')
        self._bond_length, self.current_length = bond_length, Hf1.distance(bead1(), bead2())
        self.make_change(self.current_length)

    def _function(self, length):

        """
        Energy Calculation
        :param length: the length to check the energy value of the bond (float)
        :return: Energy Value at bond length 'length'
        """

        return (self._ene_con * 0.5) * ((length - self._bond_length) ** 2)

    def try_length(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length, but does not update to this length
        :param new_length: (float)
        :return: new energy of the bond (predicted)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        if new_length < 0:
            raise ValueError('Impossible lengths')
        return self._function(new_length)

    def make_change(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length and updates bond length and energy
        :param new_length: (float)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        self.current_length = new_length
        if new_length < 0:
            raise ValueError('Impossible lengths')
        self.energy = self._function(new_length)

    def __call__(self):

        return self.energy

    def __str__(self):

        return f'{self._bead1.id} {self._bead2.id} {self._ene_con} {self._bond_length} {self.current_length} {self.type}'

    def __repr__(self):

        return f'SpringBond: {self.id}'


class AttLjBond:

    def __init__(self, bead1, bead2, energy_constant=1.0, bond_length=1.0):

        """

        :param bead1: First bead in the bond (Bead)
        :param bead2: Second bead in the bond (Bead)
        :param energy_constant: Energy Parameter for the SpringBond (float/int)
        :param bond_length: The natural bond length of the SpringBond (float)
        """

        if bond_length < 0 or type(bond_length) != float:
            raise ValueError('bond_length not possible')
        if isinstance(bead1, Bc.Bead) and isinstance(bead2, Bc.Bead):
            self._bead1, self._bead2 = bead1, bead2
            self.bond_beads = self._bead1, self._bead2
        else:
            raise TypeError(f'Bead identities should be of type Bead, but found to be {type(bead1)}, {type(bead2)}')
        if not (isinstance(energy_constant, float) or isinstance(energy_constant, int)):
            raise TypeError(f'energy constant should be type int or float, but type found was {type(energy_constant)}')
        self._ene_con = energy_constant
        self.id = (bead1.id, bead2.id)
        if not isinstance(bond_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(bond_length)}')
        self._bond_length, self.current_length = bond_length, Hf1.distance(bead1(), bead2())
        self._least_energy_length = (2 ** (1/6)) * self._bond_length
        self.make_change(self.current_length)
        self.type = '__attractive__'

    def _function(self, length):

        """
        Energy Calculation
        :param length: the length to check the energy value of the bond (float)
        :return: Energy Value at bond length 'length'
        """

        return (4 * self._ene_con * ((self._bond_length / length) ** 12 -
                                     (self._bond_length / length) ** 6)) + self._ene_con

    def try_length(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length, but does not update to this length
        :param new_length: (float)
        :return: new energy of the bond (predicted)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        if new_length < 0:
            raise ValueError('Impossible lengths')
        if new_length > self._least_energy_length:
            return self._function(new_length)
        else:
            return 0

    def make_change(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length and updates bond length and energy
        :param new_length: (float)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        self.current_length = new_length
        if new_length < 0:
            raise ValueError('Impossible lengths')
        elif new_length == 0:
            self.energy = np.inf
        elif new_length > self._least_energy_length:
            self.energy = self._function(new_length)
        else:
            self.energy = 0
        return

    def __call__(self):

        return self.energy

    def __str__(self):

        return f'{self._bead1.id} {self._bead2.id} {self._ene_con} {self._bond_length} {self.current_length} {self.type}'

    def __repr__(self):

        return f'AttLJBond: {self.id}'


class RepLjBond:

    """
    Repulsive part of the LJ Potential
    """

    def __init__(self, bead1, bead2, energy_constant=1.0, bond_length=1.0):

        """

        :param bead1: First bead in the bond (Bead)
        :param bead2: Second bead in the bond (Bead)
        :param energy_constant: Energy Parameter for the SpringBond (float/int)
        :param bond_length: The natural bond length of the SpringBond (float)
        """

        if bond_length < 0 or type(bond_length) != float:
            raise ValueError('bond_length not possible')
        if isinstance(bead1, Bc.Bead) and isinstance(bead2, Bc.Bead):
            self._bead1, self._bead2 = bead1, bead2
            self.bond_beads = self._bead1, self._bead2
        else:
            raise TypeError(f'Bead identities should be of type Bead, but found to be {type(bead1)}, {type(bead2)}')
        if not (isinstance(energy_constant, float) or isinstance(energy_constant, int)):
            raise TypeError(f'energy constant should be type int or float, but type found was {type(energy_constant)}')
        self._ene_con = energy_constant
        self.id = (bead1.id, bead2.id)
        if not isinstance(bond_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(bond_length)}')
        self._bond_length, self.current_length = bond_length, Hf1.distance(bead1(), bead2())
        self._least_energy_length = (2 ** (1/6)) * self._bond_length
        self.make_change(self.current_length)
        self.type = '__repulsive__'

    def _function(self, length):

        """
        Energy Calculation
        :param length: the length to check the energy value of the bond (float)
        :return: Energy Value at bond length 'length'
        """

        return (4 * self._ene_con * ((self._bond_length / length) ** 12 -
                                     (self._bond_length / length) ** 6)) + self._ene_con

    def try_length(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length, but does not update to this length
        :param new_length: (float)
        :return: new energy of the bond (predicted)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        if new_length < 0:
            raise ValueError('Impossible lengths')
        if new_length < self._least_energy_length:
            return self._function(new_length)
        else:
            return 0

    def make_change(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length and updates bond length and energy
        :param new_length: (float)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        self.current_length = new_length
        if new_length < 0:
            raise ValueError('Impossible lengths')
        elif new_length == 0:
            self.energy = np.inf
        elif new_length < self._least_energy_length:
            self.energy = self._function(new_length)
        else:
            self.energy = 0
        return

    def __call__(self):

        return self.energy

    def __str__(self):

        return f'{self._bead1.id} {self._bead2.id} {self._ene_con} {self._bond_length} {self.current_length} {self.type}'

    def __repr__(self):

        return f'RepLJBond: {self.id}'


class NucleusConfinement:

    def __init__(self,
                 bead,
                 centre=Bc.Bead(position=np.array([0, 0, 0], dtype=float), identity=0),
                 energy_constant=30,
                 radius=1.0):

        if radius < 0:
            raise ValueError('bond_length not possible')
        if isinstance(bead, Bc.Bead) and isinstance(centre, Bc.Bead):

            self._centre = centre
            self._bead1, self._bead2 = bead, self._centre
            self.bond_beads = self._bead1, self._bead2

        else:
            raise TypeError(f'Bead identities should be of type Bead, but found to be {type(bead)}, {type(centre)}')
        if not (isinstance(energy_constant, float) or isinstance(energy_constant, int)):
            raise TypeError(f'energy constant should be type int or float, but type found was {type(energy_constant)}')
        if not isinstance(radius, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(radius)}')
        self._radius = radius
        self.current_length = Hf1.distance(bead(), centre())
        self._ene_con = energy_constant
        self.id = (self._bead1.id, self._bead2.id)
        self.energy, self.type = 0, '__NucleusConfinement__'
        self.make_change(self.current_length)

    def _function(self, length):

        """
        Energy Calculation
        :param length: the length to check the energy value of the bond (float)
        :return: Energy Value at bond length 'length'
        """

        return (self._ene_con * 0.5) * ((length - self._radius) ** 2)

    def try_length(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length, but does not update to this length
        :param new_length: (float)
        :return: new energy of the bond (predicted)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        if new_length < 0:
            raise ValueError('Impossible lengths')
        if new_length > self._radius:
            return self._function(new_length)
        else:
            return 0

    def make_change(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length and updates bond length and energy
        :param new_length: (float)
        """

        if not isinstance(new_length, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(new_length)}')
        if new_length < 0:
            raise ValueError('Impossible lengths')
        self.current_length = new_length
        if new_length > self._radius:
            self.energy = self._function(new_length)
        else:
            self.energy = 0

    def __call__(self):

        return self.energy

    def __str__(self):

        return f'{self._bead1.id} {self._bead2.id} {self._ene_con} {self._radius} {self.current_length} {self.type}'

    def __repr__(self):

        return f'NucleusConfinementEnergy: {self.id}'


class NucleusBinding:

    def __init__(self,
                 bead,
                 centre=Bc.Bead(position=np.array([0, 0, 0], dtype=float), identity=0),
                 energy_constant=50,
                 radius=1.0):

        if radius < 0:
            raise ValueError('bond_length not possible')
        if isinstance(bead, Bc.Bead) and isinstance(centre, Bc.Bead):
            self._bead1, self._bead2 = bead, centre
            self._centre = centre
            self.bond_beads = self._bead1, self._bead2

        else:
            raise TypeError(f'Bead identities should be of type Bead, but found to be {type(bead)}, {type(centre)}')
        if not (isinstance(energy_constant, float) or isinstance(energy_constant, int)):
            raise TypeError(f'energy constant should be type int or float, but type found was {type(energy_constant)}')
        if not isinstance(radius, float):
            raise TypeError(f' Bond length should be of type float, but found to be of type {type(radius)}')
        self._radius = radius
        self.current_length = Hf1.distance(bead(), centre())
        self._ene_con = energy_constant
        self.id = (self._bead1.id, self._bead2.id)
        self.energy, self.type = 0, '__NucleusBinding__'
        self.make_change(self.current_length)

    def _function(self, length):

        """
               Energy Calculation
               :param length: the length to check the energy value of the bond (float)
               :return: Energy Value at bond length 'length'
               """

        return (self._ene_con * 0.5) * ((length - self._radius) ** 2)

    def try_length(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length, but does not update to this length
        :param new_length: (float)
        :return: new energy of the bond (predicted)
        """

        return (self._ene_con * 0.5) * ((new_length - self._radius) ** 2)

    def make_change(self, new_length):

        """
        Calculates the Energy Value of the bond at new_length and updates bond length and energy
        :param new_length: (float)
        """

        self.energy = (self._ene_con * 0.5) * ((new_length - self._radius) ** 2)
        self.current_length = self._radius - new_length

    def __call__(self):

        return self.energy

    def __str__(self):

        return f'{self._bead1.id} {self._bead2.id} {self._ene_con} {self._radius} {self.current_length} {self.type}'

    def __repr__(self):

        return f'NucleusBindingEnergy: {self.id}'

