""" Code to load all data into the simulation script."""


def load_polymer(file):

    file = open(file, 'r')
    polymers = {}
    positions = []
    lines = file.read()
    list_of_lines = lines.splitlines()
    file.close()
    for line in range(len(list_of_lines) - 1, 1, -1):
        if list_of_lines[line] == '':
            pass
        else:
            beads = int(list_of_lines[0])
            positions = list_of_lines[line - beads + 1:]
            print(len(positions))
            break
    for pos in positions[:-1]:
        bead = pos.split('    ')
        if bead[0] not in polymers.keys():
            polymers.setdefault(bead[0], [[float(bead[1]), float(bead[2]), float(bead[3])]])
        else:
            polymers[bead[0]].append([float(bead[1]), float(bead[2]), float(bead[3])])
    return polymers


a = '/home/anirudhjay/Book/YeastGenome/YeastGenomeExamplePos.xyz'
print(load_polymer(a))
