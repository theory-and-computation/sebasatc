import numpy as np

class Poscar:
    """ Python Class for VASP POSCAR """

    def __init__(self, file_path):
        self.file_path = file_path
        self.atom_types = None
        self.atom_counts = None
        self.lattice_constants = None
        self.lattice_matrix = None
        self.initial_positions = None
        self.after_rotation = None
        self.atom_id = None
        self.slab_center = None  # Center of the slab
        self.read_poscar()

    def read_poscar(self):
        """Reads the POSCAR file and extracts lattice information and atomic positions"""
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

        # Reading lattice vectors
        lattice_constants = np.array(list(map(float, lines[1].split())))
        lattice_matrix = np.array([list(map(float, line.split())) for line in lines[2:5]])

        # Reading atomic species and counts
        atom_types = lines[5].split()
        atom_counts = list(map(int, lines[6].split()))

        # Reading atomic positions (Cartesian coordinates)
        initial_positions = []
        for line in lines[8:]:
            if line.strip() == '':
                break
            initial_positions.append(list(map(float, line.split()[:3])))

        self.lattice_constants = lattice_constants
        self.lattice_matrix = lattice_matrix
        self.atom_types = atom_types
        self.atom_counts = atom_counts
        self.initial_positions = np.array(initial_positions)

    def set_slab_center(self, center):
        """Manually sets the center of the slab"""
        self.slab_center = np.array(center)

    def rotate_atoms(self, rotation_radius, rotation_angle_degrees):
        """Rotates atoms within a specified radius by a given angle in degrees"""
        rotation_angle_radians = np.radians(rotation_angle_degrees)
        rotation_matrix = np.array([[np.cos(rotation_angle_radians), -np.sin(rotation_angle_radians), 0],
                                    [np.sin(rotation_angle_radians), np.cos(rotation_angle_radians), 0],
                                    [0, 0, 1]])
        rotated_positions = []

        if self.slab_center is None:
            raise ValueError("Center of the slab is not set. Please use set_slab_center method to set the center manually.")

        for position in self.initial_positions:
            # Calculate the distance between the position and the center of the slab in the xy plane
            distance_to_center = np.linalg.norm(position[:2] - self.slab_center[:2])
            if distance_to_center <= rotation_radius:
                # Translate atom position to the center of the slab
                translated_position = position - self.slab_center
                # Rotate atom position
                rotated_position = np.dot(translated_position, rotation_matrix)
                # Translate back to the original position
                rotated_position += self.slab_center
                rotated_positions.append(rotated_position)
            else:
                rotated_positions.append(position)

        self.after_rotation = np.array(rotated_positions)

    def pair_atoms_with_type(self):
        """Pairs each XYZ position with its corresponding atom type"""
        self.atom_id = np.repeat(np.arange(len(self.atom_types)), self.atom_counts)

    def remove_overlapping_atoms(self, cutoff_distance=2.1):
        """Removes overlapping atoms based on a specified cutoff distance"""
        remaining_atoms = []
        remaining_atom_counts = np.zeros(len(self.atom_types), dtype=int)

        for i in range(len(self.after_rotation)):
            is_overlapping = False
            for j in range(i+1, len(self.after_rotation)):
                dist = np.linalg.norm(self.after_rotation[i] - self.after_rotation[j])
                if dist < cutoff_distance:
                    is_overlapping = True
                    break
            if not is_overlapping:
                remaining_atoms.append(self.after_rotation[i])
                remaining_atom_counts[self.atom_id[i]] += 1

        self.after_rotation = np.array(remaining_atoms)
        self.atom_counts = remaining_atom_counts

    def write_poscar(self, output_file):
        """Writes the modified POSCAR to a new file"""
        with open(output_file, 'w') as f:
            f.write(" ".join(self.atom_types) + "\n")
            for lattice_constant in self.lattice_constants:
                f.write(str(lattice_constant) + "\n")
            for row in self.lattice_matrix:
                f.write(" ".join(map(str, row)) + "\n")
            f.write(" ".join(map(str, self.atom_counts)) + "\n")
            f.write("Cartesian\n")
            for pos in self.after_rotation:
                f.write(" ".join(map(str, pos)) + "\n")

if __name__ == '__main__':
    rotation_radius = float(input("Enter the rotation radius: "))  # Input the rotation radius
    rotation_angle = float(input("Enter the rotation angle in degrees: "))  # Input the rotation angle in degrees

    poscar = Poscar("POSCAR1.vasp")

    # Manually set the center of the slab
    center_x = 54.6745
    center_y = 54.6745
    poscar.set_slab_center([center_x, center_y, 0])

    # Rotate atoms within the specified radius by the specified angle
    poscar.rotate_atoms(rotation_radius, rotation_angle)

    # Pair atoms with their types
    poscar.pair_atoms_with_type()

    # Remove overlapping atoms
    poscar.remove_overlapping_atoms()

    # Write the modified POSCAR to a new file
    poscar.write_poscar("rotated_and_filtered_POSCAR")
