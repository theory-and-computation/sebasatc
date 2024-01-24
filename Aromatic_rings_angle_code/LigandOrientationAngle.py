import numpy as np
import matplotlib.pyplot as plt
from sys import argv

# Written by Sebastian Casanova (1/18/2024)
# Edited by Liza Lee (1/23/2024)


#- Example Usage
#python ./LigandOrientationAngle.py 140 143
#- Input: XDATCAR
#- where 140 = index of the 1st Carbon in an aromatic ring
#      143 = index of the 2nd Carbon in an aromatic ring
#      Index starts from 0

class xdatcar:
    """ Python Class for VASP XDATCAR """

    def __init__(self, File=None):
        if File is None:
            self.xdatcar = 'XDATCAR'
        else:
            self.xdatcar = File

        self.TypeName = None
        self.ChemSymb = None
        self.Ntype = None
        self.Nions = None
        self.Nelem = None
        self.Niter = None

        # position in Direct Coordinate
        self.position = None
        # position in Cartesian Coordinate
        self.positionC = None
        self.readxdat()

    def readxdat(self):
        """ Read VASP XDATCAR """
        inp = [line for line in open(self.xdatcar) if line.strip()]
        scale = float(inp[1])

        self.cell = np.array([line.split() for line in inp[2:5]], dtype=float)
        self.cell *= scale

        ta = inp[5].split()
        tb = inp[6].split()

        if ta[0].isalpha():
            self.TypeName = ta
            self.Ntype = len(ta)
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
        else:
            print("VASP 4.X Format encountered...")
            self.Nelem = np.array(tb, dtype=int)
            self.Nions = self.Nelem.sum()
            self.Ntype = len(tb)
            self.TypeName = None

        pos = np.array([line.split() for line in inp[7:] if not line.split()[0].isalpha()], dtype=float)
        self.position = pos.ravel().reshape((-1, self.Nions, 3))
        self.Niter = self.position.shape[0]

        self.positionC = np.zeros_like(self.position)

        for i in range(self.Niter):
            self.positionC[i, :, :] = np.dot(self.position[i, :, :], self.cell)
            for j in range(0, 3):
                self.positionC[i, :, j] = np.mod(self.positionC[i, :, j], self.cell[j, j])

# Function to calculate the angle between two vectors
def calculate_angle(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    cos_theta = dot_product / (norm_v1 * norm_v2)
    angle = np.arccos(cos_theta) * 180 / np.pi  # Convert to degrees
    return angle

if __name__ == '__main__':
    inp = xdatcar()

    # Indices for carbon atoms and gold atom
    #carbon_atom_index1 = 140
    #carbon_atom_index2 = 143
    carbon_atom_index1 = int(argv[1])
    carbon_atom_index2 = int(argv[2])


    # Z-coordinate for the gold surface
    gold_surface_z = 16.0

    # Lists to store time and corresponding angles
    time_ps = []
    angles = []

    # Lists to store time and corresponding adjusted angles
    adjusted_angles = []

    # Assuming each configuration is separated by 1 ps
    time_step = 1.0

    # Loop through configurations and calculate angles
    for i in range(0, inp.Niter):
        config = inp.positionC[i, :, :]

        # Define vectors normal to planes formed by atoms
        plane_vector1 = config[carbon_atom_index2, :] - config[carbon_atom_index1, :]
        plane_vector2 = np.array([0.0, 0.0, gold_surface_z])  # Vector normal to the gold surface plane

        # Calculate the angle between planes
        angle = calculate_angle(plane_vector1, plane_vector2)

        # Adjust the angle by subtracting 90 degrees for reference angle=0 when parallel to the Au surface
        adjusted_angle = angle - 90.0

        time = i * time_step / 1000.0
        time_ps.append(time)
        angles.append(angle)
        adjusted_angles.append(adjusted_angle)

    """
    # Specify the file path where you want to save the output
    output_file_path = "output.txt"
    adjusted_output_file_path = "adjusted_output.txt"

    # Open the file in write mode
    with open(output_file_path, "w") as file:
        # Loop through configurations
        for i in range(inp.Niter):
            # Write configuration header to the file
            file.write(f"Configuration {i + 1}:\n")

            # Loop through atoms in the configuration
            for j in range(inp.Nions):
                # Write atom position to the file
                file.write(f"Atom {j + 1}: {config[j, :]}\n")

    # Print a message indicating that the data has been saved
   
    print(f"Adjusted data saved to {adjusted_output_file_path}")
    """

    # Plot the angle vs time graph
    plt.figure(figsize=(8, 6))

    #plt.subplot(2, 1, 1)
    plt.plot(time_ps, adjusted_angles, marker='o', markersize=1)
    plt.xlabel(r'$t$ (ps)', fontsize=20)
    plt.ylabel(r'$\theta$ (deg)',fontsize=20)
    plt.xlim((0, max(time_ps)))
    #plt.title('Angle between C-C axis and Gold surface plane vs Time')

    # Plot the adjusted angle vs time graph
    #plt.subplot(2, 1, 2)
    #plt.plot(time_fs, adjusted_angles, marker='o', markersize=3, color='orange')
    #plt.xlabel('Time (ps)')
    #plt.ylabel('Adjusted Angle (degrees)')
    #plt.title('Adjusted Angle between C-C plane and Gold surface plane vs Time')

    plt.tight_layout()
    plt.show()
