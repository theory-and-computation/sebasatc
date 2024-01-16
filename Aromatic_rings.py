import numpy as np
import matplotlib.pyplot as plt

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
    carbon_atom_index1 = 140
    carbon_atom_index2 = 143

    # Z-coordinate for the gold surface
    gold_surface_z = 16.0

    # Lists to store time and corresponding angles
    time_ps = []
    angles = []

    # Lists to store positions
    positions_list = []

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

        time = i * time_step
        time_ps.append(time)
        angles.append(angle)

        # Append positions to the positions_list
        positions_list.append(config)

    # Specify the file path where you want to save the output
    
    positions_file_path = "positions_list.txt"

    

    # Open the file in write mode for positions
    with open(positions_file_path, "w") as positions_file:
        # Loop through configurations
        for i, positions in enumerate(positions_list):
            # Write positions to the file
            positions_file.write(f"Configuration {i + 1} Positions:\n")
            for j, position in enumerate(positions):
                positions_file.write(f"Atom {j + 1}: {position}\n")

    # Print a message indicating that the data has been saved
   
    print(f"Positions list saved to {positions_file_path}")

    # Plot the angle vs time graph
    plt.plot(time_ps, angles, marker='o', markersize=3)
    plt.xlabel('Time (ps)')
    plt.ylabel('Angle (degrees)')
    plt.title('Angle between C-C plane and Gold surface plane vs Time')
    plt.show()
