#####Original Code (with original calculate_distance method):###
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

    # Function to calculate the distance between two atoms
    def calculate_distance(self, atom1, atom2):
        delta = atom1 - atom2
        for i in range(3):
            if delta[i] > self.cell[i, i] / 2:
                delta[i] -= self.cell[i, i]
            elif delta[i] < -self.cell[i, i] / 2:
                delta[i] += self.cell[i, i]
        return np.linalg.norm(delta)


if __name__ == '__main__':
    inp1 = xdatcar("XDATCAR_final_1") # Name of first XDATCAR
    inp2 = xdatcar("XDATCAR_final_2")  # Naem of the Second 

    # Indices for atoms C* and N
    atom_index1_1 = 35
    atom_index2_1 = 94

    # Indices for atoms C* and O
    atom_index1_2 = 32
    atom_index2_2 = 35

    # Lists to store distances
    
    # C* and N from Reactants
    distances_1_1 = []
    
    # C* and O from Reactants
    distances_2_1 = []
    
    #  C* and N from Products
    distances_1_2 = []
    
     # C* and O from Products
    distances_2_2 = []

    # Loop through configurations and calculate distances for XDATCAR_final_1
    for i in range(inp1.Niter):
        config1 = inp1.positionC[i, :, :]
        # Calculate distance between atoms in XDATCAR_final_1
        dist1_1 = inp1.calculate_distance(config1[atom_index1_1], config1[atom_index2_1])
        dist2_1 = inp1.calculate_distance(config1[atom_index1_2], config1[atom_index2_2])
        distances_1_1.append(dist1_1)
        distances_2_1.append(dist2_1)

    # Loop through configurations and calculate distances for XDATCAR_final_2
    for i in range(inp2.Niter):
        config2 = inp2.positionC[i, :, :]
        # Calculate distance between atoms in XDATCAR_final_2
        dist1_2 = inp2.calculate_distance(config2[atom_index1_1], config2[atom_index2_1])
        dist2_2 = inp2.calculate_distance(config2[atom_index1_2], config2[atom_index2_2])
        distances_1_2.append(dist1_2)
        distances_2_2.append(dist2_2)

    # Plot distances from both XDATCAR files
    plt.scatter(distances_2_1, distances_1_1, color='blue', marker='o', label='XDATCAR_final_1', s=20)
    plt.scatter(distances_2_2, distances_1_2, color='red', marker='x', label='XDATCAR_final_2', s=20)

    # Adding labels and title
    plt.xlabel('d C*-O')
    plt.ylabel('d C*-X*')
    plt.title('Comparison of Collective Variables')
    #plt.legend()

    # Displaying the plot
    plt.show()

    # Save distances to files
    np.savetxt('distances_1_1.txt', distances_1_1)
    np.savetxt('distances_2_1.txt', distances_2_1)
    np.savetxt('distances_1_2.txt', distances_1_2)  ###note should be the same since N and C* should be bonded
    np.savetxt('distances_2_2.txt', distances_2_2)

#####Modified Code (with calculate_distance_xy method):###
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

    # Function to calculate the distance between two atoms considering only x and y coordinates
    def calculate_distance_xy(self, atom1, atom2):
        delta_xy = atom1[:2] - atom2[:2]  # Consider only x and y coordinates
        for i in range(2):  # Iterate over x and y coordinates
            if delta_xy[i] > self.cell[i, i] / 2:
                delta_xy[i] -= self.cell[i, i]
            elif delta_xy[i] < -self.cell[i, i] / 2:
                delta_xy[i] += self.cell[i, i]
        return np.linalg.norm(delta_xy)


if __name__ == '__main__':
    inp1 = xdatcar("XDATCAR_final_1")
    inp2 = xdatcar("XDATCAR_final_2")

    # Indices for atoms C* and N 
    atom_index1_1 = 35
    atom_index2_1 = 94

    # Indices for atoms C* and O
    atom_index1_2 = 32
    atom_index2_2 = 35

    # Lists to store distances
    
    # C* and N from R
    distances_1_1 = []
    
    # C* and O from R
    distances_2_1 = []
    
    #  C* and N from P
    distances_1_2 = []
    
     # C* and O from P
    distances_2_2 = []

    # Loop through configurations and calculate distances for XDATCAR_final_1
    for i in range(inp1.Niter):
        config1 = inp1.positionC[i, :, :]
        # Calculate distance between atoms in XDATCAR_final_1 considering only x and y coordinates
        dist1_1 = inp1.calculate_distance_xy(config1[atom_index1_1], config1[atom_index2_1])
        dist2_1 = inp1.calculate_distance_xy(config1[atom_index1_2], config1[atom_index2_2])
        distances_1_1.append(dist1_1)
        distances_2_1.append(dist2_1)

    # Loop through configurations and calculate distances for XDATCAR_final_2
    for i in range(inp2.Niter):
        config2 = inp2.positionC[i, :, :]
        # Calculate distance between atoms in XDATCAR_final_2 considering only x and y coordinates
        dist1_2 = inp2.calculate_distance_xy(config2[atom_index1_1], config2[atom_index2_1])
        dist2_2 = inp2.calculate_distance_xy(config2[atom_index1_2], config2[atom_index2_2])
        distances_1_2.append(dist1_2)
        distances_2_2.append(dist2_2)

    # Plot distances from both XDATCAR files
    plt.scatter(distances_2_1, distances_1_1, color='blue', marker='o', label='XDATCAR_final_1', s=20)
    plt.scatter(distances_2_2, distances_1_2, color='red', marker='x', label='XDATCAR_final_2', s=20)

    # Adding labels and title
    plt.xlabel('d C*-O (xy)')
    plt.ylabel('d C*-X* (xy)')
    plt.title('Comparison of Collective Variables considering only xy-coordinates')
    #plt.legend()

    # Displaying the plot
    plt.show()

    
