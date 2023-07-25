#!/home/sebastian/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt

def parse_atomic_positions(lines):
    atomic_positions = []
    found_cartesian = False

    for line in lines:
        if line.strip() == "Cartesian":
            found_cartesian = True
            continue

        if found_cartesian and line.strip() != "":
            try:
                coords = [float(coord) for coord in line.split()[:3]]
                atomic_positions.append(coords)
            except ValueError:
                continue

    if not found_cartesian:
        raise ValueError("Keyword 'Cartesian' not found in the input file.")

    return atomic_positions


def rearrange_by_z(atomic_positions):
    return sorted(atomic_positions, key=lambda pos: pos[2])


def write_updated_file(input_file, atomic_positions):
    output_filename = "inputfil.F"
    with open(input_file, "r") as infile, open(output_filename, "w") as outfile:
        for line in infile:
            outfile.write(line)
            if line.strip() == "Cartesian":
                break
        else:
            raise ValueError("Keyword 'Cartesian' not found in the input file.")

        for position in atomic_positions:
            outfile.write("{:15.10f} {:15.10f} {:15.10f}\n".format(*position))


def main():
    input_file = input("Enter the input file name: ")

    with open(input_file, "r") as infile:
        lines = infile.readlines()

    atomic_positions = parse_atomic_positions(lines)
    sorted_atomic_positions = rearrange_by_z(atomic_positions)
    write_updated_file(input_file, sorted_atomic_positions)

    print("Updated data has been saved to 'inputfil.F'.")


if __name__ == "__main__":
    main()

