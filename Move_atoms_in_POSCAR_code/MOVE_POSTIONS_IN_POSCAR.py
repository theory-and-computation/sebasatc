#!/Users/ieatzombies/mambaforge3/bin/python
def update_positions_in_file(file_name, start_line, end_line, dx, dy, dz):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    for i in range(start_line - 1, end_line):
        parts = lines[i].split()
        if len(parts) >= 3:
            x = round(float(parts[0]) + dx, 9)
            y = round(float(parts[1]) + dy, 9)
            z = round(float(parts[2]) + dz, 9)
            lines[i] = f"{x} {y} {z}\n"
    with open(file_name, 'w') as f:
        f.writelines(lines)
def main():
    file_name = input("Enter the name of the POSCAR file: ")
    start_line = int(input("Enter the start line: "))
    end_line = int(input("Enter the end line: "))
    dx = float(input("Enter the x shift: "))
    dy = float(input("Enter the y shift: "))
    dz = float(input("Enter the z shift: "))
    update_positions_in_file(file_name, start_line, end_line, dx, dy, dz)
    print(f"Updated positions in '{file_name}'.")
if __name__ == "__main__":
    main()
