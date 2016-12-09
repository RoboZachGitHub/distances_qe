import math
import os
import sys



class Atom:

	def __init__(self, element, cartesian_coordinates_list):
		self.element = element
		self.cartesian_coordinates = cartesian_coordinates_list

	def print_atom(self):
		x = self.cartesian_coordinates[0]
		y = self.cartesian_coordinates[1]
		z = self.cartesian_coordinates[2]
		print "{:2}   {: 8.4f}   {: 8.4f}   {: 8.4f}".format(self.element, x, y, z)

	def return_atom_string(self):
		x = self.cartesian_coordinates[0]
		y = self.cartesian_coordinates[1]
		z = self.cartesian_coordinates[2]
		return "{:2}   {: 8.4f}   {: 8.4f}   {: 8.4f}".format(self.element, x, y, z)

	def set_coordinates(self, cartesian_coordinates_list):
		self.cartesian_coordinates = cartesian_coordinates_list



def get_files():
	directory = os.getcwd()
	out_files = []
	for file in os.listdir(directory):
		if file.endswith(".out"):
			out_files.append(file)

	return out_files



def reasign_negative_coordinates(atom, primitive_vectors_list):
	x = atom.cartesian_coordinates[0]
	y = atom.cartesian_coordinates[1]
	z = atom.cartesian_coordinates[2]

	ax = primitive_vectors_list[0]
	ay = primitive_vectors_list[1]
	az = primitive_vectors_list[2]

	if x < 0:
		x = x + ax
	if y < 0:
		y = y + ay
	#if x < 0:
	#	z = z + az

	atom.set_coordinates([x, y, z])



def calculate_interatomic_distance(atom_1, atom_2):
	# due to the periodic boundary conditions, we must find the shorter distances, going across the box in positive or negative direction
	# in our case, only relevant for x and y coordinates
	x1 = atom_1.cartesian_coordinates[0]
	y1 = atom_1.cartesian_coordinates[1]
	z1 = atom_1.cartesian_coordinates[2]

	x2 = atom_2.cartesian_coordinates[0]
	y2 = atom_2.cartesian_coordinates[1]
	z2 = atom_2.cartesian_coordinates[2]

	del_x = math.fabs(x2-x1)
	del_y = math.fabs(y2-y1)
	del_z = math.fabs(z2-z1)

	del_x_tmp = math.fabs((a1-x2) - x1)	
	del_y_tmp = math.fabs((a2-y2) - y1)

	if del_x_tmp < del_x:
		del_x = del_x_tmp

	if del_y_tmp < del_y:
		del_y = del_y_tmp


	dx2 = math.pow(del_x, 2)
	dy2 = math.pow(del_y, 2)
	dz2 = math.pow(del_z, 2)
			
	interatomic_distance = math.sqrt(dx2 + dy2 + dz2)

	return interatomic_distance



def return_clean_data(in_lines, primitive_vectors_list):
	index_i = 0
	index_f = 0
	for i, line in reversed(list(enumerate(in_lines))):
		if line.find('ATOMIC_POSITIONS (angstrom)') != -1:
			index_i = i
			break	
	
	for i, line in reversed(list(enumerate(in_lines))):
		if line.find('End final coordinates') != -1:
			index_f = i
			break

	parsed_lines = input_lines[index_i + 1:index_f]

	atoms_list = []
	for line in parsed_lines:
		split_line = line.split()
		atom = Atom(split_line[0], [float(split_line[1]), float(split_line[2]), float(split_line[3])])
		atoms_list.append(atom)


	for atom in atoms_list:
		reasign_negative_coordinates(atom, primitive_vectors_list)

	return atoms_list


# main code starts here
out_files = get_files()

for out_file in out_files:

	#we will be adding to the output string throughout the routine
	output_string = ""

	input_name = out_file
	output_string += "Filename: {}\n\n".format(input_name)
	if input_name.find("110") != -1:
		tag = "110"
	elif input_name.find("100") != -1:
		tag = "100"
	else:
		tag = raw_input("supply tag for {}: ".format(input_name))

	if tag == "110":
		a1 = float(3.187)
		a2 = float(4.5071)
		a3 = float(4.5071)
	elif tag == "100":
		a1 = float(6.374)
		a2 = float(6.374)
		a3 = float(6.374)
	else:
		print tag
		print "incorrect tag"
		sys.exit()


	input_file = open(input_name, 'r')
	input_lines = input_file.readlines()
	input_file.close()

	list_of_atoms = return_clean_data(input_lines, [a1, a2, a3])

	# separate the atoms into two lists based on their element type, hydrogen or tungsten
	hydrogens_list = []
	tungstens_list = []

	for atom in list_of_atoms:
		if atom.element == "H":
			hydrogens_list.append(atom)
		elif atom.element == "W":
			tungstens_list.append(atom) 


	# find tungsten with largest z value
	highest_tungsten = float(0.0)
	for w in tungstens_list:
		z = w.cartesian_coordinates[2]
		if z > highest_tungsten:
			highest_tungsten = z

	output_string += "highest tungsten {:.3}\n".format(highest_tungsten)
	# sort hydrogens based on their z-value, for cleaner output	
	hydrogens_list = sorted(hydrogens_list, key=lambda x: x.cartesian_coordinates[2])
	
	# calculate the perpendicular distance to the surface for each h atom
	for i,h in enumerate(hydrogens_list):
	
		h_z = h.cartesian_coordinates[2]
		 
		perpendicular_distance = h_z - highest_tungsten
	

		# initialize at an unreasonably high value
		nearest_tungsten = float(10.0)
		for w in tungstens_list:
			interatomic_distance_tmp = calculate_interatomic_distance(w,h)	
			if interatomic_distance_tmp < nearest_tungsten:
				nearest_tungsten = interatomic_distance_tmp


		output_string += "\nhydrogen {}\n perpendicular distance: {:.3}\n nearest tungsten: {:.3}\n".format(i,perpendicular_distance, nearest_tungsten)
	
	
	output_file = open(input_name[:-3] + 'dist_info', 'w')
	output_file.write(output_string)
	output_file.close()
	input_file.close()

			
	



