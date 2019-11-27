from mercury_interface import MercuryInterface
from ccdc import descriptors

#Set up mercury interface.
helper = MercuryInterface() #assigning a name to the mercury interface
html_file = helper.output_html_file #creating a html file to display output on (in mercury)

#Create the output html file.
f = open(html_file, "w") # Opening the html file

#Read the input cif/res/entry.
entry = helper.current_entry
crystal = entry.crystal
molecule = entry.molecule
label = entry.identifier

original_molecules = crystal.molecule #Molecule before packing shell around it.
crystal_packed = crystal.packing_shell(packing_shell_size=50) #A generated shell of sym/translated molecules around centre cell.
potential_H_atoms = [] # an empty list to put all H atoms into
H_O_cutoff_distance = 2.5 #Cutoff for H...O distance
H_O_cutoff_angle = 10 #Angle From 180

for atom in crystal_packed.atoms: #for each atom in the generated shell (+ centre molecule):
    if atom.atomic_symbol == "H": # if it's a hydrogen atom
        potential_H_atoms.append(atom) #add it to the H atom list.

def H_Bonds(atom): #function to check if H bonded.
    for h_atom in potential_H_atoms: #check each H atom
        if h_atom not in atom.neighbours: #if hydrogen not connected to target atom (OH)
            R = descriptors.MolecularDescriptors.atom_distance(atom, h_atom) #Distance of O...H contact
            if R <= H_O_cutoff_distance: # if contact within defined cut-off
                h_O_atom = h_atom.neighbours[0] #O atom from OH component
                Angle = descriptors.MolecularDescriptors.atom_angle (atom, h_atom, h_O_atom) #OH....O angle
                if abs(180-Angle)<H_O_cutoff_angle: #if this is straight(ish)
                    return h_atom #return the h_atom in question.

total_O_atoms_in_component = 0 #start with 0 oxygen atoms
total_h_bonds_in_component = 0 #start with 0 h bonds.

for atom in original_molecules.atoms: #for each atom in the centre molecule
    if atom.atomic_symbol in ['O']: #if an oxygen (can be expanded... e.g. ['O','N']
        total_O_atoms_in_component = total_O_atoms_in_component +1 #add to the O count
        H_bond_for_atom = H_Bonds(atom) #check if it has a H bond
        if H_bond_for_atom: #if it does
            total_h_bonds_in_component = total_h_bonds_in_component + 1 #add 1 to count


O_h_sat = float(total_h_bonds_in_component)/float(total_O_atoms_in_component) #fraction of o with h bond

#report to user.
f.write("%d O atoms in structure, with %d O...H bonds (%s)" %(total_O_atoms_in_component,total_h_bonds_in_component,O_h_sat))