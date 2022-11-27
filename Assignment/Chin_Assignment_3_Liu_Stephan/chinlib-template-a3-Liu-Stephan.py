import argparse
import os, sys

from rdkit import Chem



smiles = ""


# Assign a custom ID to the atoms in a given molecule
def assign_custom_atom_id(mol):
    for i, atoms in enumerate(mol.GetAtoms()):
        atoms.SetIntProp("id", i)
        #print(atoms.GetIntProp("id"))




# Return SMILES primitive for a given bond
def get_bond_symbol(bond):
    for i in bond.GetBonds():
        #typing the output of GetBondType() which is 'BondType' to string to compare and return
        #primitive for SMILES
        if str(i.GetBondType()) == "SINGLE":
            return "-"
        elif str(i.GetBondType()) == "DOUBLE": 
            return "="
        elif str(i.GetBondType()) == "TRIPLE":
            return "#"
        
"""
        if i.GetBondType() == "SINGLE":
            print("-")
        elif i.GetBondType() == "DOUBLE": 
            print("=")
        elif i.GetBondType() == "TRIPLE":
            print("#")
"""

# Return SMILES primitive for a given atom
def get_atom_symbol(atom):
    for mol in atom.GetAtoms():
        #get atomic number and return corrsponding element symbol of that atomic number
        symbol = mol.GetAtomicNum()
        return Chem.GetPeriodicTable().GetElementSymbol(symbol)
    #print(Chem.rdchem.GetPeriodicTable().GetElementSymbol(atom.GetAtomicNum())) package is showns as unknown
    



# Traverse molecular graph in depth-first order
def mol_dft(atom):
    global smiles
    visited = set()
    for atoms in atom.GetAtoms():
    # Set to keep track of visited nodes of graph.
        if atoms.GetIdx() not in visited:
            #print(atoms.GetIdx())
            visited.add(atoms.GetIdx())
            #for neighbors in atoms.GetNeighbors():
                #mol_dft(neighbors)
    
#couldnt really implement the DFS so its just the basic shell of DFS so far
# the recursion hasnt been considered yet
    


def generate_smiles(mol):
    global smiles

    # Clear SMILES string
    smiles = ""

    # Get list of atoms
    atoms = mol.GetAtoms()
    for atom in atoms:
        symbol = atom.GetAtomicNum()
        smiles+= Chem.GetPeriodicTable().GetElementSymbol(symbol)
    
        for i in atom.GetBonds():
            if str(i.GetBondType()) == "SINGLE":
                smiles+= "-"
            elif str(i.GetBondType()) == "DOUBLE": 
                smiles+= "="
            elif str(i.GetBondType()) == "TRIPLE":
                smiles+= "#"
    print(smiles)
    #what is missing:
    # couldnt consider every bond index and keep neighbors in track, therefore more 
    # more singles and double bonds are appended to the smiles string
    # also single atoms were appeneded a single bond, 
    # fix: iterate through all atoms considering their actual position in the order
    # and also their previous and next atoms so the bond do not get appended redundantly


# ----------------------------------------------------------
# Main script
# ----------------------------------------------------------


# Command-line argument parsing

parser = argparse.ArgumentParser()
parser.add_argument("i", help="SDF MOL input file")
parser.add_argument("o", help="CSV output file")
parser.add_argument("--overwrite", help="If flag is specified, an already existing output file is overwritten", action="store_true")

args = parser.parse_args()

if not os.path.isfile(args.i):
    parser.print_help(sys.stderr)
    sys.exit(1)
if not args.i.endswith(".sdf"):
    parser.print_help(sys.stderr)
    sys.exit(1)
if os.path.isfile(args.o) and not args.overwrite:
    print("Output file already exists.")
    print("To overwrite use '--overwrite'.")
    sys.exit(1)


# Read SD input file

file_i = Chem.SDMolSupplier(args.i)  

# Open output file
file_o = open(args.o, "w")
file_o.write("mol_id\tSMILES\n")

# Iterate over molecules
mol_id = 1
for mol in file_i:
    print("-- Processing molecule " + str(mol_id), end='\r')

    # Generate SMILES
    assign_custom_atom_id(mol)
    get_atom_symbol(mol)
    get_bond_symbol(mol)
    generate_smiles(mol)

    # Append SMILES to output file
    out_str = "{}\t{}".format(format(mol_id,'05d'), smiles)
    file_o.write(out_str + "\n")

    mol_id += 1



file_o.close()
"""
if __name__ == "__main__":

    file_i = Chem.SDMolSupplier("smiles_01.sdf")

    mol_id = 1
    for mol in file_i:
        print("-- Processing molecule " + str(mol_id), end='\r')

        # Generate SMILES
        assign_custom_atom_id(mol)
        get_atom_symbol(mol)
        get_bond_symbol(mol)
        mol_dft(mol)
        generate_smiles(mol)
        
        # Append SMILES to output file
        mol_id += 1
"""   
