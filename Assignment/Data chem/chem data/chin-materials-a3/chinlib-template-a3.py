import argparse
import os, sys

from rdkit import Chem


smiles = ""


# Assign a custom ID to the atoms in a given molecule
def assign_custom_atom_id(mol):
    return # TODO


# Return SMILES primitive for a given bond
def get_bond_symbol(bond):
    return # TODO


# Return SMILES primitive for a given atom
def get_atom_symbol(atom):
    return # TODO



# Traverse molecular graph in depth-first order
def mol_dft(atom):
    global smiles

    return # TODO


def generate_smiles(mol):
    global smiles

    # Clear SMILES string
    smiles = ""

    # Get list of atoms
    atoms = mol.GetAtoms()

    return # TODO




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
    generate_smiles(mol)

    # Append SMILES to output file
    out_str = "{}\t{}".format(format(mol_id,'05d'), smiles)
    file_o.write(out_str + "\n")

    mol_id += 1

file_o.close()
