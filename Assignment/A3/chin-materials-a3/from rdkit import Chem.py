from rdkit import Chem



def get_atoms():
    m = Chem.MolFromSmiles("C1OC1")
    for atom in m.GetAtoms():
        print(atom.GetAtomicNum())
    
    for b in m:
        print(m.GetBonds()[b].GetBondType())
    

if __name__ == "__main__":
    get_atoms()