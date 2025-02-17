import ase
from ase.atoms import Atoms
import numpy as np
from ase.io import read, write

from rdkit import Chem
from rdkit.Chem import AllChem

from pysoftk.linear_polymer.linear_polymer import *
from pysoftk.format_printers.format_mol import *

def generate_and_export_polymer(pdb_file="mol.pdb", output_xyz="ptb7_pol10.xyz",
                                monomer_repeats=10, rel_steps=70000, rot_steps=45):
    """
    Generates a linear polymer from a PDB file, optimizes it, and exports it to an XYZ file.

    Args:
        pdb_file (str, optional): Path to the input PDB file. Defaults to "mol.pdb".
        output_xyz (str, optional): Path to the output XYZ file. Defaults to "ptb7_pol10.xyz".
        monomer_repeats (int, optional): Number of times the monomer is repeated. Defaults to 10.
        polymer_length (int, optional): Length of the polymer optimization. Defaults to 70000.
        polymer_steps (int, optional): Number of steps in the polymer optimization. Defaults to 45.

    Returns:
        bool: True if the process was successful, False otherwise.  Returns None if there was a file error.
    """
    try:
        monomer = Chem.MolFromPDBFile(pdb_file, removeHs=False)
        if monomer is None:  # Check if the PDB file was read successfully
            print(f"Error: Could not read PDB file '{pdb_file}'.")
            return None # Indicate file error
    except Exception as e:
        print(f"An error occurred while reading the PDB file: {e}")
        return None # Indicate file error

    AllChem.EmbedMolecule(monomer)  # Embedding is crucial

    try:
        polymer = Lp(monomer, "Br", monomer_repeats, shift=1.0).linear_polymer("MMFF", rel_steps, rot_steps, no_att=True)
        Fmt(polymer).xyz_print(output_xyz)
        return True

    except Exception as e:
        print(f"An error occurred during polymer generation or optimization: {e}")
        return False


import ase
from ase.atoms import Atoms
import numpy as np

def translate_and_rotate_to_sphere(atoms, R, theta, gamma):
    """
    Translates and rotates an ASE Atoms object to a sphere's surface.
    Handles the case where theta and gamma are zero correctly, translating by R.
    """

    if not isinstance(atoms, Atoms):
        raise TypeError("Input must be an ASE Atoms object.")

    if R <= 0:
        raise ValueError("Radius R must be positive.")

    new_atoms = atoms.copy()

    # Calculate translation vector (always, even if angles are zero)
    x = R * np.sin(theta) * np.cos(gamma)
    y = R * np.sin(theta) * np.sin(gamma)
    z = R * np.cos(theta)
    translation_vector = np.array([x, y, z])

    # Center before rotation (only if rotation is needed)
    if theta != 0 or gamma != 0:
        center_of_mass = new_atoms.get_center_of_mass()
        new_atoms.translate(-center_of_mass)

        # Rotate
        new_atoms.rotate(gamma, 'z')
        new_atoms.rotate(theta, 'y')

    # Translate to sphere surface (always)
    new_atoms.translate(translation_vector)

    return new_atoms

def process_and_combine_molecules(input_file="ptb7_pol10.xyz", output_file="both_molecules.xyz", R=35.0, theta=np.pi, gamma=np.pi/3):
    """
    Reads an XYZ file, translates and rotates the atoms to a sphere, 
    combines the original and translated atoms, and writes them to a new file.

    Args:
        input_file (str, optional): Path to the input XYZ file. Defaults to "ptb7_pol10.xyz".
        output_file (str, optional): Path to the output XYZ file. Defaults to "both_molecules.xyz".
        R (float, optional): Radius of the sphere. Defaults to 35.0.
        theta (float, optional): Polar angle for translation. Defaults to np.pi.
        gamma (float, optional): Azimuthal angle for translation. Defaults to np.pi / 3.

    Returns:
        bool: True if the process was successful, False otherwise.  Returns None if there was a file error.
    """
    try:
        original_atoms = read(input_file)
        atoms = original_atoms.copy()
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found.")
        return None  # Indicate file error
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return False # Indicate other error

    translated_atoms = translate_and_rotate_to_sphere(atoms, R, theta, gamma) 

    if translated_atoms:
        all_atoms = original_atoms.copy()
        all_atoms.extend(translated_atoms)

        write(output_file, all_atoms)
        print(f"Original and translated molecules written to '{output_file}'.")
        return True
    else:
        return False


if __name__ == "__main__":

    # Example with pysoftk and LP function for polymer creation.
    success_pt = generate_and_export_polymer(pdb_file="mol.pdb", output_xyz="ptb7_pol10.xyz",
                                             monomer_repeats=10, rel_steps=50000, rot_steps=45)
    if success_pt is None:
        print("File error occurred.")
    elif success_pt:
        print(f"Polymer generated and exported.")
    else:
        print("Polymer generation failed.")
    
    theta_degrees = 1.5  # Example: 90 degrees
    gamma_degrees = 1.5  # Example: 45 degrees

    theta_user = np.deg2rad(theta_degrees)
    gamma_user = np.deg2rad(gamma_degrees)
        
    success = process_and_combine_molecules(input_file="ptb7_pol10.xyz",
                                            output_file="combined_ptb7.xyz",
                                            R=20.0, theta=theta_user,
                                            gamma= gamma_user)
    if success is None:
        print("File error occurred.")
    elif success:
        print("Molecule processing complete.")
    else:
        print("Molecule processing failed.")
