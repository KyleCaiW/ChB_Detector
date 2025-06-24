import os
import gzip
import math
import csv
import warnings
import numpy as np
from collections import defaultdict
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from tqdm import tqdm

warnings.simplefilter('ignore', PDBConstructionWarning)

def calculate_distance(atom1, atom2):
    if hasattr(atom2, 'coord'):
        coord2 = atom2.coord
    else:
        coord2 = atom2
    return np.linalg.norm(atom1.coord - coord2)

def calculate_angle(atom1, atom2, atom3):
    vector1 = np.array(atom1.coord) - np.array(atom2.coord)
    if hasattr(atom3, 'coord'):
        vector2 = np.array(atom3.coord) - np.array(atom2.coord)
    else:
        vector2 = atom3 - np.array(atom2.coord)
    norm1 = np.linalg.norm(vector1)
    norm2 = np.linalg.norm(vector2)
    cosine_angle = np.dot(vector1, vector2) / (norm1 * norm2)
    angle = math.acos(cosine_angle) * 180 / math.pi
    return angle

def find_connected_atoms(atom):
    connected_atoms = []
    for neighbor in atom.get_parent():
        if neighbor is not atom:
            distance = calculate_distance(atom, neighbor)
            if distance < 2.0:
                connected_atoms.append(neighbor)
    return connected_atoms

def is_divalent_sulfur(atom):
    connected_atoms = find_connected_atoms(atom)
    return len(connected_atoms) == 2

def calculate_centroid(atoms):
    coords = np.array([atom.coord for atom in atoms])
    centroid = np.mean(coords, axis=0)
    return centroid

def calculate_normal_vector(ring_atoms):
    if len(ring_atoms) >= 3:
        atom1 = np.array(ring_atoms[0].coord)
        atom2 = np.array(ring_atoms[1].coord)
        atom3 = np.array(ring_atoms[2].coord)
        vector1 = atom2 - atom1
        vector2 = atom3 - atom1
        normal_vector = np.cross(vector1, vector2)
        normal_vector /= np.linalg.norm(normal_vector)
        return normal_vector
    else:
        return None

def calculate_angle_with_normal(centroid, sulfur_atom, normal_vector):
    vector = np.array(sulfur_atom.coord) - np.array(centroid)
    cosine_angle = np.dot(vector, normal_vector) / np.linalg.norm(vector)
    angle = math.acos(abs(cosine_angle)) * 180 / math.pi
    return angle

def process_cif(pdb_id, cif_file, ligands, angle_range, s_centroid_distance, normal_angle_max):
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"Error parsing {pdb_id}: {e}")
        return ([], 0), ([], 0), 0

    results_unrestricted = []
    results_restricted = []
    processed_ligands = set()
    found_ligands_count = 0
    valid_ligands_count_unrestricted = 0
    valid_ligands_count_restricted = 0

    aromatic_residues = ["PHE", "TYR", "TRP", "HIS"]
    aromatic_ring_atoms = {
        "PHE": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TYR": ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
        "TRP": ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
        "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]
    }

    for chain in structure.get_chains():
        for residue in chain:
            resname = residue.get_resname().strip()
            if resname in ligands and resname not in processed_ligands:
                processed_ligands.add(resname)
                sulfur_atoms = [atom for atom in residue.get_atoms() if atom.element == "S"]

                found_ligands_count += 1
                for sulfur_atom in sulfur_atoms:
                    if not is_divalent_sulfur(sulfur_atom):
                        continue
                    connected_atoms = find_connected_atoms(sulfur_atom)

                    for other_chain in structure.get_chains():
                        for other_residue in other_chain:
                            if is_aa(other_residue) and other_residue.get_resname() in aromatic_residues:
                                ring_def = aromatic_ring_atoms[other_residue.get_resname()]
                                ring_atoms = [atom for atom in other_residue.get_atoms() if atom.name in ring_def]
                                
                                if ring_atoms:
                                    centroid = calculate_centroid(ring_atoms)
                                    centroid_residue = f"{other_residue.get_resname()}{other_residue.id[1]}"
                                    normal_vector = calculate_normal_vector(ring_atoms)
                                    
                                    if normal_vector is not None:
                                        for connected_atom in connected_atoms:
                                            centroid_distance = calculate_distance(sulfur_atom, centroid)
                                            if centroid_distance < s_centroid_distance:
                                                angle = calculate_angle(connected_atom, sulfur_atom, centroid)
                                                if angle_range[0] <= angle <= angle_range[1]:
                                                    new_angle = calculate_angle_with_normal(centroid, sulfur_atom, normal_vector)
                                                    
                                                    valid_ligands_count_unrestricted += 1
                                                    results_unrestricted.append({
                                                        "protein": pdb_id,
                                                        "ligand": resname,
                                                        "sulfur_atom": sulfur_atom.name,
                                                        "connected_atom_serial": connected_atom.serial_number,
                                                        "connected_atom_name": connected_atom.name,
                                                        "connected_atom_coordinates": connected_atom.coord.tolist(),
                                                        "centroid_coordinates": centroid.tolist(),
                                                        "centroid_distance": centroid_distance,
                                                        "angle": angle,
                                                        "centroid_residue": centroid_residue,
                                                        "normal_angle": new_angle
                                                    })
                                                    
                                                    if new_angle <= normal_angle_max:
                                                        valid_ligands_count_restricted += 1
                                                        results_restricted.append({
                                                            "protein": pdb_id,
                                                            "ligand": resname,
                                                            "sulfur_atom": sulfur_atom.name,
                                                            "connected_atom_serial": connected_atom.serial_number,
                                                            "connected_atom_name": connected_atom.name,
                                                            "connected_atom_coordinates": connected_atom.coord.tolist(),
                                                            "centroid_coordinates": centroid.tolist(),
                                                            "centroid_distance": centroid_distance,
                                                            "angle": angle,
                                                            "centroid_residue": centroid_residue,
                                                            "normal_angle": new_angle
                                                        })
    return (results_unrestricted, valid_ligands_count_unrestricted), (results_restricted, valid_ligands_count_restricted), found_ligands_count

def main():
    input_csv = "/data0/caiwenhao/PDB/cif_contain_S.csv"
    cif_folder = "/data0/caiwenhao/PDB/PDBCIF"
    output_csv_unrestricted = "/data0/caiwenhao/PDB/sulfurbond_aromatic_unrestricted_cif.csv"
    output_csv_restricted = "/data0/caiwenhao/PDB/sulfurbond_aromatic_restricted_cif.csv"

    min_angle = float(input("Enter minimum angle (°): "))
    max_angle = float(input("Enter maximum angle (°): "))
    s_centroid_distance = float(input("Enter maximum distance between sulfur atom and centroid (Å): "))
    normal_angle_max = float(input("Enter maximum normal angle (°): "))

    angle_range = (min_angle, max_angle)

    pdb_ligand_map = defaultdict(set)
    with open(input_csv, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            pdb_id = row[0].strip().upper()
            ligand = row[1].strip().upper()
            pdb_ligand_map[pdb_id].add(ligand)
    pdb_ligand_map = {k: list(v) for k, v in pdb_ligand_map.items()}

    cif_files = [os.path.join(cif_folder, f"{pdb_id}.cif.gz") for pdb_id in pdb_ligand_map.keys()]

    total_found_ligands = 0
    total_valid_ligands_unrestricted = 0
    total_valid_ligands_restricted = 0
    all_results_unrestricted = []
    all_results_restricted = []

    for cif_file in tqdm(cif_files, desc="Processing CIF files"):
        if not os.path.exists(cif_file):
            print(f"File {cif_file} does not exist, skipping...")
            continue
        
        pdb_id = os.path.basename(cif_file).split('.')[0].upper()
        try:
            with gzip.open(cif_file, "rt", errors='replace') as f:
                content = f.read().replace('\x00', '')
                with open("temp.cif", "w") as tmp:
                    tmp.write(content)
                
                with open("temp.cif", "r") as tmp:
                    (unrestricted, unrestricted_count), (restricted, restricted_count), found = process_cif(
                        pdb_id, tmp, 
                        pdb_ligand_map.get(pdb_id, []),
                        angle_range,
                        s_centroid_distance,
                        normal_angle_max
                    )
                
                all_results_unrestricted.extend(unrestricted)
                all_results_restricted.extend(restricted)
                total_found_ligands += found
                total_valid_ligands_unrestricted += unrestricted_count
                total_valid_ligands_restricted += restricted_count
        except Exception as e:
            print(f"An error occurred while processing file {cif_file}: {str(e)}")
        finally:
            if os.path.exists("temp.cif"):
                os.remove("temp.cif")

    print(f"\nStatistics:")
    print(f"Total number of ligands found: {total_found_ligands}")
    print(f"Number of ligands meeting criteria (unrestricted angle): {total_valid_ligands_unrestricted}")
    print(f"Number of ligands meeting criteria (restricted angle): {total_valid_ligands_restricted}")

    def write_output(filename, data):
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([
                "PDB", "Ligand", "Sulfur Atom", "Connected Atom Serial",
                "Connected Atom Name", "Connected Atom Coordinates",
                "Centroid Coordinates", "Centroid Distance", "Angle", 
                "Centroid Residue", "Normal Angle"
            ])
            current_pdb = None
            for row in data:
                if row["protein"] != current_pdb:
                    current_pdb = row["protein"]
                    writer.writerow([
                        row["protein"], row["ligand"], row["sulfur_atom"],
                        row["connected_atom_serial"], row["connected_atom_name"],
                        row["connected_atom_coordinates"],
                        row["centroid_coordinates"],
                        row["centroid_distance"], row["angle"], 
                        row["centroid_residue"], row["normal_angle"]
                    ])
                else:
                    writer.writerow([
                        "", row["ligand"], row["sulfur_atom"],
                        row["connected_atom_serial"], row["connected_atom_name"],
                        row["connected_atom_coordinates"],
                        row["centroid_coordinates"],
                        row["centroid_distance"], row["angle"], 
                        row["centroid_residue"], row["normal_angle"]
                    ])

    write_output(output_csv_unrestricted, all_results_unrestricted)
    write_output(output_csv_restricted, all_results_restricted)

if __name__ == "__main__":
    main()