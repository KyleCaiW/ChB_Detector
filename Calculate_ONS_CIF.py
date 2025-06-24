import os
import gzip
import math
import csv
import numpy as np
from collections import defaultdict
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from tqdm import tqdm


def calculate_distance(atom1, atom2):
    return np.linalg.norm(atom1.coord - atom2.coord)


def calculate_angle(atom1, atom2, atom3):
    vector1 = np.array(atom1.coord) - np.array(atom2.coord)
    vector2 = np.array(atom3.coord) - np.array(atom2.coord)
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


def process_cif(pdb_id, cif_file, ligands, distance_cutoffs, angle_range):
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(pdb_id, cif_file)
    except Exception as e:
        print(f"Error parsing {pdb_id}: {e}")
        return [], 0, 0

    results = []
    processed_ligands = set()
    found_ligands_count = 0
    valid_ligands_count = 0

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
                            if is_aa(other_residue, standard=True):
                                for atom in other_residue.get_atoms():
                                    if atom.element in ["O", "N", "S"]:
                                        distance = calculate_distance(sulfur_atom, atom)
                                        if distance < distance_cutoffs[atom.element]:
                                            for connected_atom in connected_atoms:
                                                angle = calculate_angle(connected_atom, sulfur_atom, atom)
                                                if angle_range[0] <= angle <= angle_range[1]:
                                                    valid_ligands_count += 1
                                                    results.append({
                                                        "protein": pdb_id,
                                                        "ligand": resname,
                                                        "sulfur_atom": sulfur_atom.get_name(),
                                                        "neighbor_atom_serial": atom.get_serial_number(),
                                                        "neighbor_atom_name": atom.get_name(),
                                                        "neighbor_residue": other_residue.get_resname(),
                                                        "neighbor_coordinates": atom.coord.tolist(),
                                                        "distance": distance,
                                                        "angle": angle,
                                                        "connected_atom_serial": connected_atom.get_serial_number(),
                                                        "connected_atom_name": connected_atom.get_name(),
                                                        "connected_atom_coordinates": connected_atom.coord.tolist()
                                                    })
    return results, found_ligands_count, valid_ligands_count


def main():
    input_csv = "/data0/caiwenhao/PDB/cif_contain_S.csv"
    cif_folder = "/data0/caiwenhao/PDB/PDBCIF/"
    output_csv = "/data0/caiwenhao/PDB/sulfurbond_ONS_cif.csv"

    s_o_distance = float(input("Enter max S-O distance (Å): "))
    s_n_distance = float(input("Enter max S-N distance (Å): "))
    s_s_distance = float(input("Enter max S-S distance (Å): "))
    min_angle = float(input("Enter minimum angle (°): "))
    max_angle = float(input("Enter maximum angle (°): "))

    distance_cutoffs = {"O": s_o_distance, "N": s_n_distance, "S": s_s_distance}
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
    total_valid_ligands = 0
    all_results = []

    for cif_file in tqdm(cif_files, desc="Processing CIF files"):
        if not os.path.exists(cif_file):
            print(f"File {cif_file} does not exist, skipping...")
            continue
        pdb_id = os.path.basename(cif_file).split('.')[0].upper()
        try:
            with gzip.open(cif_file, "rt") as f:
                results, found, valid = process_cif(pdb_id, f, pdb_ligand_map.get(pdb_id, []), distance_cutoffs, angle_range)
                all_results.extend(results)
                total_found_ligands += found
                total_valid_ligands += valid
        except Exception as e:
            print(f"An error occurred while processing file {cif_file}: {e}")
            continue

    print(f"Total number of ligands found: {total_found_ligands}")
    print(f"Number of ligands meeting the criteria: {total_valid_ligands}")

    with open(output_csv, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "PDB", "Ligand", "Sulfur Atom", "Neighbor Atom Serial",
            "Neighbor Atom Name", "Neighbor Residue",
            "Neighbor Coordinates", "Distance", "Angle",
            "Connected Atom Serial", "Connected Atom Name", "Connected Atom Coordinates"
        ])
        current_pdb = None
        for row in all_results:
            if row["protein"] != current_pdb:
                current_pdb = row["protein"]
                writer.writerow([
                    row["protein"], row["ligand"], row["sulfur_atom"],
                    row["neighbor_atom_serial"], row["neighbor_atom_name"],
                    row["neighbor_residue"], row["neighbor_coordinates"],
                    row["distance"], row["angle"],
                    row["connected_atom_serial"], row["connected_atom_name"],
                    row["connected_atom_coordinates"]
                ])
            else:
                writer.writerow([
                    "", row["ligand"], row["sulfur_atom"],
                    row["neighbor_atom_serial"], row["neighbor_atom_name"],
                    row["neighbor_residue"], row["neighbor_coordinates"],
                    row["distance"], row["angle"],
                    row["connected_atom_serial"], row["connected_atom_name"],
                    row["connected_atom_coordinates"]
                ])


if __name__ == "__main__":
    main()