import os
import gzip
import math
import csv
import numpy as np
from Bio.PDB import PDBParser
from tqdm import tqdm


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


def process_pdb(pdb_id, pdb_file, ligands, angle_range, s_centroid_distance, normal_angle_max):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)

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
            resname = residue.get_resname()
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
                            if other_residue.get_resname() in aromatic_residues:
                                ring_atoms = [atom for atom in other_residue.get_atoms() if atom.get_name() in aromatic_ring_atoms[other_residue.get_resname()]]
                                if ring_atoms:
                                    centroid = calculate_centroid(ring_atoms)
                                    centroid_residue = other_residue.get_resname() + str(other_residue.get_id()[1])
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
                                                        "sulfur_atom": sulfur_atom.get_name(),
                                                        "connected_atom_serial": connected_atom.get_serial_number(),
                                                        "connected_atom_name": connected_atom.get_name(),
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
                                                            "sulfur_atom": sulfur_atom.get_name(),
                                                            "connected_atom_serial": connected_atom.get_serial_number(),
                                                            "connected_atom_name": connected_atom.get_name(),
                                                            "connected_atom_coordinates": connected_atom.coord.tolist(),
                                                            "centroid_coordinates": centroid.tolist(),
                                                            "centroid_distance": centroid_distance,
                                                            "angle": angle,
                                                            "centroid_residue": centroid_residue,
                                                            "normal_angle": new_angle
                                                        })
    return (results_unrestricted, valid_ligands_count_unrestricted), (results_restricted, valid_ligands_count_restricted), found_ligands_count


def main():
    input_txt = "/data0/caiwenhao/PDB/sulfur_containing_ligands_2.txt"
    pdb_folder = "/data0/caiwenhao/PDB/PDB"
    output_csv_unrestricted = "/data0/caiwenhao/PDB/sulfurbond_aromatic_unrestricted.csv"
    output_csv_restricted = "/data0/caiwenhao/PDB/sulfurbond_aromatic_restricted.csv"

    min_angle = float(input("Enter minimum angle (°): "))
    max_angle = float(input("Enter maximum angle (°): "))
    s_centroid_distance = float(input("Enter maximum distance between sulfur atom and centroid (Å): "))
    normal_angle_max = float(input("Enter maximum normal angle (°): "))


    angle_range = (min_angle, max_angle)

    pdb_ligand_map = {}
    with open(input_txt, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            pdb_id = cols[0]
            ligands = cols[1::3]
            pdb_ligand_map[pdb_id] = ligands


    total_found_ligands = 0
    total_valid_ligands_unrestricted = 0
    total_valid_ligands_restricted = 0

    all_results_unrestricted = []
    all_results_restricted = []
    pdb_files = [os.path.join(pdb_folder, f"{pdb_id}.pdb.gz") for pdb_id in pdb_ligand_map.keys()]
    for pdb_file in tqdm(pdb_files, desc="Processing PDB files"):
        if not os.path.exists(pdb_file):
            print(f"File {pdb_file} does not exist, skipping...")
            continue

        pdb_id = os.path.basename(pdb_file).split('.')[0]
        with gzip.open(pdb_file, "rt") as f:
            (results_unrestricted, valid_ligands_count_unrestricted), (results_restricted, valid_ligands_count_restricted), found_ligands = process_pdb(pdb_id, f, pdb_ligand_map[pdb_id], angle_range, s_centroid_distance, normal_angle_max)
            all_results_unrestricted.extend(results_unrestricted)
            all_results_restricted.extend(results_restricted)
            total_found_ligands += found_ligands
            total_valid_ligands_unrestricted += valid_ligands_count_unrestricted
            total_valid_ligands_restricted += valid_ligands_count_restricted

    print(f"Total number of ligands found: {total_found_ligands}")
    print(f"Number of ligands meeting criteria (unrestricted angle): {total_valid_ligands_unrestricted}")
    print(f"Number of ligands meeting criteria (restricted angle): {total_valid_ligands_restricted}")

    with open(output_csv_unrestricted, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "PDB", "Ligand", "Sulfur Atom", "Connected Atom Serial",
            "Connected Atom Name", "Connected Atom Coordinates",
            "Centroid Coordinates", "Centroid Distance", "Angle", "Centroid Residue", "Normal Angle"
        ])
    
        current_pdb = None
        for row in all_results_unrestricted:
            if row["protein"] != current_pdb:
                current_pdb = row["protein"]
                writer.writerow([
                    row["protein"], row["ligand"], row["sulfur_atom"],
                    row["connected_atom_serial"], row["connected_atom_name"],
                    row["connected_atom_coordinates"],
                    row["centroid_coordinates"],
                    row["centroid_distance"], row["angle"], row["centroid_residue"], row["normal_angle"]
                ])
            else:
                writer.writerow([
                    "", row["ligand"], row["sulfur_atom"],
                    row["connected_atom_serial"], row["connected_atom_name"],
                    row["connected_atom_coordinates"],
                    row["centroid_coordinates"],
                    row["centroid_distance"], row["angle"], row["centroid_residue"], row["normal_angle"]
                ])

    with open(output_csv_restricted, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "PDB", "Ligand", "Sulfur Atom", "Connected Atom Serial",
            "Connected Atom Name", "Connected Atom Coordinates",
            "Centroid Coordinates", "Centroid Distance", "Angle", "Centroid Residue", "Normal Angle"
        ])
    
        current_pdb = None
        for row in all_results_restricted:
            if row["protein"] != current_pdb:
                current_pdb = row["protein"]
                writer.writerow([
                    row["protein"], row["ligand"], row["sulfur_atom"],
                    row["connected_atom_serial"], row["connected_atom_name"],
                    row["connected_atom_coordinates"],
                    row["centroid_coordinates"],
                    row["centroid_distance"], row["angle"], row["centroid_residue"], row["normal_angle"]
                ])
            else:
                writer.writerow([
                    "", row["ligand"], row["sulfur_atom"],
                    row["connected_atom_serial"], row["connected_atom_name"],
                    row["connected_atom_coordinates"],
                    row["centroid_coordinates"],
                    row["centroid_distance"], row["angle"], row["centroid_residue"], row["normal_angle"]
                ])


if __name__ == "__main__":
    main()