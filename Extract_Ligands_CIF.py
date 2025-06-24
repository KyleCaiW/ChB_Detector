import os
import gzip
import csv
from tqdm import tqdm
import re
import sys

def process_cif_file(file_path):
    result = []
    pending_first_three = None

    try:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            content = f.readlines()
    except Exception as e:
        print(f"File reading error {file_path}: {e}", file=sys.stderr)
        return result

    start_index, end_index = None, None
    for i, line in enumerate(content):
        if line.startswith('_chem_comp.formula_weight'):
            start_index = i
        if start_index is not None and line.startswith('#'):
            end_index = i
            break

    if start_index is None or end_index is None:
        return result

    for line in content[start_index + 1 : end_index]:
        line = line.strip()

        if pending_first_three is not None:
            rest = line
            while True:
                match = re.search(r"'([^']+)'", rest)
                if not match:
                    break
                formula = match.group(1)
                if re.match(r'^([A-Z][a-z]?\s*[+-]?\d*(\s+[A-Z][a-z]?\s*[+-]?\d*)*)$', formula):
                    result.append((
                        os.path.basename(file_path).replace('.cif.gz', ''),
                        pending_first_three,
                        formula
                    ))
                    pending_first_three = None
                    break
                rest = rest[match.end():]

        if 'non-polymer' in line:
            parts = line.split('non-polymer', 1)
            pending_first_three = parts[0].strip()
            rest = parts[1] if len(parts) > 1 else ""
            while True:
                match = re.search(r"'([^']+)'", rest)
                if not match:
                    break
                formula = match.group(1)
                if re.match(r'^([A-Z][a-z]?\s*[+-]?\d*(\s+[A-Z][a-z]?\s*[+-]?\d*)*)$', formula):
                    result.append((
                        os.path.basename(file_path).replace('.cif.gz', ''),
                        pending_first_three,
                        formula
                    ))
                    pending_first_three = None
                    break
                rest = rest[match.end():]

    return result

def main(input_directory, output_file):
    file_names = [f for f in os.listdir(input_directory) if f.endswith('.cif.gz')]

    with open(output_file, 'w', newline='', encoding='utf-8') as out_file:
        writer = csv.writer(out_file)
        writer.writerow(['File Name', 'Name', 'Chemical Formula'])

        for file_name in tqdm(file_names, desc="Processing"):
            file_path = os.path.join(input_directory, file_name)
            try:
                matches = process_cif_file(file_path)
                if not matches:
                    continue

                for match in matches:
                    writer.writerow([
                        match[0],
                        match[1],
                        match[2]
                    ])

            except Exception as e:
                print(f"Error processing file {file_name}: {e}", file=sys.stderr)

if __name__ == "__main__":
    input_dir = "/data0/caiwenhao/PDB/PDBCIF/"
    output_path = "/data0/caiwenhao/PDB/cif_extract.csv"
    main(input_dir, output_path)