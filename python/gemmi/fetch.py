
import os
import urllib.request
import urllib.error
import time

import gemmi

def get_url_for_code(pdb_id: str, use_cif: bool = True, pdb_site: str = 'R') -> str:
    pdb_site = pdb_site.upper()[0]
    if pdb_site == 'R':
        url_prefix = f"https://files.rcsb.org/download/{pdb_id.upper()}"
        if use_cif:
            return f"{url_prefix}.cif"
        else:
            return f"{url_prefix}.pdb"
    elif pdb_site == 'E':
        url_prefix = "https://www.ebi.ac.uk/pdbe/entry-files/download"
        if use_cif:
            return f"{url_prefix}/{pdb_id.lower()}.cif"
        else:
            return f"{url_prefix}/pdb{pdb_id.lower()}.ent"
    elif pdb_site == 'J':
        filetype = ('M' if use_cif else 'P')
        return f"https://data.pdbj.org/pub/pdb/data/{gemmi.path_in_pdb_dir(pdb_id, filetype)}"
    else:
        raise RuntimeError(f"{pdb_site=}; expected 'R', 'E', or 'J' (RCSB/Europe/Japan)")

# fetch and read a file for a PDB entry from RCSB/PDBe/PDBj/PDB_REDO
# if not available locally
def read_structure_with_code(pdb_id: str, use_cif: bool = True, pdb_site: str = 'R') -> gemmi.Structure:
    assert gemmi.is_pdb_code(pdb_id)
    local_path = gemmi.expand_if_pdb_code(pdb_id[0])
    if os.path.isfile(local_path):
        return gemmi.read_structure(local_path)
    url = get_url_for_code(pdb_id, use_cif, pdb_site)
    for i in range(3):
        try:
            with urllib.request.urlopen(url, timeout=10) as response:
                content = response.read()
                if use_cif:
                    doc = gemmi.cif.read_string(content)
                    return gemmi.make_structure_from_block(doc[0])
                else:
                    return gemmi.read_pdb_string(content)
        except urllib.error.URLError as e:
            if i == 3:
                raise
            print(f"Retrying {url} after error: {e}")
            time.sleep(2 * (i + 1))

if __name__ == '__main__':
    import sys
    use_cif = True
    for pdb_id in sys.argv[1:]:
        url = get_url_for_code(pdb_id, use_cif=use_cif)
        print(f'get {url}...')
        with urllib.request.urlopen(url, timeout=10) as response:
            content = response.read()
            output_path = pdb_id + ('.cif' if use_cif else '.pdb')
            print(f'write {output_path}...')
            with open(output_path, 'wb') as output_file:
                output_file.write(content)
