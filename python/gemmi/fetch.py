
import os
import urllib.request
import urllib.error
import time
from typing import Optional

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
    elif pdb_site == 'U': # also 'E', but _updated cif file
        url_prefix = "https://www.ebi.ac.uk/pdbe/entry-files/download"
        if use_cif:
            return f"{url_prefix}/{pdb_id.lower()}_updated.cif"
        else:
            raise RuntimeError(f"updated PDBe files must be in mmCIF not PDB format")
    elif pdb_site == 'J':
        filetype = ('M' if use_cif else 'P')
        return f"https://data.pdbj.org/pub/pdb/data/{gemmi.path_in_pdb_dir(pdb_id, filetype)}"
    else:
        raise RuntimeError(f"{pdb_site=}; expected 'R', 'E','U', or 'J' (RCSB/Europe/E-updated/Japan)")

# fetch and read a file for a PDB entry from RCSB/PDBe/PDBj/PDB_REDO
# if not available locally
def read_structure_with_code(pdb_id: str, use_cif: bool = True, pdb_site: str = 'R') -> Optional[gemmi.Structure]:
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
    return None

def main():
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Fetch PDB/CIF files.')
    parser.add_argument('-p', '--pdb', action='store_true',
                        help='fetch PDB format (default: mmCIF)')
    site_group = parser.add_mutually_exclusive_group()
    site_group.add_argument('-e', '--pdbe', action='store_const', const='E', dest='pdb_site',
                            help='fetch from PDBe')
    site_group.add_argument('-u', '--udpated', action='store_const', const='u', dest='pdb_site',
                            help='fetch "updated" file from PDBe')
    site_group.add_argument('-j', '--pdbj', action='store_const', const='J', dest='pdb_site',
                            help='fetch from PDBj')
    parser.set_defaults(pdb_site='R')
    parser.add_argument('pdb_ids', nargs='+', help='PDB IDs')
    args = parser.parse_args()

    use_cif = not args.pdb

    for pdb_id in args.pdb_ids:
        url = get_url_for_code(pdb_id, use_cif=use_cif, pdb_site=args.pdb_site)
        print(f'get {url}...')
        with urllib.request.urlopen(url, timeout=5) as response:
            content = response.read()
            ext = ('.cif' if use_cif else '.pdb')
            if args.pdb_site == 'u':
                ext = '_updated.cif'
            output_path = pdb_id + ext
            print(f'write {output_path}...')
            with open(output_path, 'wb') as output_file:
                output_file.write(content)

if __name__ == '__main__':
    main()
