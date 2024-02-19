#!/usr/bin/env python3

import sys
import zipfile

if len(sys.argv) < 3:
    sys.exit("specify whl files")

def get_names(zip_obj):
    return sorted(x for x in zip_obj.namelist() if not x.endswith('/'))

with zipfile.ZipFile(sys.argv[1], 'r') as z:
    print('Reference:', sys.argv[1])
    ref_lst = set(get_names(z))
for arg in sys.argv[2:]:
    print('Comparing with', arg)
    with zipfile.ZipFile(arg, 'r') as z:
        lst = get_names(z)
    for name in ref_lst:
        if name in lst:
            lst.remove(name)
        else:
            print('    - ', name)
    for name in lst:
        print('    + ', name)
