// Test if all gemmi headers can be included after Windows headers.
// Windows headers contain typedefs and macros (#define small char)
// that cause surprising effects when these names are used in the code.
#ifdef _WIN32

// this is necessary, we use std::min etc.
#define NOMINMAX

// tinydir.h defines both UNICODE and _UNICODE if only one is defined.
// If including <windows.h> before tinydir.h, the same needs to be done.
#if defined(_UNICODE) && !defined(UNICODE)
# define UNICODE
#endif

#include <windows.h>
// To check for missing files run:
// for i in include/gemmi/*.hpp; do grep -qF `basename $i` tests/windowsh.cpp || echo $i; done
#include <gemmi/addends.hpp>
#include <gemmi/align.hpp>
#include <gemmi/assembly.hpp>
#include <gemmi/asudata.hpp>
#include <gemmi/asumask.hpp>
#include <gemmi/atof.hpp>
#include <gemmi/atox.hpp>
#include <gemmi/bessel.hpp>
#include <gemmi/binner.hpp>
#include <gemmi/blob.hpp>
#include <gemmi/bond_idx.hpp>
#include <gemmi/c4322.hpp>
#include <gemmi/calculate.hpp>
#include <gemmi/ccp4.hpp>
#include <gemmi/cellred.hpp>
#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>
#include <gemmi/cif2mtz.hpp>
#include <gemmi/cifdoc.hpp>
#include <gemmi/contact.hpp>
#include <gemmi/crd.hpp>
#include <gemmi/ddl.hpp>
#include <gemmi/dencalc.hpp>
#include <gemmi/dirwalk.hpp>
#include <gemmi/ecalc.hpp>
#include <gemmi/eig3.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/enumstr.hpp>
#include <gemmi/fail.hpp>
#include <gemmi/fileutil.hpp>
#include <gemmi/floodfill.hpp>
#include <gemmi/formfact.hpp>
#include <gemmi/fourier.hpp>
#include <gemmi/fprime.hpp>
#include <gemmi/fstream.hpp>
#include <gemmi/grid.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/input.hpp>
#include <gemmi/intensit.hpp>
#include <gemmi/interop.hpp>
#include <gemmi/it92.hpp>
#include <gemmi/iterator.hpp>
#include <gemmi/json.hpp>
#include <gemmi/levmar.hpp>
#include <gemmi/linkhunt.hpp>
#include <gemmi/logger.hpp>
#include <gemmi/math.hpp>
#include <gemmi/metadata.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/mmcif_impl.hpp>
//#include <gemmi/mmdb.hpp>
#include <gemmi/mmread.hpp>
#include <gemmi/mmread_gz.hpp>
#include <gemmi/model.hpp>
#include <gemmi/modify.hpp>
#include <gemmi/monlib.hpp>
#include <gemmi/mtz.hpp>
#include <gemmi/mtz2cif.hpp>
#include <gemmi/neighbor.hpp>
#include <gemmi/neutron92.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/pdb_id.hpp>
#include <gemmi/pirfasta.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/qcp.hpp>
#include <gemmi/read_cif.hpp>
#include <gemmi/recgrid.hpp>
#include <gemmi/reciproc.hpp>
#include <gemmi/refln.hpp>
#include <gemmi/resinfo.hpp>
#include <gemmi/riding_h.hpp>
#include <gemmi/scaling.hpp>
#include <gemmi/select.hpp>
#include <gemmi/seqalign.hpp>
#include <gemmi/seqid.hpp>
#include <gemmi/seqtools.hpp>
#include <gemmi/serialize.hpp>
#include <gemmi/sfcalc.hpp>
#include <gemmi/small.hpp>
#include <gemmi/smcif.hpp>
#include <gemmi/solmask.hpp>
#include <gemmi/span.hpp>
#include <gemmi/sprintf.hpp>
#include <gemmi/stats.hpp>
#include <gemmi/symmetry.hpp>
#include <gemmi/to_chemcomp.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/to_json.hpp>
#include <gemmi/to_mmcif.hpp>
#include <gemmi/to_pdb.hpp>
#include <gemmi/topo.hpp>
#include <gemmi/twin.hpp>
#include <gemmi/unitcell.hpp>
#include <gemmi/utf.hpp>
#include <gemmi/util.hpp>
#include <gemmi/version.hpp>
#include <gemmi/xds2mtz.hpp>
#include <gemmi/xds_ascii.hpp>
#endif
