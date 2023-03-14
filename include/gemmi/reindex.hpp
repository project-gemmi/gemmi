// Copyright 2021 Global Phasing Ltd.
//
// Reindex merged or unmerged MTZ file.

#ifndef GEMMI_REINDEX_HPP_
#define GEMMI_REINDEX_HPP_

#include <ostream>
#include "mtz.hpp"  // for Mtz

namespace gemmi {

// For now it's only partly-working
inline void reindex_mtz(Mtz& mtz, const Op& op, std::ostream* out) {
  if (op.tran != Op::Tran{0, 0, 0})
    gemmi::fail("reindexing operator must not have a translation");
  mtz.switch_to_original_hkl();  // it changes hkl for unmerged data only
  Op transposed_op{op.transposed_rot(), {0, 0, 0}};
  Op real_space_op = transposed_op.inverse();
  if (out)
    *out << "Real space transformation: " << real_space_op.triplet() << '\n';
  size_t replace_row = size_t(-1);
  // change Miller indices
  for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size()) {
    Miller hkl_den = transposed_op.apply_to_hkl_without_division(mtz.get_hkl(n));
    Miller hkl = Op::divide_hkl_by_DEN(hkl_den);
    if (hkl[0] * Op::DEN == hkl_den[0] &&
        hkl[1] * Op::DEN == hkl_den[1] &&
        hkl[2] * Op::DEN == hkl_den[2]) {
      mtz.set_hkl(n, hkl);
    } else {  // fractional hkl - remove
      if (replace_row == size_t(-1))
        replace_row = n;
      mtz.data[n] = NAN;  // mark for removal
    }
  }

  // remove reflections marked for removal
  if (replace_row != size_t(-1)) {
    size_t ncol = mtz.columns.size();
    for (size_t n = replace_row + ncol; n < mtz.data.size(); n += ncol)
      if (!std::isnan(mtz.data[n])) {
        std::memcpy(mtz.data.data() + replace_row,
                    mtz.data.data() + n,
                    ncol * sizeof(float));
        replace_row += ncol;
      }
    if (out)
      *out << "Reflections removed (because of fractional indices): "
           << (mtz.data.size() - replace_row) / ncol << '\n';
    mtz.data.resize(replace_row);
    mtz.nreflections = int(replace_row / ncol);
  }

  // change space group
  const SpaceGroup* sg_before = mtz.spacegroup;
  if (sg_before) {
    GroupOps gops = sg_before->operations();
    gops.change_basis_impl(real_space_op, transposed_op);
    mtz.spacegroup = find_spacegroup_by_ops(gops);
  }
  if (mtz.spacegroup) {
    if (mtz.spacegroup != sg_before) {
      if (out)
        *out << "Space group changed from " << sg_before->xhm() << " to "
             << mtz.spacegroup->xhm() << ".\n";
      mtz.spacegroup_number = mtz.spacegroup->ccp4;
      mtz.spacegroup_name = mtz.spacegroup->hm;
    } else {
      if (out)
        *out << "Space group stays the same:" << sg_before->xhm() << ".\n";
    }
  } else {
    fail("reindexing: failed to determine new space group name");
  }

  // change unit cell parameters
  mtz.cell = mtz.cell.changed_basis_backward(transposed_op, false);
  for (Mtz::Dataset& ds : mtz.datasets)
    ds.cell = ds.cell.changed_basis_backward(transposed_op, false);
  for (Mtz::Batch& batch : mtz.batches)
    batch.set_cell(batch.get_cell().changed_basis_backward(transposed_op, false));

  // revert switch_to_original_hkl() for unmerged data
  mtz.switch_to_asu_hkl();
}

} // namespace gemmi
#endif
