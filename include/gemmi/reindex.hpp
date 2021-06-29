// Copyright 2021 Global Phasing Ltd.
//
// Reindex merged or unmerged MTZ file.

#ifndef GEMMI_REINDEX_HPP_
#define GEMMI_REINDEX_HPP_

#include "mtz.hpp"  // for Mtz

namespace gemmi {

// For now it's only partly-working
// TODO: handle operations that results in non-integral indices.
void reindex_mtz(Mtz& mtz, const Op& op, bool verbose) {
  mtz.switch_to_original_hkl();
  // transpose rotation matrix
  // TODO: handle op.rot that is not a rotation matrix (i.e. cell volume change)
  Op real_space_op{op.transposed_rot(), {0, 0, 0}};
  if (verbose)  { // for now we use mtz.warn for output
    if (!mtz.warnings)
      mtz.warnings = stderr;
    mtz.warn("real space transformation: " + real_space_op.triplet());
  }
  for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
    mtz.set_hkl(n, real_space_op.apply_to_hkl(mtz.get_hkl(n)));
  // hand change requires data modification
  if (op.det_rot() < 0)
    for (Mtz::Column& column : mtz.columns) {
      // negate anomalous difference
      if (column.type == 'D') {
        for (float& value : column)
          value = -value;
        if (verbose)
          mtz.warn("Column " + column.label + ": anomalous difference negated.");
        continue;
      }
      // swap (+) and (-)
      size_t pos = column.label.find("(+)");
      if (pos != std::string::npos) {
        std::string minus_label = column.label;
        minus_label[pos+1] = '-';
        Mtz::Column* minus_column =
            mtz.column_with_label(minus_label, &column.dataset());
        if (minus_column) {
          for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size())
            std::swap(mtz.data[n + column.idx],
                      mtz.data[n + minus_column->idx]);
          if (verbose)
            mtz.warn("Swapped values from " + column.label + " and " + minus_label + ".");
        } else {
          mtz.warn("Warning: matching pair not found for: " + column.label);
        }
      }
    }
  // change space group
  mtz.spacegroup = find_spacegroup_by_change_of_basis(mtz.spacegroup, op.inverse());
  if (mtz.spacegroup) {
    if (verbose)
      mtz.warn("Space group changed from " + mtz.spacegroup_name + " to " +
               mtz.spacegroup->xhm() + ".");
  } else {
    mtz.warn("WARNING: new space group name could not be determined.");
  }
  // change unit cell
  // TODO: change each unit cell separately, including batch headers
  UnitCell new_cell = mtz.cell.change_basis(real_space_op, false);
  mtz.set_cell_for_all(new_cell);

  if (mtz.is_merged())
    mtz.ensure_asu();
  else
    mtz.switch_to_asu_hkl();
}

} // namespace gemmi
#endif
