// Copyright 2021 Global Phasing Ltd.
//
// Reindex merged or unmerged MTZ file.

#ifndef GEMMI_REINDEX_HPP_
#define GEMMI_REINDEX_HPP_

#include <ostream>  // for Mtz
#include "mtz.hpp"  // for Mtz

namespace gemmi {

// For now it's only partly-working
void reindex_mtz(Mtz& mtz, const Op& op, bool verbose, std::ostream* out) {
  mtz.switch_to_original_hkl();
  Op real_space_op{op.transposed_rot(), {0, 0, 0}};
  if (verbose && out)
    *out << "real space transformation: " << real_space_op.triplet() << '\n';
  size_t replace_row = size_t(-1);
  // change Miller indices
  for (size_t n = 0; n < mtz.data.size(); n += mtz.columns.size()) {
    Miller hkl_den = real_space_op.apply_to_hkl_without_division(mtz.get_hkl(n));
    Miller hkl = Op::divide_hkl_by_DEN(hkl_den);
    if (hkl[0] * Op::DEN == hkl_den[0] &&
        hkl[1] * Op::DEN == hkl_den[1] &&
        hkl[2] * Op::DEN == hkl_den[2]) {
      mtz.set_hkl(n, hkl);
    } else {
      if (replace_row == size_t(-1))
        replace_row = n;
      mtz.data[n] = NAN;  // mark for removal
    }
  }
  // remove reflections that came out fractional
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
  // hand change requires data modification
  if (op.det_rot() < 0)
    for (Mtz::Column& column : mtz.columns) {
      // negate anomalous difference
      if (column.type == 'D') {
        for (float& value : column)
          value = -value;
        if (verbose && out)
          *out << "Column " << column.label << ": anomalous difference negated.\n";
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
          if (verbose && out)
            *out << "Swapped columns " << column.label << " and " << minus_label << ".\n";
        } else {
          if (out)
            *out << "WARNING: matching pair not found for " << column.label << '\n';
        }
      }
    }
  // change space group
  mtz.spacegroup = find_spacegroup_by_change_of_basis(mtz.spacegroup, op.inverse());
  if (mtz.spacegroup) {
    if (verbose && out)
      *out << "Space group changed from " << mtz.spacegroup_name << " to "
           << mtz.spacegroup->xhm() << ".\n";
  } else {
    if (out)
      *out << "WARNING: new space group name could not be determined.\n";
  }
  // change unit cell parameters
  mtz.cell = mtz.cell.change_basis(real_space_op, false);
  for (Mtz::Dataset& ds : mtz.datasets)
    ds.cell = ds.cell.change_basis(real_space_op, false);
  for (Mtz::Batch& batch : mtz.batches)
    batch.set_cell(batch.get_cell().change_basis(real_space_op, false));

  if (mtz.is_merged())
    mtz.ensure_asu();
  else
    mtz.switch_to_asu_hkl();
}

} // namespace gemmi
#endif
