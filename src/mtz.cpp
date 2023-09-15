// Copyright 2019-2023 Global Phasing Ltd.

#include <gemmi/mtz.hpp>
#include <gemmi/gz.hpp>
#include <gemmi/sprintf.hpp>

namespace gemmi {

double wrap_degrees(double phi) {
  if (phi >= 0 && phi < 360.)
    return phi;
  return phi - std::floor(phi / 360.) * 360.;
}

void shift_phase(float& phi, double shift, bool negate=false) {
  double phi_ = phi + deg(shift);
  phi = float(wrap_degrees(negate ? -phi_ : phi_));
}

// apply phase shift to Hendrickson–Lattman coefficients HLA, HLB, HLC and HLD
void shift_hl_coefficients(float& a, float& b, float& c, float& d,
                           double shift, bool negate=false) {
  double sinx = std::sin(shift);
  double cosx = std::cos(shift);
  double sin2x = 2 * sinx * cosx;
  double cos2x = sq(cosx)- sq(sinx);
  // a sin(x+y) + b cos(x+y) = a sin(x) cos(y) - b sin(x) sin(y)
  //                         + a cos(x) sin(y) + b cos(x) cos(y)
  float a_ = float(a * cosx - b * sinx);
  float b_ = float(a * sinx + b * cosx);
  float c_ = float(c * cos2x - d * sin2x);
  float d_ = float(c * sin2x + d * cos2x);
  a = a_;                 // cos(phi)
  b = negate ? -b_ : b_;  // sin(phi)
  c = c_;                 // cos(2 phi)
  d = negate ? -d_ : d_;  // sin(2 phi)
}

// for probing/testing individual reflections, no need to optimize it
size_t Mtz::find_offset_of_hkl(const Miller& hkl, size_t start) const {
  if (!has_data() || columns.size() < 3)
    fail("No data.");
  if (start != 0)
    start -= (start % columns.size());
  for (size_t n = start; n + 2 < data.size(); n += columns.size())
    if (get_hkl(n) == hkl)
      return n;
  return (size_t)-1;
}

void Mtz::ensure_asu(bool tnt_asu) {
  if (!is_merged())
    fail("Mtz::ensure_asu() is for merged MTZ only");
  if (!spacegroup)
    return;
  GroupOps gops = spacegroup->operations();
  ReciprocalAsu asu(spacegroup, tnt_asu);
  std::vector<int> phase_columns = positions_of_columns_with_type('P');
  std::vector<int> abcd_columns = positions_of_columns_with_type('A');
  std::vector<int> dano_columns = positions_of_columns_with_type('D');
  std::vector<std::pair<int,int>> plus_minus_columns = positions_of_plus_minus_columns();
  bool no_special_columns = phase_columns.empty() && abcd_columns.empty() &&
                            plus_minus_columns.empty() && dano_columns.empty();
  bool centric = no_special_columns || gops.is_centrosymmetric();
  for (size_t n = 0; n < data.size(); n += columns.size()) {
    Miller hkl = get_hkl(n);
    if (asu.is_in(hkl))
      continue;
    auto result = asu.to_asu(hkl, gops);
    // cf. impl::move_to_asu() in asudata.hpp
    set_hkl(n, result.first);
    if (no_special_columns)
      continue;
    int isym = result.second;
    if (!phase_columns.empty() || !abcd_columns.empty()) {
      const Op& op = gops.sym_ops[(isym - 1) / 2];
      double shift = op.phase_shift(hkl);
      bool negate = (isym % 2 == 0);
      for (int col : phase_columns)
        shift_phase(data[n + col], shift, negate);
      for (auto i = abcd_columns.begin(); i+3 < abcd_columns.end(); i += 4)
        // we expect coefficients HLA, HLB, HLC and HLD - in this order
        shift_hl_coefficients(data[n + *(i+0)], data[n + *(i+1)],
                              data[n + *(i+2)], data[n + *(i+3)],
                              shift, negate);
    }
    if (isym % 2 == 0 && !centric &&
        // usually, centric reflections have empty F(-), so avoid swapping it
        !gops.is_reflection_centric(hkl)) {
      for (std::pair<int,int> cols : plus_minus_columns)
        std::swap(data[n + cols.first], data[n + cols.second]);
      for (int col : dano_columns)
        data[n + col] = -data[n + col];
    }
  }
}

void Mtz::reindex(const Op& op, std::ostream* out) {
  if (op.tran != Op::Tran{0, 0, 0})
    gemmi::fail("reindexing operator must not have a translation");
  if (op.det_rot() < 0)
    gemmi::fail("reindexing operator must preserve the hand of the axes");
  switch_to_original_hkl();  // changes hkl for unmerged data only
  Op transposed_op{op.transposed_rot(), {0, 0, 0}};
  Op real_space_op = transposed_op.inverse();
  if (out)
    *out << "Real space transformation: " << real_space_op.triplet() << '\n';
  bool row_removal = false;
  // change Miller indices
  for (size_t n = 0; n < data.size(); n += columns.size()) {
    Miller hkl_den = transposed_op.apply_to_hkl_without_division(get_hkl(n));
    Miller hkl = Op::divide_hkl_by_DEN(hkl_den);
    if (hkl[0] * Op::DEN == hkl_den[0] &&
        hkl[1] * Op::DEN == hkl_den[1] &&
        hkl[2] * Op::DEN == hkl_den[2]) {
      set_hkl(n, hkl);
    } else {  // fractional hkl - remove
      row_removal = true;
      data[n] = NAN;  // mark for removal
    }
  }

  // remove reflections marked for removal
  if (row_removal) {
    int n_before = nreflections;
    remove_rows_if([](float* h) { return std::isnan(*h); });
    if (out)
      *out << "Reflections removed (because of fractional indices): "
           << (n_before - nreflections) << '\n';
  }

  switch_to_asu_hkl();  // revert switch_to_original_hkl() for unmerged data

  // change space group
  if (spacegroup) {
    GroupOps gops = spacegroup->operations();
    gops.change_basis_impl(real_space_op, transposed_op);
    const SpaceGroup* new_sg = find_spacegroup_by_ops(gops);
    if (!new_sg)
      fail("reindexing: failed to determine new space group name");
    if (new_sg != spacegroup) {
      if (out)
        *out << "Space group changed from " << spacegroup->xhm() << " to "
             << new_sg->xhm() << ".\n";
      set_spacegroup(new_sg);
    } else {
      if (out)
        *out << "Space group stays the same:" << spacegroup->xhm() << ".\n";
    }
  }

  // change unit cell parameters
  cell = cell.changed_basis_backward(transposed_op, false);
  for (Mtz::Dataset& ds : datasets)
    ds.cell = ds.cell.changed_basis_backward(transposed_op, false);
  for (Mtz::Batch& batch : batches)
    batch.set_cell(batch.get_cell().changed_basis_backward(transposed_op, false));
}

void Mtz::expand_to_p1() {
  if (!spacegroup || !has_data())
    return;
  std::vector<int> phase_columns = positions_of_columns_with_type('P');
  std::vector<int> abcd_columns = positions_of_columns_with_type('A');
  bool has_phases = (!phase_columns.empty() || !abcd_columns.empty());
  GroupOps gops = spacegroup->operations();
  data.reserve(gops.sym_ops.size() * data.size());
  size_t orig_size = data.size();
  std::vector<Miller> hkl_copies;
  for (size_t n = 0; n < orig_size; n += columns.size()) {
    hkl_copies.clear();
    Miller hkl = get_hkl(n);
    // no reallocations because of reserve() above
    auto orig_iter = data.begin() + n;
    for (auto op = gops.sym_ops.begin() + 1; op < gops.sym_ops.end(); ++op) {
      Miller new_hkl = op->apply_to_hkl(hkl);
      Op::Miller negated{{-new_hkl[0], -new_hkl[1], -new_hkl[2]}};
      if (new_hkl != hkl && !in_vector(new_hkl, hkl_copies) &&
          negated != hkl && !in_vector(negated, hkl_copies)) {
        hkl_copies.push_back(new_hkl);
        size_t offset = data.size();
        data.insert(data.end(), orig_iter, orig_iter + columns.size());
        set_hkl(offset, new_hkl);
        if (has_phases) {
          double shift = op->phase_shift(hkl);
          if (shift != 0) {
            for (int col : phase_columns)
              shift_phase(data[offset + col], shift);
            for (auto i = abcd_columns.begin(); i+3 < abcd_columns.end(); i += 4)
              // we expect coefficients HLA, HLB, HLC and HLD - in this order
              shift_hl_coefficients(data[offset + *(i+0)], data[offset + *(i+1)],
                                    data[offset + *(i+2)], data[offset + *(i+3)], shift);
          }
        }
      }
    }
  }
  nreflections = int(data.size() / columns.size());
  sort_order = {{0, 0, 0, 0, 0}};
  set_spacegroup(&get_spacegroup_p1());
}

bool Mtz::switch_to_original_hkl() {
  if (indices_switched_to_original)
    return false;
  if (!has_data())
    fail("switch_to_original_hkl(): data not read yet");
  if (nreflections == 0) {
    // This function can be called before the data is populated
    // to set indices_switched_to_original, which is not exposed in Python.
    indices_switched_to_original = true;
    return true;
  }
  const Column* col = column_with_label("M/ISYM");
  if (col == nullptr || col->type != 'Y' || col->idx < 3)
    return false;
  std::vector<Op> inv_symops;
  inv_symops.reserve(symops.size());
  for (const Op& op : symops)
    inv_symops.push_back(op.inverse());
  for (size_t n = 0; n + col->idx < data.size(); n += columns.size()) {
    int isym = static_cast<int>(data[n + col->idx]) & 0xFF;
    const Op& op = inv_symops.at((isym - 1) / 2);
    Miller hkl = op.apply_to_hkl(get_hkl(n));
    int sign = (isym & 1) ? 1 : -1;
    for (int i = 0; i < 3; ++i)
      data[n+i] = static_cast<float>(sign * hkl[i]);
  }
  indices_switched_to_original = true;
  return true;
}

bool Mtz::switch_to_asu_hkl() {
  if (!indices_switched_to_original)
    return false;
  if (!has_data())
    fail("switch_to_asu_hkl(): data not read yet");
  const Column* col = column_with_label("M/ISYM");
  if (col == nullptr || col->type != 'Y' || col->idx < 3 || !spacegroup)
    return false;
  size_t misym_idx = col->idx;
  UnmergedHklMover hkl_mover(spacegroup);
  for (size_t n = 0; n + col->idx < data.size(); n += columns.size()) {
    Miller hkl = get_hkl(n);
    int isym = hkl_mover.move_to_asu(hkl);  // modifies hkl
    set_hkl(n, hkl);
    float& misym = data[n + misym_idx];
    misym = float(((int)misym & ~0xff) | isym);
  }
  indices_switched_to_original = false;
  return true;
}

void Mtz::read_file_gz(const std::string& path, bool with_data) {
  try {
    read_input(MaybeGzipped(path), with_data);
  } catch (std::runtime_error& e) {
    // append path to the error like in read_file(), but shouldn't the path go first?
    fail(std::string(e.what()) + ": " + path);
  }
}

#define WRITE(...) do { \
    int len = snprintf_z(buf, 81, __VA_ARGS__); \
    if (len < 80) \
      std::memset(buf + len, ' ', 80 - len); \
    if (write(buf, 80, 1) != 1) \
      sys_fail("Writing MTZ file failed"); \
  } while(0)

template<typename Write>
void Mtz::write_to_stream(Write write) const {
  // uses: data, spacegroup, nreflections, batches, cell, sort_order,
  //       valm, columns, datasets, history
  if (!has_data())
    fail("Cannot write Mtz which has no data");
  if (!spacegroup)
    fail("Cannot write Mtz which has no space group");
  char buf[81] = {'M', 'T', 'Z', ' ', '\0'};
  std::int64_t real_header_start = (int64_t) columns.size() * nreflections + 21;
  std::int32_t header_start = (int32_t) real_header_start;
  if (real_header_start > std::numeric_limits<int32_t>::max()) {
    header_start = -1;
  } else {
    real_header_start = 0;
  }
  std::memcpy(buf + 4, &header_start, 4);
  std::int32_t machst = is_little_endian() ? 0x00004144 : 0x11110000;
  std::memcpy(buf + 8, &machst, 4);
  std::memcpy(buf + 12, &real_header_start, 8);
  if (write(buf, 80, 1) != 1 ||
      write(data.data(), 4, data.size()) != data.size())
    fail("Writing MTZ file failed");
  WRITE("VERS MTZ:V1.1");
  WRITE("TITLE %s", title.c_str());
  WRITE("NCOL %8zu %12d %8zu", columns.size(), nreflections, batches.size());
  if (cell.is_crystal())
    WRITE("CELL  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f",
          cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma);
  WRITE("SORT  %3d %3d %3d %3d %3d", sort_order[0], sort_order[1],
        sort_order[2], sort_order[3], sort_order[4]);
  GroupOps ops = spacegroup->operations();
  char lat_type = spacegroup->ccp4_lattice_type();
  WRITE("SYMINF %3d %2d %c %5d %*s'%c%s' PG%s",
        ops.order(),               // number of symmetry operations
        (int) ops.sym_ops.size(),  // number of primitive operations
        lat_type,                  // lattice type
        spacegroup->ccp4,          // space group number
        20 - (int) std::strlen(spacegroup->hm), "",
        lat_type,                  // space group name (first letter)
        spacegroup->hm + 1,        // space group name (the rest)
        spacegroup->point_group_hm()); // point group name
  // If we have symops that are the same as spacegroup->operations(),
  // write symops to preserve the order of SYMM records.
  if (!symops.empty() && ops.is_same_as(split_centering_vectors(symops)))
    for (Op op : symops)
      WRITE("SYMM %s", to_upper(op.triplet()).c_str());
  else
    for (Op op : ops)
      WRITE("SYMM %s", to_upper(op.triplet()).c_str());
  auto reso = calculate_min_max_1_d2();
  WRITE("RESO %-20.12f %-20.12f", reso[0], reso[1]);
  if (std::isnan(valm))
    WRITE("VALM NAN");
  else
    WRITE("VALM %f", valm);
  auto format17 = [](float f) {
    char buffer[18];
    int len = snprintf_z(buffer, 18, "%.9f", f);
    return std::string(buffer, len > 0 ? std::min(len, 17) : 0);
  };
  for (const Column& col : columns) {
    auto minmax = calculate_min_max_disregarding_nans(col.begin(), col.end());
    const char* label = !col.label.empty() ? col.label.c_str() : "_";
    WRITE("COLUMN %-30s %c %17s %17s %4d",
          label, col.type,
          format17(minmax[0]).c_str(), format17(minmax[1]).c_str(),
          col.dataset_id);
    if (!col.source.empty())
      WRITE("COLSRC %-30s %-36s  %4d", label, col.source.c_str(), col.dataset_id);
  }
  WRITE("NDIF %8zu", datasets.size());
  for (const Dataset& ds : datasets) {
    WRITE("PROJECT %7d %s", ds.id, ds.project_name.c_str());
    WRITE("CRYSTAL %7d %s", ds.id, ds.crystal_name.c_str());
    WRITE("DATASET %7d %s", ds.id, ds.dataset_name.c_str());
    const UnitCell& uc = (ds.cell.is_crystal() && ds.cell.a > 0 ? ds.cell : cell);
    WRITE("DCELL %9d %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f",
          ds.id, uc.a, uc.b, uc.c, uc.alpha, uc.beta, uc.gamma);
    WRITE("DWAVEL %8d %10.5f", ds.id, ds.wavelength);
    for (size_t i = 0; i < batches.size(); i += 12) {
      std::memcpy(buf, "BATCH ", 6);
      int pos = 6;
      for (size_t j = i; j < std::min(batches.size(), i + 12); ++j, pos += 6)
        snprintf_z(buf + pos, 7, "%6zu", j + 1);
      std::memset(buf + pos, ' ', 80 - pos);
      if (write(buf, 80, 1) != 1)
        fail("Writing MTZ file failed");
    }
  }
  WRITE("END");
  if (!history.empty()) {
    // According to mtzformat.html the file can have only up to 30 history
    // lines, but we don't enforce it here.
    WRITE("MTZHIST %3zu", history.size());
    for (const std::string& line : history)
      WRITE("%s", line.c_str());
  }
  if (!batches.empty()) {
    WRITE("MTZBATS");
    for (const Batch& batch : batches) {
      // keep the numbers the same as in files written by libccp4
      WRITE("BH %8d %7zu %7zu %7zu",
            batch.number, batch.ints.size() + batch.floats.size(),
            batch.ints.size(), batch.floats.size());
      WRITE("TITLE %.70s", batch.title.c_str());
      if (batch.ints.size() != 29 || batch.floats.size() != 156)
        fail("wrong size of binaries batch headers");
      write(batch.ints.data(), 4, batch.ints.size());
      write(batch.floats.data(), 4, batch.floats.size());
      WRITE("BHCH  %7.7s %7.7s %7.7s",
            batch.axes.size() > 0 ? batch.axes[0].c_str() : "",
            batch.axes.size() > 1 ? batch.axes[1].c_str() : "",
            batch.axes.size() > 2 ? batch.axes[2].c_str() : "");
    }
  }
  WRITE("MTZENDOFHEADERS");
  if (!appended_text.empty()) {
    if (write(appended_text.data(), appended_text.size(), 1) != 1)
      fail("Writing MTZ file failed");
  }
}

#undef WRITE

void Mtz::write_to_cstream(std::FILE* stream) const {
  write_to_stream([&](const void *ptr, size_t size, size_t nmemb) {
      return std::fwrite(ptr, size, nmemb, stream);
  });
}

void Mtz::write_to_string(std::string& str) const {
  write_to_stream([&](const void *ptr, size_t size, size_t nmemb) {
      str.append(static_cast<const char*>(ptr), size * nmemb);
      return nmemb;
  });
}

void Mtz::write_to_file(const std::string& path) const {
  fileptr_t f = file_open(path.c_str(), "wb");
  try {
    return write_to_cstream(f.get());
  } catch (std::runtime_error& e) {
    fail(std::string(e.what()) + ": " + path);
  }
}

} // namespace gemmi
