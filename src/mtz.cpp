// Copyright 2019-2023 Global Phasing Ltd.

#include <gemmi/mtz.hpp>
#include <gemmi/sprintf.hpp>

namespace gemmi {

#define WRITE(...) do { \
    int len = gstb_snprintf(buf, 81, __VA_ARGS__); \
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
    int len = gstb_snprintf(buffer, 18, "%.9f", f);
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
        gstb_snprintf(buf + pos, 7, "%6zu", j + 1);
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
