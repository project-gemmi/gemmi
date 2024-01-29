// Copyright Global Phasing Ltd.

#include <gemmi/intensit.hpp>
#include <gemmi/atof.hpp>  // for fast_from_chars

namespace gemmi {

namespace {

bool parse_voigt_notation(const char* start, const char* end, SMat33<double>& b) {
  bool first = true;
  for (double* u : {&b.u11, &b.u22, &b.u33, &b.u23, &b.u13, &b.u12}) {
    if (*start != (first ? '(' : ','))
      return false;
    first = false;
    auto result = fast_from_chars(++start, end, *u);
    if (result.ec != std::errc())
      return false;
    start = skip_blank(result.ptr);
  }
  return *start == ')';
}

template<typename Source>
void copy_metadata(const Source& source, Intensities& intensities) {
  intensities.unit_cell = source.cell;
  intensities.spacegroup = source.spacegroup;
  if (!intensities.spacegroup)
    fail("unknown space group");
}

inline void add_if_valid(Intensities& intensities,
                         const Miller& hkl, short isign, double value, double sigma) {
  // XDS marks rejected reflections with negative sigma.
  // Sigma 0.0 is also problematic - it rarely happens (e.g. 5tkn).
  if (!std::isnan(value) && sigma > 0)
    intensities.data.push_back({hkl, isign, /*nobs=*/0, value, sigma});
}

template<typename DataProxy>
void read_data(Intensities& intensities, const DataProxy& proxy,
               size_t value_idx, size_t sigma_idx) {
  for (size_t i = 0; i < proxy.size(); i += proxy.stride())
    add_if_valid(intensities, proxy.get_hkl(i), 0,
                 proxy.get_num(i + value_idx),
                 proxy.get_num(i + sigma_idx));
}

template<typename DataProxy>
void read_anomalous_data(Intensities& intensities, const DataProxy& proxy,
                         int mean_idx, size_t (&value_idx)[2], size_t (&sigma_idx)[2]) {
  GroupOps gops = intensities.spacegroup->operations();
  for (size_t i = 0; i < proxy.size(); i += proxy.stride()) {
    Miller hkl = proxy.get_hkl(i);
    bool centric = gops.is_reflection_centric(hkl);

    // sanity check
    if (mean_idx >= 0 && !std::isnan(proxy.get_num(i + mean_idx)) && !centric) {
      if (std::isnan(proxy.get_num(i + value_idx[0])) &&
          std::isnan(proxy.get_num(i + value_idx[1])))
        fail(miller_str(hkl), " has <I>, but I(+) and I(-) are both null");
    }

    add_if_valid(intensities, hkl, 1, proxy.get_num(i + value_idx[0]),
                                      proxy.get_num(i + sigma_idx[0]));
    if (!centric)  // ignore I(-) of centric reflections
      add_if_valid(intensities, hkl, -1, proxy.get_num(i + value_idx[1]),
                                         proxy.get_num(i + sigma_idx[1]));
  }
}

bool read_staraniso_b_from_mmcif(const cif::Block& block, SMat33<double>& output) {
  // read what is written by MtzToCif::write_staraniso_b()
  cif::Table table = const_cast<cif::Block*>(&block)->find(
      "_reflns.pdbx_aniso_B_tensor_eigen",
      {"value_1", "value_2", "value_3",
       "vector_1_ortho[1]", "vector_1_ortho[2]", "vector_1_ortho[3]",
       "vector_2_ortho[1]", "vector_2_ortho[2]", "vector_2_ortho[3]",
       "vector_3_ortho[1]", "vector_3_ortho[2]", "vector_3_ortho[3]"});
  if (!table.ok())
    return false;
  cif::Table::Row row = table.one();
  using cif::as_number;
  double eigval[3] = {as_number(row[0]), as_number(row[1]), as_number(row[2])};
  double min_val = std::min(std::min(eigval[0], eigval[1]), eigval[2]);
  Mat33 mat(as_number(row[3]), as_number(row[6]), as_number(row[9]),
            as_number(row[4]), as_number(row[7]), as_number(row[10]),
            as_number(row[5]), as_number(row[8]), as_number(row[11]));
  Vec3 diag(eigval[0] - min_val, eigval[1] - min_val, eigval[2] - min_val);
  // If the columns of mat are an orthonomal basis, mat^−1==mat^T.
  // But just in case it isn't, we use mat^-1 here.
  Mat33 t = mat.multiply_by_diagonal(diag).multiply(mat.inverse());
  // t is a symmetric tensor, so we return only 6 numbers
  output = {t[0][0], t[1][1], t[2][2], t[0][1], t[0][2], t[1][2]};
  return true;
}

// returns STARANISO version or empty string
std::string read_staraniso_b_from_mtz(const Mtz& mtz, SMat33<double>& output) {
  std::string version;
  size_t hlen = mtz.history.size();
  for (size_t i = 0; i != hlen; ++i)
    if (mtz.history[i].find("STARANISO") != std::string::npos) {
      size_t version_pos = mtz.history[i].find("version:");
      if (version_pos != std::string::npos)
        version = read_word(mtz.history[i].c_str() + version_pos + 8);
      else
        version = "?";
      // StarAniso 2.3.74 (24-Apr-2021) and later write B tensor in history
      for (size_t j = i+1; j < std::min(i+4, hlen); ++j) {
        const std::string& line = mtz.history[j];
        if (starts_with(line, "B=(")) {
          if (!parse_voigt_notation(line.c_str() + 2,
                                    line.c_str() + line.size(),
                                    output))
            fail("failed to parse tensor Voigt notation: " + line);
          break;
        }
      }
      break;
    }
  return version;
}

} // anonymous namespace

std::array<double,2> Intensities::resolution_range() const {
  double min_1_d2 = INFINITY;
  double max_1_d2 = 0;
  for (const Refl& x : data) {
    double a_1_d2 = unit_cell.calculate_1_d2(x.hkl);
    if (a_1_d2 < min_1_d2)
      min_1_d2 = a_1_d2;
    if (a_1_d2 > max_1_d2)
      max_1_d2 = a_1_d2;
  }
  return {{ 1 / std::sqrt(min_1_d2), 1 / std::sqrt(max_1_d2) }};
}

// cf. calculate_hkl_value_correlation()
Correlation Intensities::calculate_correlation(const Intensities& other) const {
  Correlation corr;
  auto r1 = data.begin();
  auto r2 = other.data.begin();
  while (r1 != data.end() && r2 != other.data.end()) {
    if (r1->hkl == r2->hkl && r1->isign == r2->isign) {
      corr.add_point(r1->value, r2->value);
      ++r1;
      ++r2;
    } else if (*r1 < *r2) {
      ++r1;
    } else {
      ++r2;
    }
  }
  return corr;
}

void Intensities::merge_in_place(DataType data_type) {
  type = data_type;
  if (data.empty())
    return;
  if (data_type == DataType::Mean)
    // discard signs so that merging produces Imean
    for (Refl& refl : data)
      refl.isign = 0;
  sort();
  std::vector<Refl>::iterator out = data.begin();
  double sum_wI = 0.;
  double sum_w = 0.;
  int nobs = 0;
  for (auto in = data.begin(); in != data.end(); ++in) {
    if (out->hkl != in->hkl || out->isign != in->isign) {
      out->value = sum_wI / sum_w;
      out->sigma = 1.0 / std::sqrt(sum_w);
      out->nobs = nobs;
      sum_wI = sum_w = 0.;
      nobs = 0;
      ++out;
      out->hkl = in->hkl;
      out->isign = in->isign;
    }
    double w = 1. / (in->sigma * in->sigma);
    sum_wI += w * in->value;
    sum_w += w;
    ++nobs;
  }
  out->value = sum_wI / sum_w;
  out->sigma = 1.0 / std::sqrt(sum_w);
  out->nobs = nobs;
  data.erase(++out, data.end());
}

void Intensities::switch_to_asu_indices(bool merged) {
  GroupOps gops = spacegroup->operations();
  ReciprocalAsu asu(spacegroup);
  for (Refl& refl : data) {
    if (asu.is_in(refl.hkl)) {
      if (!merged) {
        // isign is 0 for original hkl (e.g. from XDS file)
        if (refl.isign == 0)
          refl.isign = 1;  // since it's in asu - I+ or centric
        // when reading asu hkl from MTZ file - count centrics always as I+
        else if (refl.isign == -1 && gops.is_reflection_centric(refl.hkl))
          refl.isign = 1;
      }
      continue;
    }
    bool sign;
    std::tie(refl.hkl, sign) = asu.to_asu_sign(refl.hkl, gops);
    if (!merged)
      refl.isign = sign ? 1 : -1;
  }
}

void Intensities::read_unmerged_intensities_from_mtz(const Mtz& mtz) {
  if (mtz.batches.empty())
    fail("expected unmerged file");
  const Mtz::Column* isym_col = mtz.column_with_label("M/ISYM");
  if (!isym_col || isym_col->idx != 3)
    fail("unmerged file should have M/ISYM as 4th column");
  const Mtz::Column& col = mtz.get_column_with_label("I");
  size_t value_idx = col.idx;
  size_t sigma_idx = mtz.get_column_with_label("SIGI").idx;
  unit_cell = mtz.get_average_cell_from_batch_headers(unit_cell_rmsd);
  spacegroup = mtz.spacegroup;
  if (!spacegroup)
    fail("unknown space group");
  wavelength = mtz.dataset(col.dataset_id).wavelength;
  for (size_t i = 0; i < mtz.data.size(); i += mtz.columns.size()) {
    short isign = ((int)mtz.data[i + 3] % 2 == 0 ? -1 : 1);
    add_if_valid(*this, mtz.get_hkl(i), isign,
                 mtz.data[i + value_idx], mtz.data[i + sigma_idx]);
  }
  // Aimless >=0.7.6 (from 2021) has an option to output unmerged file
  // with original indices instead of reduced indices, with all ISYM = 1.
  switch_to_asu_indices();
  type = DataType::Unmerged;
}

void Intensities::read_mean_intensities_from_mtz(const Mtz& mtz) {
  if (!mtz.batches.empty())
    fail("expected merged file");
  const Mtz::Column* col = mtz.imean_column();
  if (!col)
    fail("Mean intensities (IMEAN, I, IOBS or I-obs) not found");
  size_t sigma_idx = mtz.get_column_with_label("SIG" + col->label).idx;
  copy_metadata(mtz, *this);
  wavelength = mtz.dataset(col->dataset_id).wavelength;
  read_data(*this, MtzDataProxy{mtz}, col->idx, sigma_idx);
  type = DataType::Mean;
}

void Intensities::read_anomalous_intensities_from_mtz(const Mtz& mtz, bool check_complete) {
  if (!mtz.batches.empty())
    fail("expected merged file");
  const Mtz::Column* colp = mtz.iplus_column();
  const Mtz::Column* colm = mtz.iminus_column();
  if (!colp || !colm)
    fail("anomalous intensities not found");
  size_t value_idx[2] = {colp->idx, colm->idx};
  size_t sigma_idx[2] = {mtz.get_column_with_label("SIG" + colp->label).idx,
                         mtz.get_column_with_label("SIG" + colm->label).idx};
  int mean_idx = -1;
  if (check_complete)
    if (const Mtz::Column* mean_col = mtz.imean_column())
      mean_idx = (int) mean_col->idx;
  copy_metadata(mtz, *this);
  wavelength = mtz.dataset(colp->dataset_id).wavelength;
  read_anomalous_data(*this, MtzDataProxy{mtz}, mean_idx, value_idx, sigma_idx);
  type = DataType::Anomalous;
}

static
void read_simple_intensities_from_mmcif(Intensities& intensities, const ReflnBlock& rb,
                                        const char* inten, const char* sigma) {
  size_t value_idx = rb.get_column_index(inten);
  size_t sigma_idx = rb.get_column_index(sigma);
  copy_metadata(rb, intensities);
  intensities.wavelength = rb.wavelength;
  read_data(intensities, ReflnDataProxy(rb), value_idx, sigma_idx);
}

void Intensities::read_unmerged_intensities_from_mmcif(const ReflnBlock& rb) {
  read_simple_intensities_from_mmcif(*this, rb, "intensity_net", "intensity_sigma");
  switch_to_asu_indices();
  type = DataType::Unmerged;
}

void Intensities::read_mean_intensities_from_mmcif(const ReflnBlock& rb) {
  read_simple_intensities_from_mmcif(*this, rb, "intensity_meas", "intensity_sigma");
  type = DataType::Mean;
}

void Intensities::read_anomalous_intensities_from_mmcif(const ReflnBlock& rb,
                                                        bool check_complete) {
  size_t value_idx[2] = {rb.get_column_index("pdbx_I_plus"),
                         rb.get_column_index("pdbx_I_minus")};
  size_t sigma_idx[2] = {rb.get_column_index("pdbx_I_plus_sigma"),
                         rb.get_column_index("pdbx_I_minus_sigma")};
  int mean_idx = -1;
  if (check_complete)
    mean_idx = rb.find_column_index("intensity_meas");
  copy_metadata(rb, *this);
  wavelength = rb.wavelength;
  read_anomalous_data(*this, ReflnDataProxy(rb), mean_idx, value_idx, sigma_idx);
  type = DataType::Anomalous;
}

void Intensities::read_f_squared_from_mmcif(const ReflnBlock& rb) {
  int value_idx = rb.find_column_index("F_meas");
  if (value_idx == -1)
    value_idx = rb.find_column_index("F_meas_au");
  if (value_idx == -1)
    fail("Column F_meas[_au] not found.");
  int sigma_idx = rb.find_column_index("F_meas_sigma");
  if (sigma_idx == -1)
    sigma_idx = rb.find_column_index("F_meas_sigma_au");
  if (sigma_idx == -1)
    fail("Column F_meas_sigma[_au] not found.");
  copy_metadata(rb, *this);
  wavelength = rb.wavelength;
  read_data(*this, ReflnDataProxy(rb), value_idx, sigma_idx);
  for (Refl& r : data) {
    r.value *= r.value;
    r.sigma *= 2 * r.value;
  }
  type = DataType::Mean;
}

void Intensities::read_unmerged_intensities_from_xds(const XdsAscii& xds) {
  unit_cell = xds.unit_cell;
  spacegroup = find_spacegroup_by_number(xds.spacegroup_number);
  wavelength = xds.wavelength;
  data.reserve(xds.data.size());
  for (const XdsAscii::Refl& in : xds.data)
    add_if_valid(*this, in.hkl, 0, in.iobs, in.sigma);
  switch_to_asu_indices();
  type = DataType::Unmerged;
}

std::string Intensities::take_staraniso_b_from_mtz(const Mtz& mtz) {
  return read_staraniso_b_from_mtz(mtz, staraniso_b.b);
}

bool Intensities::take_staraniso_b_from_mmcif(const cif::Block& block) {
  return read_staraniso_b_from_mmcif(block, staraniso_b.b);
}

Mtz Intensities::prepare_merged_mtz(bool with_nobs) {
  gemmi::Mtz mtz(/*with_base=*/true);
  mtz.spacegroup = spacegroup;
  mtz.set_cell_for_all(unit_cell);
  mtz.add_dataset("unknown").wavelength = wavelength;
  if (type == DataType::Mean) {
    mtz.add_column("IMEAN", 'J', -1, -1, false);
    mtz.add_column("SIGIMEAN", 'Q', -1, -1, false);
    if (with_nobs)
      mtz.add_column("NOBS", 'I', -1, -1, false);
  } else if (type == DataType::Anomalous) {
    mtz.add_column("I(+)", 'K', -1, -1, false);
    mtz.add_column("SIGI(+)", 'M', -1, -1, false);
    mtz.add_column("I(-)", 'K', -1, -1, false);
    mtz.add_column("SIGI(-)", 'M', -1, -1, false);
    if (with_nobs) {
      mtz.add_column("NOBS(+)", 'I', -1, -1, false);
      mtz.add_column("NOBS(-)", 'I', -1, -1, false);
    }
  } else {
    fail("prepare_merged_mtz(): data is not merged");
  }
  mtz.data.resize(data.size() * mtz.columns.size(), NAN);
  gemmi::Miller prev_hkl = data[0].hkl;
  mtz.set_hkl(0, prev_hkl);
  size_t offset = 0;
  for (const Intensities::Refl& refl : data) {
    if (refl.hkl != prev_hkl) {
      offset += mtz.columns.size();
      mtz.set_hkl(offset, refl.hkl);
      prev_hkl = refl.hkl;
    }
    size_t value_offset = offset + (refl.isign >= 0 ? 3 : 5);
    mtz.data[value_offset] = (float) refl.value;
    mtz.data[value_offset + 1] = (float) refl.sigma;
    if (with_nobs) {
      size_t nobs_offset = offset + 5;  // for "NOBS"
      if (type == DataType::Anomalous)
        nobs_offset += (refl.isign >= 0 ? 2 : 3);
      mtz.data[nobs_offset] = (float) refl.nobs;
    }
  }
  mtz.data.resize(offset + mtz.columns.size());
  mtz.nreflections = int(mtz.data.size() / mtz.columns.size());
  return mtz;
}

} // namespace gemmi

#if WITH_TEST
// Quick test for not exported function parse_voigt_notation().
// Compile with: c++ src/intensit.cpp -DWITH_TEST -Iinclude -lgemmi_cpp

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../third_party/doctest.h"
#include "gemmi/mtz2cif.hpp"  // for write_staraniso_b_in_mmcif
#include "gemmi/cif.hpp"  // for read_string

TEST_CASE("aniso_b_tensor_eigen") {
  std::string line = "(0.486, 17.6, 0.981, 3.004, -0.689, -1.99)";
  std::array<double,6> bval{0.486, 17.6, 0.981, 3.004, -0.689, -1.99};
  gemmi::SMat33<double> b1;
  bool ok = gemmi::parse_voigt_notation(line.c_str(), line.c_str() + line.size(), b1);
  CHECK(ok);
  CHECK_EQ(b1.elements_voigt(), bval);
  std::ostringstream os;
  os << "data_a\n";
  char buf[256];
  gemmi::write_staraniso_b_in_mmcif(b1, "xxxx", buf, os);
  gemmi::cif::Document doc = gemmi::cif::read_string(os.str());
  gemmi::SMat33<double> b2{0, 0, 0, 0, 0, 0};
  ok = gemmi::read_staraniso_b_from_mmcif(doc.blocks[0], b2);
  CHECK(ok);
  auto b2_elem = b2.elements_voigt();
  for (int i = 0; i < 6; ++i)
    CHECK(std::fabs(bval[i] - b2_elem[i]) < 1e-3);
}
#endif
