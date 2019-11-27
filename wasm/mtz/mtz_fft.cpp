// Copyright 2019 Global Phasing Ltd.

#define POCKETFFT_CACHE_SIZE 0
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#include <cstdlib>            // for free
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/fourier.hpp>  // for update_cif_block
#include <gemmi/math.hpp>     // for Variance
#include <emscripten/bind.h>

using gemmi::Mtz;

class MtzFft {
public:
  // char* or void* cannot be used
  // https://github.com/emscripten-core/emscripten/issues/9448
  MtzFft(int32_t data, size_t size)
    : data_((char*) data) {
    try {
      mtz_.read_stream(gemmi::MemoryStream(data_, data_ + size), false);
    } catch (std::runtime_error& e) {
      (void) e;
    }
  }

  ~MtzFft() {
    std::free(data_);
  }

  int32_t calculate_map_from_columns(const Mtz::Column* f_col,
                                     const Mtz::Column* phi_col) {
    try {
      const float* raw_data = (const float*)(data_ + 80);
      gemmi::MtzExternalDataProxy proxy(mtz_, raw_data);
      auto size = gemmi::get_size_for_hkl(proxy, {{0, 0, 0}}, 3.);
      gemmi::Grid<std::complex<float>> coefs
        = gemmi::get_f_phi_on_grid<float>(proxy, f_col->idx, phi_col->idx,
                                          size, /*half_l=*/true,
                                          gemmi::HklOrient::LKH);
      gemmi::transform_f_phi_grid_to_map_(std::move(coefs), grid_);
    } catch (std::runtime_error& e) {
      std::puts(e.what());
      return 0;
    }
    gemmi::Variance grid_variance(grid_.data.begin(), grid_.data.end());
    rmsd_ = std::sqrt(grid_variance.for_population());
    return (int32_t) grid_.data.data();
  }

  int32_t calculate_map_from_labels(const std::string& f_label,
                                    const std::string& phi_label) {
    if (const Mtz::Column* f_col = mtz_.column_with_label(f_label))
      if (const Mtz::Column* phi_col = mtz_.column_with_label(phi_label))
        return calculate_map_from_columns(f_col, phi_col);
    std::puts("Requested labels not found in the MTZ file.");
    return 0;
  }

  int32_t calculate_map(bool diff_map) {
    static const char* normal_labels[] = {
      "FWT", "PHWT",
      "2FOFCWT", "PH2FOFCWT",
    };
    static const char* diff_labels[] = {
      "DELFWT", "PHDELWT",
      "FOFCWT", "PHFOFCWT",
    };
    const char** labels = diff_map ? diff_labels : normal_labels;
    for (int i = 0; i < 4; i += 2)
      if (const Mtz::Column* f_col = mtz_.column_with_label(labels[i]))
        if (const Mtz::Column* phi_col = mtz_.column_with_label(labels[i+1]))
          return calculate_map_from_columns(f_col, phi_col);
    std::puts("Default map coefficient labels not found");
    return 0;
  }

  int get_nx() const { return grid_.nu; }
  int get_ny() const { return grid_.nv; }
  int get_nz() const { return grid_.nw; }

  double get_rmsd() const { return rmsd_; }

  double get_cell_param(int n) const {
    switch (n) {
      case 0: return mtz_.cell.a;
      case 1: return mtz_.cell.b;
      case 2: return mtz_.cell.c;
      case 3: return mtz_.cell.alpha;
      case 4: return mtz_.cell.beta;
      case 5: return mtz_.cell.gamma;
    }
    return 0.;
  }

private:
  gemmi::Mtz mtz_;
  char* data_;
  gemmi::Grid<float> grid_;
  double rmsd_;
};

using emscripten::class_;
using emscripten::allow_raw_pointers;

EMSCRIPTEN_BINDINGS(GemmiMtz) {
  class_<MtzFft>("Mtz")
    .constructor<int32_t, size_t>(allow_raw_pointers())
    .function("calculate_map", &MtzFft::calculate_map, allow_raw_pointers())
    .property("nx", &MtzFft::get_nx)
    .property("ny", &MtzFft::get_ny)
    .property("nz", &MtzFft::get_nz)
    .property("rmsd", &MtzFft::get_rmsd)
    .function("cell_param", &MtzFft::get_cell_param)
    ;
}
