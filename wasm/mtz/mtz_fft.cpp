// Copyright 2019 Global Phasing Ltd.

#define POCKETFFT_CACHE_SIZE 0
#define POCKETFFT_NO_MULTITHREADING
#define POCKETFFT_NO_VECTORS
#include <cstdlib>            // for free
#include <gemmi/mtz.hpp>      // for Mtz
#include <gemmi/fourier.hpp>  // for update_cif_block
#include <gemmi/math.hpp>     // for Variance
#include <emscripten/bind.h>
#include <emscripten/val.h>

using gemmi::Mtz;
using emscripten::val;

class MtzFft {
public:
  MtzFft() {}
  // char* or void* cannot be used
  // https://github.com/emscripten-core/emscripten/issues/9448
  bool read(int32_t data, size_t size) {
    data_ = (char*) data;
    try {
      mtz_.read_stream(gemmi::MemoryStream(data_, size), false);
    } catch (std::runtime_error& e) {
      last_error_ = "Failed to read MTZ file: ";
      last_error_ += e.what();
      return false;
    }
    return true;
  }

  ~MtzFft() {
    std::free(data_);
  }

  val calculate_map_from_columns(const Mtz::Column* f_col,
                                 const Mtz::Column* phi_col) {
    try {
      const float* raw_data = (const float*)(data_ + 80);
      gemmi::MtzExternalDataProxy proxy(mtz_, raw_data);
      auto size = gemmi::get_size_for_hkl(proxy, {{0, 0, 0}}, 3.);
      using namespace gemmi;
      FPhiProxy<MtzExternalDataProxy> fphi(proxy, f_col->idx, phi_col->idx);
      FPhiGrid<float> coefs = get_f_phi_on_grid<float>(fphi, size, /*half_l=*/true,
                                                       gemmi::AxisOrder::ZYX);
      gemmi::transform_f_phi_grid_to_map_(std::move(coefs), grid_);
    } catch (std::runtime_error& e) {
      last_error_ = e.what();
      return val::null();
    }
    const std::vector<float>& data = grid_.data;
    gemmi::Variance grid_variance(data.begin(), data.end());
    rmsd_ = std::sqrt(grid_variance.for_population());
    return val(emscripten::typed_memory_view(data.size(), data.data()));
    //return (int32_t) data.data();
  }

  val calculate_map_from_labels(const std::string& f_label,
                                const std::string& phi_label) {
    if (const Mtz::Column* f_col = mtz_.column_with_label(f_label))
      if (const Mtz::Column* phi_col = mtz_.column_with_label(phi_label))
        return calculate_map_from_columns(f_col, phi_col);
    last_error_ = "Requested labels not found in the MTZ file.";
    return val::null();
  }

  val calculate_map(bool diff_map) {
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
    last_error_ = "Default map coefficient labels not found";
    return val::null();
  }

  // we used AxisOrder::ZYX
  int get_nx() const { return grid_.nw; }
  int get_ny() const { return grid_.nv; }
  int get_nz() const { return grid_.nu; }

  double get_rmsd() const { return rmsd_; }

  std::string get_last_error() const { return last_error_; }

  gemmi::UnitCell get_cell() const { return mtz_.cell; }

private:
  gemmi::Mtz mtz_;
  char* data_ = nullptr;
  gemmi::Grid<float> grid_;
  double rmsd_ = 0.;
  std::string last_error_;
};

using emscripten::class_;
using emscripten::allow_raw_pointers;

EMSCRIPTEN_BINDINGS(GemmiMtz) {
  class_<gemmi::UnitCell>("UnitCell")
    .property("a", &gemmi::UnitCell::a)
    .property("b", &gemmi::UnitCell::b)
    .property("c", &gemmi::UnitCell::c)
    .property("alpha", &gemmi::UnitCell::alpha)
    .property("beta", &gemmi::UnitCell::beta)
    .property("gamma", &gemmi::UnitCell::gamma)
    ;

  class_<MtzFft>("Mtz")
    .constructor<>()
    .function("read", &MtzFft::read, allow_raw_pointers())
    .function("calculate_map", &MtzFft::calculate_map)
    .function("calculate_map_from_labels", &MtzFft::calculate_map_from_labels)
    .property("nx", &MtzFft::get_nx)
    .property("ny", &MtzFft::get_ny)
    .property("nz", &MtzFft::get_nz)
    .property("rmsd", &MtzFft::get_rmsd)
    .property("last_error", &MtzFft::get_last_error)
    .property("cell", &MtzFft::get_cell)
    ;
}
