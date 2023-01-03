#include <gemmi/mtz.hpp>

int main() {
  gemmi::Mtz mtz(/*with_base=*/true);
  mtz.title = "My new MTZ file";
  mtz.history.push_back("Created from scratch");
  mtz.spacegroup = gemmi::find_spacegroup_by_name("P 21 21 21");
  mtz.set_cell_for_all(gemmi::UnitCell(77.7, 149.5, 62.4, 90, 90, 90));
  mtz.add_dataset("synthetic");
  mtz.datasets.back().wavelength = 0.8;
  mtz.add_column("F", 'F', -1, -1, false);
  mtz.add_column("SIGF", 'Q', -1, -1, false);
  const float data[] = { 2, 3, 4, 200.4f, 10.5f,
                         2, 3, 5, 596.1f, 7.35f };
  mtz.set_data(data, 2*5);
  mtz.write_to_file("my_new.mtz");
}
