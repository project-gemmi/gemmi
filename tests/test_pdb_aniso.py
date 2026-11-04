import unittest

import gemmi


def structure_with_atom(aniso, pos):
    st = gemmi.Structure()
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")
    res = gemmi.Residue()
    res.name = "GLY"
    res.seqid = gemmi.SeqId("1")
    atom = gemmi.Atom()
    atom.name = "CA"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(*pos)
    atom.aniso = gemmi.SMat33f(*aniso)
    atom.occ = 1.0
    atom.b_iso = 20.0
    res.add_atom(atom)
    chain.add_residue(res)
    model.add_chain(chain)
    st.add_model(model)
    return st


def aniso_fields(st):
    line = next(line for line in st.make_pdb_string().splitlines()
                if line.startswith("ANISOU"))
    fields = line[28:70]
    return [fields[i:i + 7] for i in range(0, 42, 7)]


class TestPdbFixedColumns(unittest.TestCase):
    def test_aniso_negative_zero_stays_in_its_column(self):
        st = structure_with_atom((0.01, 0.01, 0.01, -0.00003, 0, 0),
                                 (1.0, 2.0, 3.0))
        fields = aniso_fields(st)
        self.assertEqual(fields[3], "      0")
        self.assertEqual(len("".join(fields)), 42)

    def test_aniso_overflow_does_not_shift_the_next_field(self):
        st = structure_with_atom((-1307, 0.02, 0.03, 0.001, 0.002, 0.004),
                                 (1.0, 2.0, 3.0))
        fields = aniso_fields(st)
        self.assertEqual(fields[0], "-999999")
        self.assertEqual(fields[1], "    200")

    def test_coordinate_exactly_at_the_rounding_boundary(self):
        st = structure_with_atom((0.01, 0.01, 0.01, 0.001, 0, 0),
                                 (1.0, 2.0, -0.0005))
        atom = next(line for line in st.make_pdb_string().splitlines()
                    if line.startswith("ATOM"))
        self.assertEqual(atom[46:54], "   0.000")

    def test_mtrix_negative_zero(self):
        st = gemmi.Structure()
        mat = gemmi.Mat33([[-1e-8, 0, 0], [0, 1, 0], [0, 0, 1]])
        st.ncs.append(gemmi.NcsOp(gemmi.Transform(mat, gemmi.Vec3(0, 0, -1e-7)),
                                  "1", True))
        kept = gemmi.Mat33([[-1e-6, 0, 0], [0, 1, 0], [0, 0, 1]])
        st.ncs.append(gemmi.NcsOp(gemmi.Transform(kept, gemmi.Vec3(0, 0, 0)),
                                  "2", True))
        lines = [line for line in st.make_pdb_string().splitlines()
                 if line.startswith("MTRIX")]
        self.assertIn("  0.000000", lines[0])
        self.assertNotIn("-0.000000", lines[0])
        self.assertNotIn("-0.00000", lines[0])
        self.assertIn("-0.000001", lines[3])


if __name__ == "__main__":
    unittest.main()
