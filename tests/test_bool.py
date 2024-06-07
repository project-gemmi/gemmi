import gemmi 
import unittest

class TestBool(unittest.TestCase):

    def test_all(self):
        s = gemmi.read_structure("5wkd.pdb")

        for m in s:
            self.assertTrue(m, "Model is empty") 
            for c in m: 
                self.assertTrue(c, "Chain is empty") 
                for r in c:
                    self.assertTrue(r, "Residue is empty") 


    def test_structure(self):
        s = gemmi.read_structure("5wkd.pdb")
        self.assertTrue(s, "Structure is empty") 
        self.assertTrue(s[0], "Model is empty") 
        self.assertTrue(s[0][0], "Chain is empty") 

        s.add_model(gemmi.Model("A"))
        self.assertFalse(s[-1], "Model is not empty")
    
    def test_model(self):
        m = gemmi.Model("A")
        self.assertFalse(m, "Model is not empty")
        m.add_chain(gemmi.Chain("A"))
        self.assertTrue(m, "Model is empty")
        del m[0]
        self.assertFalse(m, "Model is not empty")

    def test_chain(self):
        c = gemmi.Chain("A")
        self.assertFalse(c, "Chain is not empty")
        c.add_residue(gemmi.Residue())
        self.assertTrue(c, "Chain is empty")
        del c[0]
        self.assertFalse(c, "Chain is not empty")

    def test_residue(self):
        r = gemmi.Residue()
        self.assertFalse(r, "Residue is not empty")
        r.add_atom(gemmi.Atom())
        self.assertTrue(r, "Residue is empty")
        del r[0]
        self.assertFalse(r, "Residue is not empty")


