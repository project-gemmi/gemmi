#!/usr/bin/env python

import gc
import os
import unittest
from gemmi import cif

class TestDoc(unittest.TestCase):
    def test_slice(self):
        doc = cif.read_string("""
            data_a
            _one 1 _two 2 _three 3
            data_b
            _four 4
            data_c
            _two 2 _four 4 _six 6
        """)
        self.assertEqual([b.name for b in doc[:1]], ['a'])
        self.assertEqual([b.name for b in doc[1:]], ['b', 'c'])
        self.assertEqual([b.name for b in doc[:]], ['a', 'b', 'c'])
        self.assertEqual([b.name for b in doc[1:-1]], ['b'])
        self.assertEqual([b.name for b in doc[1:1]], [])

    def test_contains(self):
        doc = cif.read_string("""
            data_a
            _one 1 _two 2 _three 3
            data_b
            _four 4
            data_c
            _two 2 _four 4 _six 6
        """)
        self.assertEqual('a' in doc, True)
        self.assertEqual('d' in doc, False)

class TestBlock(unittest.TestCase):
    def test_find(self):
        block = cif.read_string("""
            data_test
            _one 1 _two 2 _three 3
            _nonloop_a alpha
            _nonloop_b beta
            loop_ _la _lb _ln A B 1 C D 2
        """).sole_block()
        gc.collect()
        values = block.find_values('_la')
        self.assertEqual(list(values), ['A', 'C'])

        rows = list(block.find(['_lb', '_la']))  # changed order
        gc.collect()
        self.assertEqual(list(rows[0]), ['B', 'A'])
        self.assertEqual(rows[1][1], 'C')

        self.assertEqual(block.get_index('_nonloop_b'), 4)
        self.assertIsNone(block[4].loop)
        self.assertEqual(block[4].pair[0], '_nonloop_b')
        self.assertEqual(block.get_index('_lb'), 5)
        self.assertEqual(block[5].loop.tags[1], '_lb')
        self.assertIsNone(block[5].pair)

        rows = list(block.find('_nonloop_', ['a', 'b']))
        gc.collect()
        self.assertEqual([list(r) for r in rows], [['alpha', 'beta']])

        tab = block.find('_l', ['a', 'n'])
        self.assertEqual(tab.find_row('A')[1], '1')
        self.assertRaises(RuntimeError, tab.find_row, 'B')
        self.assertEqual(tab.find_row('C')[1], '2')
        self.assertEqual(tab.column(0)[1], 'C')
        self.assertEqual(tab.column(1)[0], '1')
        self.assertEqual(tab.find_column('a')[1], 'C')
        self.assertEqual(tab.find_column('_ln')[0], '1')

        tab = block.find(['_lb', '_ln'])
        self.assertRaises(RuntimeError, tab.find_row, 'A')
        self.assertEqual(tab.find_row('B')[1], '1')
        self.assertEqual(tab.find_row('D')[1], '2')

        tab = block.find(['_la', '?_nop', '?_ln'])
        self.assertEqual(len(tab), 2)
        self.assertEqual(tab.width(), 3)
        row = tab[0]
        self.assertIsNotNone(row.get(0))
        self.assertEqual(row[0], 'A')
        self.assertIsNone(row.get(1))
        self.assertEqual(row[2], '1')
        self.assertTrue("None" in repr(row))

    def test_find_values(self):
        v1 = cif.read_string('data_a _v1 one')[0].find_values('_v1')
        gc.collect()
        self.assertListEqual(list(v1), ['one'])
        v2 = cif.read_string('data_a loop_ _v2 a b')[0].find_values('_v2')
        gc.collect()
        self.assertListEqual(list(v2), ['a', 'b'])

    def test_set_pair(self):
        doc = cif.read_string('data_a _a 1 _b 2 _c 3')
        block = doc[0]
        block.set_pair('_d', '9')
        block.set_pair('_b', '8')
        self.assertEqual(block.find_pair('_a'), ('_a', '1'))
        self.assertEqual(block.find_value('_b'), '8')
        self.assertEqual(block.find_value('_c'), '3')
        self.assertEqual(block.find_value('_d'), '9')

    def test_set_pairs(self):
        doc = cif.read_string('data_a _zza 0 _z_a 1 _Z_b 2 _z_C 3 _zzb 4')
        block = doc[0]
        expected = [('_zza', '0'), ('_z_a', '1'), ('_Z_b', '2'),
                    ('_z_C', '3'), ('_zzb', '4')]
        self.assertEqual([v.pair for v in block], expected)
        block.set_pairs('_z_', {'1': 5})
        expected.insert(4, ('_z_1', '5'))
        self.assertEqual([v.pair for v in block], expected)
        block.set_pairs('_z_', {'B': 10, 'c' : 11, 'a': 9})
        expected[1:4] = [('_z_a', '9'), ('_z_B', '10'), ('_z_c', '11')]
        self.assertEqual([v.pair for v in block], expected)

    def test_set_loop(self):
        block = cif.read_string('data_a _c.a 1 _c.b 2 _c.c 3 loop_ _cx.b 3')[0]
        block.init_loop('_cx.', ['b']).add_row(['x'])
        self.assertEqual(block.find_value('_c.a'), '1')
        block.init_mmcif_loop('_c.', ['c']).add_row(['y'])
        self.assertEqual(block.find_value('_c.a'), None)
        loop = block.init_loop('_c.', ['c', 'd'])
        loop.set_all_values([list('one'), list('two')])
        self.assertEqual(loop.width(), 2)
        self.assertEqual(loop.length(), 3)
        self.assertEqual(list(block.find(['_c.c', '_c.d'])[1]), ['n', 'w'])
        tab = block.find_or_add('_c.', ['d'])
        self.assertEqual((tab.width(), len(tab)), (1, 3))
        tab = block.find_or_add('_c.', ['d', 'c'])
        self.assertEqual((tab.width(), len(tab)), (2, 3))
        tab = block.find_or_add('_u.', ['d', 'c'])
        self.assertEqual((tab.width(), len(tab)), (2, 0))

    def test_erase(self):
        block = cif.read_string('data_a _a_a 1 _a_b 2 _a_c 3 loop_ _x 4')[0]
        block.find('_a_', ['a', '?x', '?c', '?bb']).erase()
        self.assertEqual(block.as_string().split(),
                         ['data_a', '_a_b', '2', 'loop_', '_x', '4'])

    def test_setitem(self):
        block = cif.read_string('data_a _a 1 _b 2 _c 3')[0]
        self.assertEqual(block.find_value('_b'), '2')
        col_b = block.find_values('_b')
        col_b[0] = '20'
        self.assertEqual(block.find_value('_b'), '20')
        bc = block.find(['_b', '_a'])
        self.assertEqual(bc[0]['_a'], '1')
        bc[0][0] = '30'
        self.assertEqual(block.find_value('_b'), '30')
        bc[0][1] = '40'
        self.assertEqual(block.find_value('_a'), '40')
        self.assertEqual(bc[0]['_a'], '40')
        bc[0]['_a'] = '44'
        self.assertEqual(block.find_value('_a'), '44')
        self.assertEqual(block.find_value('_a'), '44')
        self.assertEqual(block.find_value('_b'), '30')

    def test_add_row(self):
        block = cif.read_string('data_a loop_ _x _y 1 2 3 4')[0]
        loop = block.find_loop('_x').get_loop()
        loop.add_row(['5', '6'])
        self.assertEqual(loop.values, ['1', '2', '3', '4', '5', '6'])
        loop.add_row(['?', '0'], 0)
        self.assertEqual(list(block.find_values('_y')), '0 2 4 6'.split())
        self.assertEqual(loop.length(), 4)
        block.find(['_x']).append_row(['xa'])
        block.find(['_y']).append_row(['ya'])
        block.find(['_y', '_x']).append_row(['yb', 'xb'])
        loop.add_columns(column_names=['_a', '_z'], value='A')
        loop.remove_column('_a')
        self.assertEqual(loop.tags, ['_x', '_y', '_z'])
        block.find(['_x', '_y']).append_row(['xc', 'yc'])
        self.assertEqual(loop.length(), 8)
        self.assertEqual(list(block.find_values('_x')),
                         '? 1 3 5 xa . xb xc'.split())
        self.assertEqual(list(block.find_values('_y')),
                         '0 2 4 6 . ya yb yc'.split())
        self.assertEqual(list(block.find_values('_z')),
                         'A A A A A A A .'.split())

    def test_set_mmcif_category(self):
        doc = cif.Document()
        block = doc.add_new_block('b')
        block.set_mmcif_category('_c', {
            'one': ('?', 'ab', ';text field\n;'),
            'two': [-1, 4./3, '"double quoted"']}, raw=True)
        self.assertEqual(block.find_values('_c.one')[0], '?')
        self.assertEqual(block.find_values('_c.one').str(1), 'ab')
        self.assertEqual(block.find_values('_c.one').str(2), 'text field')
        self.assertEqual(block.find_values('_c.two').str(0), '-1')
        self.assertEqual(block.find_values('_c.two').str(2), 'double quoted')
        block.set_mmcif_category('_d', {
            'one': (None, 'a b', 'text\nfield'),
            'two': [-1, '?', False]})
        def check_d():
            self.assertEqual(block.find_values('_d.one')[0], '?')
            self.assertEqual(block.find_values('_d.one').str(1), 'a b')
            self.assertEqual(block.find_values('_d.one').str(2), 'text\nfield')
            self.assertEqual(block.find_values('_d.one')[2], ';text\nfield\n;')
            self.assertEqual(block.find_values('_d.two').str(0), '-1')
            self.assertEqual(block.find_values('_d.two').str(1), '?')
            self.assertEqual(block.find_values('_d.two')[2], '.')
        check_d()
        block.set_mmcif_category('_D', {
            'one': ('?', "'a b'", ';text\nfield\n;'),
            'two': ['-1', "'?'", '.']},
            raw=True)
        self.assertEqual(set(block.find_mmcif_category('_d').tags),
                         {'_D.one', '_D.two'})
        check_d()
        block.set_mmcif_category('_d', {
            'one': (None, "'a b'", ';text\nfield\n;'),
            'two': [-1, "'?'", False]},
            raw=True)
        check_d()
        block.set_mmcif_category('_d', block.get_mmcif_category('_D'))
        check_d()
        block.set_mmcif_category('_d', block.get_mmcif_category('_d', raw=True),
                                 raw=True)
        check_d()
        block.set_mmcif_category('_d', {})
        self.assertEqual(block.find_mmcif_category('_d').width(), 0)

    def test_mmcif_file(self):
        path = os.path.join(os.path.dirname(__file__), '5i55.cif')
        block = cif.read(path).sole_block()
        self.assertEqual(len(block.get_mmcif_category_names()), 54)
        entry_cat = block.get_mmcif_category('_entry')
        self.assertEqual(entry_cat, {'id': ['5I55']})
        drw_cat = block.get_mmcif_category('_diffrn_radiation_wavelength.')
        self.assertEqual(drw_cat, {'id': ['1', '2', '3'],
                                   'wavelength': ['0.9792', '0.9794', '0.9796'],
                                   'wt': ['1.0']*3})
        cc_cat = block.get_mmcif_category('_chem_comp.')
        self.assertEqual(cc_cat['mon_nstd_flag'][:2], [False, 'y'])
        self.assertEqual(cc_cat['pdbx_synonyms'][:2], [None, None])

        # test __delitem__
        del block.find(['_entry.id'])[0]
        entry_cat = block.get_mmcif_category('_entry')
        self.assertEqual(entry_cat, {'id': []})
        def nums():
            return list(block.find_values('_entity_poly_seq.num'))
        tab = block.find(['_entity_poly_seq.mon_id'])
        self.assertEqual(nums(), [str(i) for i in range(1, 23)])
        del tab[1::2]
        self.assertEqual(nums(), [str(i) for i in range(1, 23, 2)])
        del tab[3:]
        self.assertEqual(nums(), ['1', '3', '5'])
        del tab[:-1]
        self.assertEqual(nums(), ['5'])

    def test_reading_gzipped_file(self):
        path = os.path.join(os.path.dirname(__file__), '1pfe.cif.gz')
        cif_doc = cif.read(path)
        block = cif_doc.sole_block()
        categories = block.get_mmcif_category_names()
        self.assertEqual(categories[0], '_entry.')
        self.assertEqual(len(categories), 72)
        exptl = block.find_mmcif_category('_exptl')
        self.assertEqual(exptl.get_prefix(), '_exptl.')
        self.assertEqual(list(exptl.tags), ['_exptl.entry_id', '_exptl.method',
                                            '_exptl.crystals_number'])
        self.assertEqual(len(exptl), 1)
        self.assertEqual(exptl.width(), 3)
        exptl = block.find_mmcif_category('_exptl')
        self.assertEqual(len(exptl), 1)
        self.assertEqual(exptl.width(), 3)
        self.assertEqual(exptl[0].str(1), 'X-RAY DIFFRACTION')
        struct_asym = block.find_mmcif_category('_struct_asym')
        self.assertEqual(len(struct_asym), 7)
        self.assertEqual(struct_asym.width(), 5)
        self.assertListEqual(list(struct_asym[3]), ['D', 'N', 'N', '4', '?'])
        nonexistent = block.find_mmcif_category('_nonexistent')
        self.assertRaises(RuntimeError, nonexistent.get_prefix)
        self.assertEqual(len(nonexistent), 0)
        self.assertEqual(nonexistent.width(), 0)

    def test_file_not_found(self):
        with self.assertRaises(IOError):
            cif.read('file-that-does-not-exist.cif')

    def test_syntax_error(self):
        with self.assertRaises(ValueError):
            cif.read_string('data_a boom')

    def test_line_endings(self):
        lines = ['data_a', '_a_field', ';', 'text line', '2nd line', ';']
        for eol in ('\r\n', '\n'):
            input_str = eol.join(lines) + eol
            doc = cif.read_string(input_str)
            output_str = doc.as_string()
            self.assertEqual(input_str.replace(eol, '\n'), output_str)

    def test_text_field_eol(self):
        # use fragment of ma-bak-cepc-0017.cif
        path = os.path.join(os.path.dirname(__file__), 'eol-test.cif')
        with open(path) as f:
            cif_string = f.read()
        doc = cif.read_string(cif_string)
        output = doc.as_string()
        self.assertEqual([s.rstrip() for s in cif_string.splitlines()],
                         output.splitlines())

    def test_write_style(self):
        doc = cif.read_string('data_one _x y')
        self.assertEqual(doc.as_string(), doc[0].as_string())
        self.assertEqual(doc.as_string().splitlines(), ['data_one', '_x y'])
        options = cif.WriteOptions()
        options.align_pairs = 33
        self.assertEqual(doc.as_string(options).splitlines(),
                         ['data_one', '_x                                y'])

    def test_relion_syntax_exception(self):
        block = cif.read_string("""\
            data_
            loop_
            _rlnImageName
            _rlnMicrographName
            _rlnDefocusU
            _rlnDefocusV
            mic1/img000001.spi mic1 10000 10500
            mic1/img000002.spi mic1 10000 10501
            """).sole_block()
        self.assertEqual(block.name, ' ')
        item = block[0]
        self.assertEqual(item.line_number, 2)
        table = block.item_as_table(item)
        self.assertEqual(table[1][3], '10501')

    def test_case_sensitivity(self):
        block = cif.read_string("""
            daTA_test
            _One 1 _two 2 _thrEE 3
            _NonLoop_a alpha
            _nonloop_B beta
            loop_ _laAa _lbBb _ln A B 1 C D 2
        """).sole_block()
        values = block.find_values('_laaA')
        self.assertEqual(list(values), ['A', 'C'])

        rows = list(block.find(['_lBbb', '_lAaa']))  # changed order
        self.assertEqual(list(rows[0]), ['B', 'A'])

        self.assertEqual(block.get_index('_nonlOOp_b'), 4)
        self.assertEqual(block.get_index('_lBBB'), 5)
        self.assertEqual(block.find_pair('_Three'), ('_thrEE', '3'))

        values.erase()
        block.find_values('_two').erase()
        block.find_values('_nonloop_b').erase()
        expected = """\
            data_test
            _One 1 _thrEE 3
            _NonLoop_a alpha
            loop_ _lbBb _ln  B 1  D 2"""
        self.assertEqual(block.as_string().split(), expected.split())

class TestQuote(unittest.TestCase):
    def test_quote(self):
        self.assertEqual(cif.quote('a.b-c'), 'a.b-c')
        self.assertEqual(cif.quote('loop_'), "'loop_'")
        self.assertEqual(cif.quote('a"\nb'), ';a"\nb\n;')
        self.assertEqual(cif.quote(u'\u0394'), u"'\u0394'")

def full_path(filename):
    return os.path.join(os.path.dirname(__file__), filename)

class TestDictinary(unittest.TestCase):
    def test_frame_reading(self):
        block = cif.read(full_path('mmcif_pdbx_v50_frag.dic')).sole_block()
        self.assertIsNone(block.find_frame('heyho'))
        frame = block.find_frame('_atom_site.auth_atom_id')
        code = frame.find_value('_item_type.code')
        self.assertEqual(code, 'atcode')

if __name__ == '__main__':
    unittest.main()
