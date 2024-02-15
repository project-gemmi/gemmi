# Compares serialization macros from the SERIALIZE_PATH file
# with libclang's AST of the TU_PATH translation unit.

from clang import cindex

SERIALIZE_PATH = 'include/gemmi/serialize.hpp'
TU_PATH = 'include/gemmi/model.hpp'

def read_macros_from_file():
    macros = {}
    with open(SERIALIZE_PATH) as f:
        name = None
        for line in f:
            line = line.strip()
            if line.startswith('SERIALIZE'):
                assert name is None
                name = line[line.index('(')+1 : line.index(',')]
                macros[name] = line
            elif name:
                macros[name] += ' ' + line
            if name and ')' in line:
                name = None
    return macros

def compare_struct(name, node, line_from_file):
    print(name, end=' ')
    struct_base = None
    fields = []
    for child in node.get_children():
        if child.kind == cindex.CursorKind.CXX_BASE_SPECIFIER:
            (base_ref,) = tuple(child.get_children())
            assert base_ref.kind == cindex.CursorKind.TYPE_REF
            struct_base = base_ref.referenced.spelling
        elif child.kind == cindex.CursorKind.FIELD_DECL:
            # public = (child.access_specifier == cindex.AccessSpecifier.PUBLIC)
            fields.append(child.spelling)
    if struct_base:
        beginning = f'SERIALIZE_P({node.spelling}, {struct_base}'
    else:
        beginning = f'SERIALIZE({node.spelling}'
    expected = beginning + ''.join(', o.' + field for field in fields) + ')'
    if expected == line_from_file:
        print('OK')
    else:
        print('differs:')
        print(' ', expected)
        print(' ', line_from_file)

def main():
    macros = read_macros_from_file()
    print(f'{SERIALIZE_PATH} has {len(macros)} SERIALIZE* macros')

    index = cindex.Index.create()
    print(f'Parsing {TU_PATH} ...')
    tu = index.parse(TU_PATH)
    for g1 in tu.cursor.get_children():
        if g1.kind == cindex.CursorKind.NAMESPACE and g1.spelling == 'gemmi':
            for g2 in g1.get_children():
                name = g2.spelling
                if g2.kind == cindex.CursorKind.STRUCT_DECL and name in macros:
                    compare_struct(name, g2, macros[name])
                    del macros[name]
    print('Unchecked:', ' '.join(macros))

if __name__ == '__main__':
    main()
