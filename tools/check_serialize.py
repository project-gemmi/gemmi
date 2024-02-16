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
    template_params = []
    fields = []
    for child in node.get_children():
        if child.kind == cindex.CursorKind.TEMPLATE_TYPE_PARAMETER:
            template_params.append('typename')
        elif child.kind == cindex.CursorKind.TEMPLATE_NON_TYPE_PARAMETER:
            first_token = next(child.get_tokens())
            template_params.append(first_token.spelling)
        elif child.kind == cindex.CursorKind.CXX_BASE_SPECIFIER:
            (base_ref,) = tuple(child.get_children())
            assert base_ref.kind == cindex.CursorKind.TYPE_REF
            struct_base = base_ref.referenced.spelling
        elif child.kind == cindex.CursorKind.FIELD_DECL:
            # public = (child.access_specifier == cindex.AccessSpecifier.PUBLIC)
            fields.append(child.spelling)
    if node.kind == cindex.CursorKind.CLASS_TEMPLATE:
        assert len(template_params) == 1
        beginning = f'SERIALIZE_T1({name}, {template_params[0]}'
    else:
        if struct_base:
            beginning = f'SERIALIZE_P({name}, {struct_base}'
        else:
            beginning = f'SERIALIZE({name}'
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
    checked_kinds = [
        cindex.CursorKind.STRUCT_DECL,
        cindex.CursorKind.CLASS_DECL,
        cindex.CursorKind.CLASS_TEMPLATE,
    ]
    for g1 in tu.cursor.get_children():
        if g1.kind == cindex.CursorKind.NAMESPACE and g1.spelling == 'gemmi':
            for g2 in g1.get_children():
                name = g2.spelling
                if g2.kind in checked_kinds and name in macros:
                    compare_struct(name, g2, macros[name])
                    del macros[name]
                    for g3 in g2.get_children():
                        nested = f'{name}::{g3.spelling}'
                        if g3.kind in checked_kinds and nested in macros:
                            compare_struct(nested, g3, macros[nested])
                            del macros[nested]
    if macros:
        print('\nUNCHECKED:', ' '.join(macros))

if __name__ == '__main__':
    main()
