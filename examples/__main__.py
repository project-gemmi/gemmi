from __future__ import print_function
import os

print('Run directly one of the example scripts:\n')
example_dir = os.path.dirname(__file__)
for d in sorted(os.listdir(example_dir)):
    if d[0] != '_' and d.endswith('.py'):
        print(os.path.join(example_dir, d))
print('\nCheck Gemmi documentation for more details,')
print('or try running examples with the option --help (-h for short).')
