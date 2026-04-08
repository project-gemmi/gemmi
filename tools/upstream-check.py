#!/usr/bin/env python3

# Check for new versions of third-party libraries

from __future__ import print_function
import datetime
import json
from pathlib import Path
import re
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
TAGS_ONLY_REPOS = {
    'emscripten-core/emscripten',
    'mrdoob/three.js',
}

TAGGED_REPOS = {
    'wjakob/nanobind': 'v2.6.1',
    'scikit-build/scikit-build-core': 'v0.11.5',
    'taocpp/PEGTL': '2.8.3',
    'cxong/tinydir': '1.2.6',
    'fastfloat/fast_float': 'v8.2.4',
    'madler/zlib': 'v1.3.2',
    'LLNL/shroud': 'v0.14.0',
}

NOT_TAGGED_REPOS = {
    'chadaustin/sajson': 'include/sajson.h',
    'sgorsten/linalg': 'linalg.h',
    'nothings/stb': 'stb_sprintf.h',
}

def load_json(url):
    with urlopen(url) as response:
        return json.load(response)

def load_latest_tag(repo):
    if repo in TAGS_ONLY_REPOS:
        data = load_json('https://api.github.com/repos/%s/tags' % repo)
        return data[0]['name']
    try:
        url = 'https://api.github.com/repos/%s/releases/latest' % repo
        data = load_json(url)
        return data['tag_name']
    except IOError:
        data = load_json('https://api.github.com/repos/%s/tags' % repo)
        return data[0]['name']

def extract_first_group(path, pattern):
    if not path.is_file():
        return '-'
    match = re.search(pattern, path.read_text(encoding='utf-8'), re.MULTILINE)
    return match.group(1) if match else '?'

def get_emscripten_version():
    return extract_first_group(ROOT / 'wasm' / 'Makefile',
                               r'^# currently using Emscripten (\S+)$')

def get_threejs_version():
    revision = extract_first_group(Path.home() / 'fresh' / 'three.js' /
                                   'package.json',
                                   r'"version"\s*:\s*"0\.(\d+)\.[^"]+"')
    if revision not in ('-', '?'):
        return 'r' + revision
    src_dir = ROOT.parent / 'gemmimol' / 'src'
    if not src_dir.is_dir():
        return '-'
    best = None
    for path in src_dir.iterdir():
        if not path.is_dir():
            continue
        match = re.match(r'^three-(r\d+)$', path.name)
        if not match:
            continue
        revision = int(match.group(1)[1:])
        if best is None or revision > best[0]:
            best = (revision, match.group(1))
    return best[1] if best is not None else '-'

def local_tagged_versions():
    for repo, version in TAGGED_REPOS.items():
        yield repo[-18:], repo, version
    yield 'Emscripten', 'emscripten-core/emscripten', get_emscripten_version()
    yield 'three.js', 'mrdoob/three.js', get_threejs_version()

def check_tags():
    for label, repo, version in local_tagged_versions():
        latest_tag = load_latest_tag(repo)
        mark = ('   !!!' if version not in ('-', '?') and version != latest_tag
                else '')
        print('%-18s %10s %10s%s' % (label, version, latest_tag, mark))

def check_recent_commits():
    since = datetime.datetime.now() - datetime.timedelta(days=90)
    for repo, filename in NOT_TAGGED_REPOS.items():
        url = 'https://api.github.com/repos/%s/commits?path=%s&since=%s' % (
              repo, filename, since.isoformat())
        data = load_json(url)
        info = '-'
        if data:
            info = data[0]['commit']['committer']['date'][:10]
        print('%-18s  modified: %s' % (filename.split('/')[-1], info))

def main():
    try:
        check_tags()
        check_recent_commits()
    except IOError as e:
        print(e)

if __name__ == '__main__':
    main()
