#!/usr/bin/env python3

# Check for new versions of third-party libraries

from __future__ import print_function
import datetime
import json
from urllib.request import urlopen

TAGGED_REPOS = {
    'wjakob/nanobind': 'v2.5.0',
    'scikit-build/scikit-build-core': 'v0.11.1',
    'taocpp/PEGTL': '2.8.3',
    'cxong/tinydir': '1.2.6',
    'fastfloat/fast_float': 'v6.1.6',
    'madler/zlib': 'v1.3.1',
    'LLNL/shroud': 'v0.13.0',
}

NOT_TAGGED_REPOS = {
    'chadaustin/sajson': 'include/sajson.h',
    'sgorsten/linalg': 'linalg.h',
    'nothings/stb': 'stb_sprintf.h',
}

def load_json(url):
    with urlopen(url) as response:
        return json.load(response)

def check_tags():
    for repo, version in TAGGED_REPOS.items():
        try:
            url = 'https://api.github.com/repos/%s/releases/latest' % repo
            data = load_json(url)
            latest_tag = data['tag_name']
        except IOError:
            data = load_json('https://api.github.com/repos/%s/tags' % repo)
            latest_tag = data[0]['name']
        mark = ('   !!!' if version != latest_tag else '')
        print('%-18s %10s %10s%s' % (repo[-18:], version, latest_tag, mark))

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

try:
    check_tags()
    check_recent_commits()
except IOError as e:
    print(e)
