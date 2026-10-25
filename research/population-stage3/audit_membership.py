"""Compare proposed tutorial identities with LVK Table 1's BBH subset.

This checks identifiers, not PE sample identity or full selection equivalence.
"""
from html.parser import HTMLParser
from pathlib import Path
import hashlib
import json
import re

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
ALIASES = {'GW190412': 'GW190412_053044', 'GW190521': 'GW190521_030229',
           'GW190425': 'GW190425_081805', 'GW190814': 'GW190814_211039'}
# Sections III, V (Table 3) and VI of arXiv:2111.03634v3.
EXCLUDED = {'GW170817', 'GW190425_081805', 'GW190426_152155',
            'GW190814_211039', 'GW190917_114630', 'GW200105_162426',
            'GW200115_042309'}

class TableNames(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_table_figure = False
        self.in_cell = False
        self.cells = []
        self.names = []

    def handle_starttag(self, tag, attrs):
        if dict(attrs).get('id') == 'S2.T1':
            self.in_table_figure = True
        if self.in_table_figure and tag == 'td':
            self.in_cell = True
            self.cells = []

    def handle_data(self, data):
        if self.in_cell:
            self.cells.append(data)

    def handle_endtag(self, tag):
        if tag == 'td' and self.in_cell:
            name = ''.join(self.cells).strip()
            if re.fullmatch(r'GW\d{6}(?:_\d{6})?', name):
                self.names.append(name)
            self.in_cell = False
        if tag == 'figure':
            self.in_table_figure = False

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def main():
    source = BASE / 'sources/population-2111.03634v3.html'
    parser = TableNames()
    parser.feed(source.read_text(encoding='utf-8'))
    table = set(parser.names)
    assert len(parser.names) == len(table) == 76
    assert EXCLUDED <= table
    assert 'GW190531' not in table  # Caption only; no PE, not a table row.
    official = table - EXCLUDED
    audit = json.loads((ROOT / 'data/sample_audit.json').read_text())
    selected = [x for x in audit if x['included']]
    local = {ALIASES.get(x['name'], x['name']) for x in selected}
    assert len(local) == len(selected) == 69
    assert len(official) == 69
    result = dict(schema='blackhole-population-event-membership/1',
        source_url='https://arxiv.org/html/2111.03634v3#S2.T1',
        source_sha256=digest(source), local_audit_sha256=digest(ROOT / 'data/sample_audit.json'),
        aliases=ALIASES, official_table_count=len(table),
        official_non_bbh_exclusions=sorted(EXCLUDED), official_bbh_names=sorted(official),
        local_selected_count=len(local), missing_from_local=sorted(official-local),
        additional_in_local=sorted(local-official), exact_identifier_set_match=(local == official),
        local_identity_basis='summary matching; sample-level evidence is reported separately',
        scope='Identifier membership only; does not establish matching waveforms, priors or detection efficiency.',
        local_rows=[dict(index=x['index'], proposed_name=x['name'],
                        canonical_name=ALIASES.get(x['name'], x['name'])) for x in selected])
    (BASE / 'event-membership.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print(json.dumps({k: v for k,v in result.items() if k not in ['local_rows','official_bbh_names']}, indent=2))

if __name__ == '__main__':
    main()
