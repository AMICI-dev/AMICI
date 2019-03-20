#!/usr/bin/env python3
"""Create AMICI publication list by year"""

import biblib.bib
import biblib.messages
import biblib.algo
import sys
import subprocess


def get_keys_by_year(bibfile):
    """Get bibtex entry keys as dict by year"""

    with open(bibfile, 'r') as f:
        db = biblib.bib.Parser().parse(f, log_fp=sys.stderr).get_entries()
        recoverer = biblib.messages.InputErrorRecoverer()
        by_year = {}
        for ent in db.values():
            with recoverer:
                if 'year' in ent:
                    try:
                        by_year[ent['year']].append(ent.key)
                    except KeyError:
                        by_year[ent['year']] = [ent.key]
                else:
                    print("Missing year for entry", ent.key)
        recoverer.reraise()
    return by_year


def get_sub_bibliography(year, by_year):
    """Get HTML bibliography for the given year"""

    entries = ','.join(['@' + x for x in by_year[year]])
    input = '---\n' \
            'bibliography: amici_refs.bib\n' \
        f'nocite: "{entries}"\n...\n' \
        f'# {year}'

    out = subprocess.run(['pandoc', '--filter=pandoc-citeproc',
                          '-f', 'markdown'],
                         input=input, capture_output=True,
                         encoding='utf-8')
    return out.stdout


def main():
    bibfile = 'amici_refs.bib'
    outfile = 'references.md'

    by_year = get_keys_by_year(bibfile)

    with open(outfile, 'w') as f:
        for year in reversed(sorted(by_year.keys())):
            cur_bib = get_sub_bibliography(year, by_year)
            f.write(cur_bib)


if __name__ == '__main__':
    main()
