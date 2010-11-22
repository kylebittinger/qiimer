import optparse
import re
import sys


class Assignment(object):
    ranks = [
        'root', 'domain', 'phylum', 'class', 'subclass', 'order', 'suborder',
        'family', 'subfamily', 'supergenus', 'genus',
        ]
    _suffixes = {
        'idae': 'subclass',
        'ales': 'order',
        'ineae': 'suborder',
        'aceae': 'family',
        'Chloroplast': 'order',
        'Family [XVI]+': 'family',
        'Caldilineacea': 'family',
        'Anaerolinaeles': 'order',
        'Anaerolinaeceea': 'family',
        }
    _next_taxon_guess = {
        'class': 'order',
        'subclass': 'order',
        'order': 'family',
        'suborder': 'family',
        'family': 'subfamily',
        'subfamily': 'supergenus',
        }

    def __init__(self, taxa=None):
        self.taxa = taxa or {}

    @property
    def lineage(self):
        for rank in self.ranks:
            yield (rank, self.taxa.get(rank))

    @classmethod
    def table_headers(cls):
        return ['taxon.%s' % s for r in cls.ranks]

    def table_vals(self):
        return [t for _, t in self.lineage]

    @classmethod
    def parse_rdp(cls, rdp_str):
        taxa = {}
        toks = [t.strip('"') for t in rdp_str.split(';')]

        try:
            taxa['root'] = toks.pop(0)
            taxa['domain'] = toks.pop(0)
            taxa['phylum'] = toks.pop(0)
            taxa['class'] = toks.pop(0)
            taxa['genus'] = toks.pop()
        except IndexError:
            return cls(taxa)

        last_taxon = 'class'
        supergenus_match = re.match("(.+) [a-z]$", taxa['genus'])

        mystery_toks = []
        for tok in toks:
            if 'family' in taxa:
                if re.search(taxa['family'] + " [0-9XVI]", tok):
                    taxa['subfamily'] = tok
                    last_taxon = 'subfamily'
                    continue
            if re.search('Incertae Sedis', tok, re.IGNORECASE):
                taxon = cls._next_taxon_guess[last_taxon]
                taxa[taxon] = tok
                last_taxon = taxon
                continue
            if supergenus_match:
                if tok == supergenus_match.group(1):
                    taxa['supergenus'] = tok
                    last_taxon = 'supergenus'
            for suffix, taxon in cls._suffixes.items():
                found_suffix = False
                if tok.endswith(suffix):
                    if taxon in taxa:
                        raise VauleError('Taxon %s assigned twice (%s, %s)' % (
                            taxon, taxa[taxon], tok))
                    taxa[taxon] = tok
                    found_suffix = True
                    break
            if not found_suffix:
                mystery_toks.append(tok)

        retval = cls(taxa)
        retval.mystery_toks = mystery_toks
        return retval


class OtuTableReformatter(object):
    def __init__(self, otu_table_file):
        self.file = otu_table_file

    @staticmethod
    def split_on_last_field(line):
        begin, _, end = line.rpartition('\t')
        return begin, end.strip()

    def reformat(self):
        comment_line = self.file.next()
        header_line = self.file.next()
        _, last_header = self.split_on_last_field(header_line)
        if last_header != 'Consensus Lineage':
            raise ValueError(
                'This OTU table does not contain a lineage field. '
                'Last header: %s' % last_header
                )
        yield header_line
        for line in self.file:
            begin, rdp_str = self.split_on_last_field(line)
            assignment = Assignment.parse_rdp(rdp_str)
            label = self.make_label(assignment)
            yield begin + '\t' + label + '\n'

    def make_label(self, assignment):
        taxa = assignment.taxa
        if 'phylum' in taxa:
            return '%s %s' % (taxa['phylum'], taxa.get('family', 'Unassigned'))
        elif 'domain' in taxa:
            return '%s Unknown' % taxa['domain']
        else:
            return 'Unassigned'


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-i', '--input_otu_table_fp', help='Input OTU table')
    parser.add_option('-o', '--output_otu_table_fp', help='Output OTU table')
    opts, args = parser.parse_args()

    input_file = open(opts.input_otu_table_fp)
    formatter = OtuTableReformatter(input_file)

    if opts.output_otu_table_fp:
        output_file = open(opts.output_otu_table_fp, 'w')
    else:
        output_file = sys.stdout
    for line in formatter.reformat():
        output_file.write(line)

