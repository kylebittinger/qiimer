from unittest import TestCase, main
from relabel_otu_table import *

class AssignmentTest(TestCase):
    def test_parse_rdp(self):
        for rdp_str, exp_taxa in zip(rdp_output.split("\n"), rdp_assigned_taxa):
            self.assertEqual(Assignment.parse_rdp(rdp_str).taxa, exp_taxa)


rdp_output = '''\
Root;Bacteria;Firmicutes;"Clostridia";Clostridiales;Incertae Sedis XIII;Anaerovorax
Root;Bacteria;Bacteroidetes;Bacteroidetes;Bacteroidales;Bacteroidaceae;Bacteroides
Root;Bacteria;Proteobacteria
Root;Bacteria;Actinobacteria;Actinobacteria;Actinobacteridae;Actinomycetales;Actinomycineae;Actinomycetaceae;Actinomyces
'''

rdp_assigned_taxa = [
    {'class': 'Clostridia', 'phylum': 'Firmicutes',
     'genus': 'Anaerovorax', 'order': 'Incertae Sedis XIII'},
    {'family': 'Bacteroidaceae', 'class': 'Bacteroidetes',
     'phylum': 'Bacteroidetes', 'genus': 'Bacteroides',
     'order': 'Bacteroidales'},
    {'phylum': 'Proteobacteria'},
    {'family': 'Actinomycetaceae', 'subclass': 'Actinobacteridae',
     'order': 'Actinomycetales', 'phylum': 'Actinobacteria',
     'suborder': 'Actinomycineae', 'genus': 'Actinomyces',
     'class': 'Actinobacteria'},
    ]
for t in rdp_assigned_taxa:
    t['root'] = 'Root'
    t['domain'] = 'Bacteria'


if __name__ == '__main__':
    main() 
