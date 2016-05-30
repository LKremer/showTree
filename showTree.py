#!/usr/bin/python
'''
Visualization of gene trees produced by Carsten Kemenas geneSearch and
calcTree tools. Next to the gene trees, a MSA including protein domains
is displayed.
Works with python2.7 and python3.3+ :)
by Lukas PM Kremer, 2016
'''


from __future__ import print_function
import argparse
from ete3 import Tree, TreeStyle, SeqMotifFace, NodeStyle
from Bio import SeqIO


COLORS = [
    'Khaki', 'DarkSeaGreen', 'LightSteelBlue', 'orange', 'yellow',
    'cyan', 'pink', 'magenta', 'green', 'blue', 'Moccasin', 'brown',
    'gray', 'black', 'white',
] * 3


class DomTreeMaker(object):
    def __init__(self, msa_path, config_path, tree_path, root=None,
                 out=False, scale_factor=1.0, highlight=None,
                 no_alignment=False, hide_nodes=False):
        self.msa_path = msa_path
        self.config_path = config_path
        self.tree_path = tree_path
        self.out = out
        self.scale_factor = scale_factor
        self.termnodes_to_highlight = self.parse_highlight_option(highlight)
        self.outgroup_node_for_rooting = root
        self.no_alignment = no_alignment
        self.hide_nodes = hide_nodes

    def parse_highlight_option(self, hl_list):
        termnodes_to_highlight = []
        if hl_list:
            for hl_str in hl_list:

                comma_sep = ',' in hl_str
                plus_sep = '+' in hl_str
                if plus_sep:
                    sep = '+'
                elif comma_sep:
                    sep = ','
                elif comma_sep and plus_sep:
                    raise Exception('Cannot use "+" and "," in '
                        'the same highlight option.')
                else:
                    sep = 'doesnt_matter'

                colon_split = hl_str.split(':')
                bg_color = colon_split[0]
                prot_ids = set(colon_split[1].split(sep))
                assert len(colon_split) == 2 and '' not in prot_ids, (
                    '\n  Invalid argument format for parameter '
                    '-hl/--highlight\n  '
                    'The correct format is for instance:\n    '
                    '--highlight yellow:prot01,prot02,prot03 red:prot04\nor\n'
                    '--highlight yellow:prot01+prot02+prot03 red:prot04')
                termnodes_to_highlight.append(
                    (bg_color, prot_ids, plus_sep)
                )
        return termnodes_to_highlight

    def parse_config(self):
        '''
        get the paths to the domain annotations
        from the geneSearch config file.
        '''
        self.dom_files = set()
        with open(self.config_path, 'r') as cfg:
            for line in cfg:
                if line.startswith('#') or \
                        line.startswith('ortho') or \
                        line.startswith('augustus') or \
                        line.startswith('!targets')or \
                        not line.strip():
                    continue
                dom_file_path = line.strip().split()[2]
                self.dom_files.add(dom_file_path)
        print('Found these domain annotations:')
        print('\t' + '\n\t'.join(self.dom_files) + '\n')

    def get_gene_list(self):
        t = Tree(self.tree_path)
        gene_set = set()
        for leaf in t:
            gene_set.add(leaf.name)
        return gene_set

    def parse_msa(self):
        '''
        store the MSA multi-fasta as a dict with header as keys and
        sequence as values.
        '''
        self.msa_fasta_dict = {}
        with open(self.msa_path, 'rU') as fa:
            for record in SeqIO.parse(fa, 'fasta'):
                header = str(record.id)
                seq = str(record.seq)
                seq_id = header.split()[0]
                self.msa_fasta_dict[seq_id] = seq

    def collect_domain_info(self):
        '''
        parse all domain files that were gathered from the config file
        and collect relevant domain info.
        '''
        self.domain_dict = {}
        self.domains = set()
        for dom_annot_path in self.dom_files:
            self.parse_domain_annotation(dom_annot_path)
        self.domains = list(self.domains)

    def parse_domain_annotation(self, dom_annot_path):
        '''
        Parses a pfam_scan domain annotation file.
        Collects domain names and start/end borders for all genes that
        are contained in the MSA file.
        '''
        this_annot = set()
        gene_set = self.get_gene_list()
        with open(dom_annot_path, 'r') as df:
            for line in df:
                if line.startswith('#') or not line.strip():
                    continue
                fields = line.split()
                gene_id = fields[0].replace(':', '_')
                if gene_id in self.domain_dict and gene_id not in this_annot:
                    print('\033[31m Warning! The ID {} is listed in '
                          'multiple domain annotations!\033[0m'.format(gene_id))
                this_annot.add(gene_id)

                if gene_id in gene_set:
                    ungapped_start, ungapped_stop = \
                        int(fields[1]), int(fields[2])
                    ungapped_start -= 1  # adjusting borders to Pythons 0-indexing
                    ungapped_stop -= 1   # (Pfam domains are 1-indexed)
                    domname = str(fields[6])
                    self.domains.add(domname)
                    if self.msa_fasta_dict:
                        gapped_seq = self.msa_fasta_dict[gene_id]


                        try:
                            start, stop = self.correct_borders_for_gaps(
                                gapped_seq, ungapped_start, ungapped_stop
                            )
                        except KeyError:
                            raise DomainIndexError(
                                dom=domname, annot=dom_annot_path,
                                plen=len(gapped_seq.replace('-', '')),
                                b1=ungapped_start, b2=ungapped_stop,
                                prot=gene_id
                            )
                        if gene_id not in self.domain_dict:
                            self.domain_dict[gene_id] = []
                        self.domain_dict[gene_id].append(
                            (domname, start, stop))
                    else:
                        if gene_id not in self.domain_dict:
                            self.domain_dict[gene_id] = []
                        self.domain_dict[gene_id].append(
                            (domname, 1, 1))

    def add_background_color_to_nodes(self, t, whole_clade=True):
        for bg_color, prot_ids, whole_clade in self.termnodes_to_highlight:
            for prot_id in prot_ids:
                assert prot_id in t, ('The node "{}" that you want to '
                    'highlight was not found in your tree!'.format(prot_id))
            nstyle = NodeStyle()
            nstyle['bgcolor'] = bg_color
            if not whole_clade or len(prot_ids) == 1:
                # don't colour the whole clade containing the proteins,
                # only color the proteins individually:
                for prot_id in prot_ids:
                    (t & prot_id).set_style(nstyle)
            else:
                # colour the whole clade that contains the proteins:
                node = t.get_common_ancestor(*prot_ids)
                node.set_style(nstyle)

    def root_tree(self, t):
        if '+' in self.outgroup_node_for_rooting:
            node_ids = self.outgroup_node_for_rooting.split('+')
            ancestor = t.get_common_ancestor(*node_ids)
            t.set_outgroup(ancestor)
            print('Rooted tree at the clade that contains {}' \
                '.\n'.format(' and '.join(node_ids)))
        else:
            t.set_outgroup(t&self.outgroup_node_for_rooting)
            print('Rooted tree at node "{}"' \
                '.\n'.format(self.outgroup_node_for_rooting))

    def show_dom_MSA_tree(self):
        '''
        Plot the gene tree with a MSA+domains aligned next to it,
        or alternatively dump the tree to a file (if self.out was defined)
        '''
        t = Tree(self.tree_path)
        ts = TreeStyle()
        ts.show_leaf_name = True
        ts.show_branch_length = True
        ts.show_branch_support = True
        ts.scale = 120  # 120 pixels per branch length unit

        if self.outgroup_node_for_rooting:
            self.root_tree(t)

        if not self.no_alignment:
            if self.msa_path:
                self.add_background_color_to_nodes(t)
                self.add_msa_with_domains_to_tree(t)
            else:
                self.add_background_color_to_nodes(t)
                self.add_domains_to_tree(t)

        if self.hide_nodes:
            for word2hide in self.hide_nodes:
                nodes_hidden = []
                print('\nRemoving these terminal nodes because they contain the '\
                    'word "{}":'.format(word2hide))
                for leaf in t.traverse():
                    if word2hide in leaf.name:
                        nodes_hidden.append(leaf.name)
                        leaf.delete()
                print(', '.join(nodes_hidden))

        if self.out:
            if not self.out.upper().endswith('.PDF'):
                self.out += '.pdf'
            t.render(self.out, tree_style=ts)
            print('Wrote output image to', self.out)
        else:
            t.show(tree_style=ts)

    def add_domains_to_tree(self, t):
        for leaf in t:
            gene_id = leaf.name
            domains = self.domain_dict.get(gene_id, [])
            # if no domains are annotated, 'domains' is an empty list and
            # no motifs are added to the sequence (for loop won't iterate)
            motifs = []
            dom_n = 0
            for domain in domains:
                domname, start, stop = domain
                start = dom_n * 55
                stop = dom_n * 55 + 50
                dom_n += 1
                color_n = self.domains.index(domname)
                try:
                    dom_color = COLORS[color_n]
                except IndexError:
                    dom_color = 'gray'
                motifs.append([
                    start, stop, '()', None, 10, None,
                    'rgradient:{}'.format(dom_color),
                    'arial|6|black|{}'.format(domname)
                ])
            if len(motifs) > 0:
                domface = SeqMotifFace(None, motifs=motifs, gap_format='line')
                (t & gene_id).add_face(domface, column=0, position='aligned')


    def add_msa_with_domains_to_tree(self, t):
        for leaf in t:
            gene_id = leaf.name
            gapped_seq = self.msa_fasta_dict[gene_id]
            domains = self.domain_dict.get(gene_id, [])
            # if no domains are annotated, 'domains' is an empty list and
            # no motifs are added to the sequence (for loop won't iterate)
            motifs = []
            for domain in domains:
                domname, start, stop = domain
                color_n = self.domains.index(domname)
                try:
                    dom_color = COLORS[color_n]
                except IndexError:
                    dom_color = 'gray'
                motifs.append([
                    start, stop, '()', None, 10, None,
                    'rgradient:{}'.format(dom_color),
                    'arial|6|black|{}'.format(domname)
                ])
            seqface = SeqMotifFace(gapped_seq, gapcolor='gray',
                                    seq_format='compactseq',
                                    scale_factor=self.scale_factor)
            (t & gene_id).add_face(seqface, column=0, position='aligned')
            domface = SeqMotifFace(gapped_seq, motifs=motifs,
                                    seq_format='blank', gap_format='blank',
                                    scale_factor=self.scale_factor)
            (t & gene_id).add_face(domface, column=0, position='aligned')
        

    def map_gapped_to_ungapped_position(self, seq):
        position_map = {}
        corrected_pos = 0
        for original_pos, char in enumerate(seq):
            position_map[corrected_pos] = original_pos
            if char == '-':
                continue
            corrected_pos += 1
        position_map[corrected_pos] = original_pos
        assert corrected_pos == len(seq.replace('-', ''))
        return position_map

    def correct_borders_for_gaps(self, seq, start, stop):
        '''
        seq: sequence with gaps, e.g. A----VTIAGISVML---HMVSPMSQ---
        start/stop: start/stop position of a domain, but according to
            an ungapped sequence e.g. AVTIAGISVMLHMVSPMSQ
        returns: corrected start/stop (according to the gapped sequence)
        '''
        gapped_to_ungapped_pos = self.map_gapped_to_ungapped_position(seq)
        corrected_start = gapped_to_ungapped_pos[start]
        corrected_stop = gapped_to_ungapped_pos[stop]
        dom_seq = seq.replace('-', '')[start:stop+1]
        dom_seq_from_new_borders = \
            seq[corrected_start:corrected_stop+1]
        if dom_seq_from_new_borders.endswith('-'):
            corrected_stop = gapped_to_ungapped_pos[stop]
            dom_seq_from_new_borders = \
                seq[corrected_start:corrected_stop+1]
        assert dom_seq == dom_seq_from_new_borders.replace('-', ''), 'Could '\
            'not determine the domain borders for the gapped MSA'
        return corrected_start, corrected_stop


class DomainIndexError(Exception):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        e_msg = '''\n
  The domain annotation
  {annot}
  lists domain borders that are out of range
  of the annotated protein {prot} :

  protein length: {plen}
  domain:         {dom}
  domain borders: {b1} - {b2}
        '''.format(**kwargs)
        Exception.__init__(self, e_msg)


def main(msa_path, config_path, tree_path, out=False,
         scale_factor=1.0, highlight=None, root=None,
         no_alignment=False, hide_nodes=False):
    d = DomTreeMaker(
         msa_path=msa_path,
         config_path=config_path,
         tree_path=tree_path,
         out=out,
         scale_factor=scale_factor,
         highlight=highlight,
         root=root,
         no_alignment=no_alignment,
         hide_nodes=hide_nodes
    )
    d.parse_config()
    if msa_path != None:
        d.parse_msa()
    else:
        d.msa_fasta_dict = None
    d.collect_domain_info()
    d.show_dom_MSA_tree()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required arguments')
    required.add_argument('-c', '--config', required=True,
                          help='Path to the geneSearch/calcTree '
                          'configuration file')
    required.add_argument('-m', '--msa', required=False,
                          help='Path to the untrimmed T-Coffee multiple '
                          'sequence alignment produced by calcTree, '
                          'usually ends with "_aln.fa"', default=None)
    required.add_argument('-t', '--tree', required=True,
                          help='Path to the best-scoring RAxML tree with '
                          'support values (not as branch labels) produced '
                          'by calcTree, usually named '
                          '"RAxML_bipartitions.(...)_aln.tree"')
    parser.add_argument('-o', '--output_path', default=False,
                        help='Path to the output image file (.PDF). '
                        'Tree will be shown in a window if omitted')
    parser.add_argument('-s', '--scale_factor', default=1.0, type=float,
                        help='Horizontal scaling factor of the MSA '
                        '(default: 1.0). '
                        'Decrease this value if your image is too wide!')
    parser.add_argument('-hl', '--highlight', nargs='*',
                        help='Highlight specific clades or terminal nodes '
                        'with a background color. '
                        'Comma-separated IDs will colour each terminal '
                        'node separately. '
                        'Plus-separated IDs will colour the whole clade '
                        'that contains those IDs. '
                        'Example: --highlight red:prot01,prot02,prot03 '
                        'green:prot04 lightyellow:prot05+prot06')
    parser.add_argument('-r', '--root', help='Set an outgroup node at which '
                        'the tree will be rooted. If multiple IDs are '
                        'specified (plus-separated), the whole clade that '
                        'contains these IDs will be used for rooting')
    parser.add_argument('-n', '--no_alignment', help='Do not draw the '
                        'multiple sequence alignment or domains',
                        action='store_true')
    parser.add_argument('-hn', '--hide_nodes', help='Hide terminal nodes '
                        'that contain the specified text.', nargs='+')
    
    args = parser.parse_args()

    main(
         msa_path=args.msa,
         config_path=args.config,
         tree_path=args.tree,
         out=args.output_path,
         scale_factor=args.scale_factor,
         highlight=args.highlight,
         root=args.root,
         no_alignment=args.no_alignment,
         hide_nodes=args.hide_nodes,
    )
