#!/usr/bin/python
'''
Visualization of gene trees produced by Carsten Kemenas geneSearch and
calcTree tools. Next to the gene trees, a MSA including protein domains
is displayed.
Works with python2.7 and python3.3+ :)
by Lukas PM Kremer, 2016/17
'''


from __future__ import print_function
import argparse
import os
from ete3 import Tree, TreeStyle, SeqMotifFace, NodeStyle
from Bio import SeqIO


COLORS = [
    'Khaki', 'DarkSeaGreen', 'LightSteelBlue', 'orange', 'yellow',
    'cyan', 'pink', 'magenta', 'green', 'blue', 'Moccasin', 'brown',
    'gray', 'black', 'white',
] * 3


class DomTreeMaker(object):
    def __init__(self, msa_path, tree_path=None, config_path=None, root=None,
                 out=False, scale_factor=1.0, highlight=None, hide_nodes=False,
                 domain_annotation=None, gff_path=None, lineseq=False):
        self.msa_path = msa_path
        self.gff_path = gff_path
        self.domain_annotation = domain_annotation
        self.config_path = config_path
        self.tree_path = tree_path
        self.out = out
        self.scale_factor = scale_factor
        self.termnodes_to_highlight = self.parse_highlight_option(highlight)
        self.outgroup_node_for_rooting = root
        self.hide_nodes = hide_nodes
        if lineseq:
            print("line")
            self.seq_style = 'line'
            self.gap_style = 'blank'
        else:
            print("compactseq")
            self.seq_style = 'compactseq'
            self.gap_style = 'line'
        return

    def parse_highlight_option(self, hl_list):
        '''
        parse the --highlight option to determine
        which nodes / clades should be marked in
        which colour.
        '''
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
                    raise Exception(
                        'Cannot use "+" and "," in the same highlight option.'
                    )
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

    def get_domain_annot_paths(self):
        '''
        get the paths to the domain annotation(s) from the geneSearch
        config file, and/or directly from the --domain_annotation param.
        '''
        self.dom_files = set()
        if self.config_path:
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
        if self.domain_annotation:
            print('You manually specified the domain annotation(s):')
            for dom_annot in self.domain_annotation:
                print('\t{}'.format(dom_annot))
                self.dom_files.add(dom_annot)
            print()
        if not self.domain_annotation and not self.config_path:
            print('\nYou did not specify a geneSearch config or '
                  'domain annotation file, so domains will not be '
                  'displayed.\n')
        return

    def get_gene_list(self):
        '''
        extracts the protein IDs from the fasta alignment file
        or from the tree (whatever is available). Otherwise,
        the domain annotation itself will be used to get the IDs.
        '''
        if not hasattr(self, '_gene_list'):
            self._gene_list = set()
            if self.msa_path:
                with open(self.msa_path, 'r') as msa_file:
                    for line in msa_file:
                        if line.startswith('>'):
                            ID = line.split()[0][1:]
                            self._gene_list.add(ID)
            elif self.tree_path:
                t = Tree(self.tree_path)
                self._gene_list = set()
                for leaf in t:
                    self._gene_list.add(leaf.name)
            elif self.domain_annotation:
                if self.domain_dict:
                    self._gene_list = self.domain_dict.keys()
                else:
                    self._gene_list = None
            else:
                raise Exception('This should never happen')
        return self._gene_list

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
        return

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
        return

    def parse_domain_annotation(self, dom_annot_path):
        '''
        Parses a pfam_scan domain annotation file.
        Collects domain names and start/end borders for all genes that
        are contained in the MSA file.
        '''
        this_annot = set()
        gene_set = self.get_gene_list()  # will be None if neither MSA nor tree were specified
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

                if not gene_set or gene_id in gene_set:  # if no MSA or tree were specified, all pfam_scan entries will be plotted
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
        return

    def parse_gff(self, gff_path):
        ''' parse a GFF file to determine intron positions '''
        if not hasattr(self, 'cds_length_dict'):
            self.cds_length_dict = {}
        gene_list = self.get_gene_list()
        found = False
        n_weird_lines = 0  # how many lines could not be parsed/are not in GFF format
        with open(gff_path, 'r') as gff_file:
            for line in gff_file:
                if line.startswith('#') or not line.strip():
                    continue
                values = line.strip().split('\t')
                if not is_valid_gff_line(values):
                    n_weird_lines += 1
                    continue
                feature = values[2]
                if feature.upper() != 'CDS':
                    continue
                gff_comment = values[8]
                for query_gene in gene_list:
                    if query_gene in gff_comment:
                        if query_gene not in self.cds_length_dict:
                            self.cds_length_dict[query_gene] = []
                        cds_len = (int(values[4]) - int(values[3])) / 3.0
                        if cds_len > 3:  # some GFFs have CDSs of length 0, why?!
                            self.cds_length_dict[query_gene].append(cds_len)
                        found = True
        if n_weird_lines:
            print('\nWarning! The GFF {} contained {} lines that could not be '
                  'interpreted! Check the file.\n'.format(gff_path, n_weird_lines))
        assert found, '\n\nDidn\'t find any relevant proteins in {}.\n\n'\
            'Looked for these:\n{}'.format(gff_path, '\t'.join(gene_list))
        return

    def add_background_color_to_nodes(self, t, whole_clade=True):
        for bg_color, prot_ids, whole_clade in self.termnodes_to_highlight:
            for prot_id in prot_ids:
                assert prot_id in t, (
                    'The node "{}" that you want to '
                    'highlight was not found in your tree!'.format(prot_id)
                )
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
        return

    def root_tree(self, t):
        if '+' in self.outgroup_node_for_rooting:
            node_ids = self.outgroup_node_for_rooting.split('+')
            ancestor = t.get_common_ancestor(*node_ids)
            t.set_outgroup(ancestor)
            print('Rooted tree at the clade that contains {}'
                  '.\n'.format(' and '.join(node_ids)))
        else:
            t.set_outgroup(t & self.outgroup_node_for_rooting)
            print('Rooted tree at node "{}"'
                  '.\n'.format(self.outgroup_node_for_rooting))
        return

    def show_dom_MSA_tree(self):
        '''
        Plot the gene tree with a MSA+domains aligned next to it,
        or alternatively dump the tree to a file (if self.out was defined)
        '''
        if not self.tree_path:
            # if no tree was specified, use a "star-shaped" newick tree
            # as a substitute (all species have the same relatedness)
            newick_star = '({});'.format(','.join(sorted(self.get_gene_list())))
            t = Tree(newick_star)
            ts = TreeStyle()
            ts.show_branch_length = False
            ts.show_branch_support = False
            ts.branch_vertical_margin = 10
            ts.scale = 0  # 0 pixels per branch length unit
            # make the tree disappear:
            nstyle = NodeStyle()
            nstyle['size'] = 0
            nstyle['hz_line_color'] = 'white'
            nstyle['vt_line_color'] = 'white'
            for n in t.traverse():
                n.set_style(nstyle)
        else:
            t = Tree(self.tree_path)
            ts = TreeStyle()
            ts.show_branch_length = True
            ts.show_branch_support = True
            ts.scale = 120  # 120 pixels per branch length unit
        ts.show_leaf_name = True
        ts.show_scale = False

        if self.outgroup_node_for_rooting:
            assert self.tree_path, ('Rooting doesn\'t make sense if '
                                    'you don\'t specify a tree file!')
            self.root_tree(t)

        if self.msa_path:
            self.add_background_color_to_nodes(t)
            self.add_msa_with_domains_to_tree(t)
        else:
            self.add_background_color_to_nodes(t)
            self.add_domains_to_tree(t)

        if self.hide_nodes:
            for word2hide in self.hide_nodes:
                nodes_hidden = []
                print('\nRemoving these terminal nodes because they contain the '
                      'word "{}":'.format(word2hide))
                for leaf in t.traverse():
                    if word2hide in leaf.name:
                        nodes_hidden.append(leaf.name)
                        leaf.delete()
                print(', '.join(nodes_hidden))

        if self.out:
            file_ext = os.path.splitext(self.out.upper())[1]
            if file_ext not in ('.SVG', '.PDF', '.PNG'):
                raise Exception('Output image file name must end with '
                    '".SVG", ".PDF" or ".PNG".')
            t.render(self.out, tree_style=ts)
            print('Wrote output image to', self.out)
        else:
            t.show(tree_style=ts)
        return

    def add_domains_to_tree(self, t):
        '''
        displaying domains without a sequence / MSA
        '''
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
        return

    def add_msa_with_domains_to_tree(self, t):
        '''
        iterating over all sequences in the tree and adding
        the sequence, domain and intron position visualizations.
        '''
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
                ])  # domain markers

            if hasattr(self, 'cds_length_dict'):
                if gene_id not in self.cds_length_dict:
                    print('Warning: No GFF entry found for {}'.format(gene_id))
                else:
                    cds_lengths = self.cds_length_dict[gene_id]
                    current_pos = 1
                    for cds_len in cds_lengths[:-1]:
                        current_pos += cds_len

                        gapped_seq = self.msa_fasta_dict[gene_id]
                        try:
                            gapped_pos, __ = self.correct_borders_for_gaps(
                                gapped_seq, int(round(current_pos)), 0
                            )
                        except KeyError:
                            raise Exception('The protein sequence of {} is shorter '
                                'in the MSA than in the GFF!'.format(gene_id))

                        motifs.append([
                            (gapped_pos - 1), (gapped_pos + 1),
                            '[]', None, 10, None, 'black', None
                        ])  # black line that marks the intron positions

            seqface = SeqMotifFace(gapped_seq, gapcolor='gray', motifs=motifs,
                                   seq_format=self.seq_style, gap_format=self.gap_style,
                                   scale_factor=self.scale_factor)
            (t & gene_id).add_face(seqface, column=0, position='aligned')
        return

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

def is_valid_gff_line(splitline):
    ''' check if a file line (split) is in valid GFF format '''
    if not len(splitline) == 9:
        return False
    try:
        feat_len = int(splitline[4]) - int(splitline[3])
        assert feat_len >= 0
    except (ValueError, AssertionError):
        return False
    return True


def main(msa_path, config_path, tree_path, gff_paths=None,
         out=False, scale_factor=1.0, highlight=None,
         root=None, hide_nodes=False, domain_annotation=None,
         lineseq=False):
    d = DomTreeMaker(
         msa_path=msa_path,
         config_path=config_path,
         tree_path=tree_path,
         gff_path=gff_paths,
         out=out,
         scale_factor=scale_factor,
         highlight=highlight,
         root=root,
         hide_nodes=hide_nodes,
         domain_annotation=domain_annotation,
         lineseq=lineseq,
    )
    d.get_domain_annot_paths()
    if msa_path is not None:
        d.parse_msa()
    else:  # why was this required again? ...
        d.msa_fasta_dict = None
    d.collect_domain_info()
    if gff_paths is not None:
        for gff_path in gff_paths:
            d.parse_gff(gff_path)
    d.show_dom_MSA_tree()
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main_args = parser.add_argument_group('main arguments')
    add_args = parser.add_argument_group('additional optional arguments')
    main_args.add_argument('-m', '--msa', default=None,
                           help='Path to an untrimmed multiple sequence '
                           'alignment in FASTA format (the one produced '
                           'by calcTree usually ends with "_aln.fa")')
    main_args.add_argument('-t', '--tree', default=None,
                           help='Path to a gene tree, i.e. the best-'
                           'scoring gene tree with bootstrap values '
                           '(not as branch labels) produced by RAxML, '
                           'usually the file name starts with '
                           '"RAxML_bipartitions."')
    main_args.add_argument('-d', '--domain_annotations', help='Path(s) '
                           'to one or more Pfam_scan domain annotation(s) '
                           'of the proteins (standard method to display '
                           'protein domains)', nargs='+', default=[])
    main_args.add_argument('-c', '--config', default=None,
                           help='Path to the geneSearch/calcTree '
                           'configuration file (alternative method to '
                           'display protein domains)')
    add_args.add_argument('-g', '--gffs', default=None,
                           help='Path to one or more GFF files '
                           'containing CDS features of the proteins. '
                           'Used to mark intron positions. Only works '
                           'if the 9th (=last) column of the GFF contains '
                           'the protein IDs', nargs='+')
    add_args.add_argument('-o', '--output_path', default=False,
                          help='Path to the output image file. Must end with '
                          '.PDF, .PNG or .SVG. Tree will be shown in a window '
                          'if omitted')
    add_args.add_argument('-s', '--scale_factor', default=1.0, type=float,
                          help='Horizontal scaling factor of the MSA '
                          '(default: 1.0). '
                          'Decrease this value if your image is too wide!')
    add_args.add_argument('-hl', '--highlight', nargs='*',
                          help='Highlight specific clades or terminal nodes '
                          'with a background color. '
                          'Comma-separated IDs will colour each terminal '
                          'node separately. '
                          'Plus-separated IDs will colour the whole clade '
                          'that contains those IDs. '
                          'Example: --highlight red:prot01,prot02,prot03 '
                          'green:prot04 lightyellow:prot05+prot06')
    add_args.add_argument('-r', '--root', help='Set an outgroup node at which '
                          'the tree will be rooted. If multiple IDs are '
                          'specified (plus-separated), the whole clade that '
                          'contains these IDs will be used for rooting')
    add_args.add_argument('-hn', '--hide_nodes', help='Hide terminal nodes '
                          'that contain the specified substring; e.g. use '
                          '"FBpp" to hide all Drosophila proteins', nargs='+')
    add_args.add_argument('--lineseq', action='store_true',
                          help='Draw a simple line instead of the amino acid / '
                          'nucleotide sequence')

    args = parser.parse_args()

    assert any([args.msa, args.tree, args.domain_annotations]), '''
\nFor visualization, you need to specify at least one of these files:
 - a multiple sequence alignment (--msa)
 - a gene/protein tree (--tree)
 - one or more pfam_scan protein domain annotations (--domain_annotations OR --config)

Read the help (--help) for further information.'''

    main(
         msa_path=args.msa,
         config_path=args.config,
         tree_path=args.tree,
         domain_annotation=args.domain_annotations,
         gff_paths=args.gffs,
         out=args.output_path,
         scale_factor=args.scale_factor,
         highlight=args.highlight,
         root=args.root,
         hide_nodes=args.hide_nodes,
         lineseq=args.lineseq,
    )
