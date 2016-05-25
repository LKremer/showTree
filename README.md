showTree
=========

A visualization tool that for protein families.
Displays a gene tree with the corresponding MSA and protein domains next to it.


Requirements
------------

showTree requires Python2.7 or Python3.4+ and the following Python packages:
- ete3

`pip install 'ete3==3.0.0b35'`
- Biopython

`pip install biopython`


Usage
------------

```
usage: showTree [-h] -c CONFIG -m MSA -t TREE [-o OUTPUT_PATH]
                [-s SCALE_FACTOR] [-hl [HIGHLIGHT [HIGHLIGHT ...]]] [-r ROOT]
                [-n] [-hn HIDE_NODES [HIDE_NODES ...]]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT_PATH, --output_path OUTPUT_PATH
                        Path to the output image file (.PDF). Tree will be
                        shown in a window if omitted
  -s SCALE_FACTOR, --scale_factor SCALE_FACTOR
                        Horizontal scaling factor of the MSA (default: 1.0).
                        Decrease this value if your image is too wide!
  -hl [HIGHLIGHT [HIGHLIGHT ...]], --highlight [HIGHLIGHT [HIGHLIGHT ...]]
                        Highlight specific clades or terminal nodes with a
                        background color. Comma-separated IDs will colour each
                        terminal node separately. Plus-separated IDs will
                        colour the whole clade that contains those IDs.
                        Example: --highlight red:prot01,prot02,prot03
                        green:prot04 lightyellow:prot05+prot06
  -r ROOT, --root ROOT  Set an outgroup node at which the tree will be rooted.
                        If multiple IDs are specified (plus-separated), the
                        whole clade that contains these IDs will be used for
                        rooting
  -n, --no_alignment    Do not draw the multiple sequence alignment or domains
  -hn HIDE_NODES [HIDE_NODES ...], --hide_nodes HIDE_NODES [HIDE_NODES ...]
                        Hide terminal nodes that contain the specified text.

required arguments:
  -c CONFIG, --config CONFIG
                        Path to the geneSearch/calcTree configuration file
  -m MSA, --msa MSA     Path to the untrimmed T-Coffee multiple sequence
                        alignment produced by calcTree, usually ends with
                        "_aln.fa"
  -t TREE, --tree TREE  Path to the best-scoring RAxML tree with support
                        values (not as branch labels) produced by calcTree,
                        usually named "RAxML_bipartitions.(...)_aln.tree"
```