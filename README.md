showTree, a visualization tool for gene/protein families.
=========

showTree can visualize the following files in a single figure:
- a multiple sequence alignment (MSA)
- a gene/protein tree
- a protein domain annotation


The inputs above can be combined freely, e.g. you may provide any combination 
such as
- tree + MSA + domains
- tree + MSA
- tree + domains
- MSA + domains
- ...



Requirements
------------

The input files (MSA and gene tree) can be computed with [calcTree](https://ebbgit.uni-muenster.de/ckeme_01/geneSearch/wikis/calcTree).
showTree can also parse the [calcTree](https://ebbgit.uni-muenster.de/ckeme_01/geneSearch/wikis/calcTree)
config file to automatically find the domain annotation of proteins in the MSA.
Alternatively, you can provide your own (pfam_scan) domain annotation.


showTree requires Python2.7 or Python3.4+ and the following Python packages:
- ete3

`pip install 'ete3==3.0.0b35'`
- Biopython

`pip install biopython`


Usage
------------

```
usage: showTree.py [-h] [-m MSA] [-t TREE]
                   [-d DOMAIN_ANNOTATION [DOMAIN_ANNOTATION ...]] [-c CONFIG]
                   [-o OUTPUT_PATH] [-s SCALE_FACTOR]
                   [-hl [HIGHLIGHT [HIGHLIGHT ...]]] [-r ROOT]
                   [-hn HIDE_NODES [HIDE_NODES ...]]

optional arguments:
  -h, --help            show this help message and exit

main arguments:
  -m MSA, --msa MSA     Path to an untrimmed multiple sequence alignment in
                        FASTA format (the one produced by calcTree usually
                        ends with "_aln.fa")
  -t TREE, --tree TREE  Path to the best-scoring gene tree with bootstrap
                        values (not as branch labels) produced by RAxML,
                        usually the file name starts with"RAxML_bipartitions."
  -d DOMAIN_ANNOTATION [DOMAIN_ANNOTATION ...], --domain_annotation DOMAIN_ANNOTATION [DOMAIN_ANNOTATION ...]
                        Path(s) to one or more Pfam_scan domain annotation(s)
                        of the proteins (standard method to display protein
                        domains)
  -c CONFIG, --config CONFIG
                        Path to the geneSearch/calcTree configuration file
                        (alternative method to display protein domains)

additional arguments:
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
  -hn HIDE_NODES [HIDE_NODES ...], --hide_nodes HIDE_NODES [HIDE_NODES ...]
                        Hide terminal nodes that contain the specified
                        substring; e.g. use "FBpp" to hide all Drosophila
                        proteins
```

Example output
------------

An example output tree that was generated with the command
`showTree -m ORTHO_ALL1421_aln.fa -c geneSearchConfig.txt -t RAxML_bipartitions.ORTHO_ALL1421_aln.tree -s 0.5 -o tree.pdf`

![example_tree](http://i.imgur.com/k52BxkR.png)
