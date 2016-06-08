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
usage: showTree [-h] [-m MSA] [-t TREE] [-c CONFIG] [-d DOMAIN_ANNOTATION]
                [-o OUTPUT_PATH] [-s SCALE_FACTOR]
                [-hl [HIGHLIGHT [HIGHLIGHT ...]]] [-r ROOT] [-n]
                [-hn HIDE_NODES [HIDE_NODES ...]]

optional arguments:
  -h, --help            show this help message and exit

main arguments:
  -m MSA, --msa MSA     Path to an untrimmed multiple sequence alignment in
                        FASTA format (the one produced by calcTree usually
                        ends with "_aln.fa")
  -t TREE, --tree TREE  Path to the best-scoring RAxML tree with support
                        values (not as branch labels) produced by calcTree,
                        usually named "RAxML_bipartitions.(...)_aln.tree"
  -c CONFIG, --config CONFIG
                        Path to the geneSearch/calcTree configuration file
                        (not required if you specify a domain annotation with
                        "--domain_annotation")
  -d DOMAIN_ANNOTATION, --domain_annotation DOMAIN_ANNOTATION
                        Path to a Pfam_scan domain annotation of the proteins
                        (not required if all domain annotations are in the
                        geneSearch config)

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
  -n, --no_alignment    Do not draw the multiple sequence alignment or domains
  -hn HIDE_NODES [HIDE_NODES ...], --hide_nodes HIDE_NODES [HIDE_NODES ...]
                        Hide terminal nodes that contain the specified text.
```

Example output
------------

An example output tree that was generated with the command
`showTree -m ORTHO_ALL1421_aln.fa -c geneSearchConfig.txt -t RAxML_bipartitions.ORTHO_ALL1421_aln.tree -s 0.5 -o tree.pdf`

![example_tree](http://i.imgur.com/k52BxkR.png)