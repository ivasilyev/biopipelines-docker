# Some important output files:

* Pipeline launch briefing: `state.json`
* Processed reads: `remove_hg/unmapped`
* Processed full WGS assemblies: `plasmid_merger/*/*.fna`
* Mobile elements: `mgefinder/*/03.results/*/*.tsv`
* Drug resistance: `rgi/*/*.txt`
* Concatenated datasets:
  * `SRST2`: `merge_srst2_results/srst2_concatenated_results.tsv`
  * Reference: `merge_blast_results/blast_concatenated_results.tsv`
* Roary 
  * Raw tree: `roary/out/accessory_binary_genes.fa.newick`
  * Tree with strains only (more compact): `roary/roary_tree_with_strains.newick`
  * Tree with full taxa: `roary/roary_tree_with_taxa.newick`
* Results of reference nucleotide database mapping: `nbee_with_annotation/*/annotated_coverages/*_coverage_annotated_filtered.tsv`
