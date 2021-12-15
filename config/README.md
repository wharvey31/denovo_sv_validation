# Configuration options

A description of every valid option in `config.yaml`.

e.g.

A prefix for outputs.

```yaml
manifest: # Manifest file location (OPTIONAL)
READS: # Dictionary for indexed alignments of raw reads
  ONT: # Example name, can be anything
ASM: # Dictionary for indexed alignments of assemblies
  HIFIASM: # Example name, could be anything
SVPOP: # Dictionary for intersect.tab produced by svpop (https://github.com/EichlerLab/svpop)
  PBSV: # Example name, could be anything
CALLABLE: # Location of callable regions
  PAV: # Bed file locations of valid regions per haplotype in the parents

```
