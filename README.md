# Bat1K

Scripts used to add HQ isoforms providing unique splice information

## TOGA
TOGA (Tool to infer Orthologs from Genome Alignments) was used to provide high quality gene projections which we heavily relied upon in the Bat1K. 

[The code](https://github.com/hillerlab/TOGA) can be found at the H [Hiller lab](https://www.mpi-cbg.de/research-groups/current-groups/michael-hiller/research-focus/) [github](https://github.com/hillerlab/)

## Bat1K Utils.pm
- Perl module with necessary subroutines
- Add this to your perl library path

`mkdir /your/perl/library/path/Bat1K && cp Utils.pm /your/perl/library/path/Bat1K`

## addIsoformsFromEvidence.pl
```
USAGE: addIsoformsFromEvidence.pl --genes(-g) [BED file of predicted genes] --trans(-t) [BED file of potential isoform transcripts]

other optional arguments
	--intersect(-i) | File of overlapping genes and transcripts created with `bedtools intersect -wao -a GENES -b TRANSCRIPTS
	--out(-o)       | Output file names, default is standard out
	--strict(-s)    | Only retain isoforms which have a significantly different exon structure
	                  This removes isofroms which have different CDS start and end but all splice sites are shared.
	                  This should help mitigate incorporating transcripts which only appear different due to misassembly or 5' degredation.
					  Also removes transcripts which are likely targets for NMD (per Colombo 2017), these are also probably misassemblies.
	--help(-h)      | Print this kinda helpful message
 ```
  
 ## updateIsoforms.pl
 ```
 USAGE: updateIsoforms.pl --genes(-g) [BED file of predicted genes] --updates(-u) [BED file of potential mdoels to incorporate]

other optional arguments
	--intersect(-i) | File of overlapping genes and transcripts created with `bedtools intersect -wao -a GENES -b TRANSCRIPTS
	--out(-o)       | Output file name, default is standard out

	--help(-h)      | Print this kinda helpful message
  ```
  
 ### dependencies
 - [Bedtools](https://bedtools.readthedocs.io/en/latest/)
 
 ### perl modules
- File::Slurp
