#!/usr/bin/env perl
use warnings;
use strict ;
use File::Slurp ;
use Getopt::Long ;
use Data::Dumper ;
use Bat1K::Utils ;
use Carp ;


my $genes='';
my $isoforms='' ;
my $outfile='' ;
my $strict ;
my $help ;
my $tmpfile=`mktemp` ;

GetOptions (
	"genes=s" => \$genes,
	"transcripts=s" => \$isoforms,
	"intersect=s" => \$tmpfile,
	"out=s" => \$outfile,
	"strict" => \$strict,
	"help" => \$help
) ;

my $usage="USAGE: $0 --genes(-g) [BED file of predicted genes] --trans(-t) [BED file of potential isoform transcripts]

other optional arguments
	--intersect(-i) | File of overlapping genes and transcripts created with `bedtools intersect -wao -a GENES -b TRANSCRIPTS
	--out(-o)       | Output file names, default is standard out
	--strict(-s)    | Only retain isoforms which have a significantly different exon structure
	                  This removes isofroms which have different CDS start and end but all splice sites are shared.
	                  This should help mitigate incorporating transcripts which only appear different due to misassembly or 5' degredation.
					  Also removes transcripts which are likely targets for NMD (per Colombo 2017), these are also probably misassemblies.
	--help(-h)      | Print this kinda helpful message		   
" ;

if ($help){
	die $usage
}

print STDERR "No genes and/or isoforms provided, will use intersected file $tmpfile\n\n" if (($genes eq "")||($isoforms eq "")) ;

my $fh ;
if ($outfile eq ""){
	$fh = *STDOUT
}
else{
	open $fh, ">$outfile"
}

#my $logfile;


chomp $tmpfile ;
my $call = "bedtools intersect -s -wao -a $genes -b $isoforms > $tmpfile" ;
my @models ;
if ($genes eq ""){
	#$logfile = $tmpfile.".log" ;
	@models=split(/\n/,`cut -f4 $tmpfile|sort -u`) ;
}
else{
	#$logfile = $genes.".log" ;
	system("$call") == 0 || die "Could not intersect $genes and $isoforms" ;
	@models=split(/\n/,`cut -f4 $genes`) ;
}
die "ERROR could not get gene model labels\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
#open LOG, ">$logfile" ;
my @final_isoforms ;
foreach my $gene (@models){
	
	my @good_indexes ;
	my @isoforms ;
	my $hits = `grep -P "\t$gene\t" $tmpfile` ;
	die "ERROR could not grep $gene from $tmpfile\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
	next if (! length($hits)) ;
	my @matches=split(/\n/,$hits) ;
	if (scalar(@matches)>20){
		print STDERR "skipping $gene as more than twenty potential isoforms\n" ;
		next
	}
	
	my ($start,$end,$strand,$exnum,$size,$block,$ichrom) = (split(/\t/,$matches[0]))[1,2,5,9,10,11,12] ;
	next if ($ichrom eq ".") ;
	my @sizes=split(/,/,$size) ;
	my @blocks=split(/,/,$block) ;
	my $gmod = makeCDSHash(split(/\t/,$matches[0])) ;
	#my %gmod ;
	#for (my $i=0;$i<scalar(@sizes);$i++){
	#	my $s = $start + $blocks[$i] ;
	#	my $e = $start + $blocks[$i] + $sizes[$i] ;
	#	$gmod{$i}{left}=$s ;
	#	$gmod{$i}{right}=$e ;
	#}
	foreach my $m (@matches){
		my ($iname,$istrand,$excount)=(split(/\t/,$m))[15,17,21] ;
		next if ($excount == 1) ;
		next if ($strand ne $istrand) ;
		my @matchdeets = split(/\t/,$m) ;
		@matchdeets = splice @matchdeets,12,-1 ;
		#my $nmd = checkNMDTarget(@matchdeets) if ($strict) ;
		if (checkNMDTarget(@matchdeets)){
			print STDERR "$iname, potential isoform of $gene, skipped as likely NMD target\n";
			next
		}
		if ((countFivePrimeUtrExons(@matchdeets)>2)&&($strict)){
			print STDERR "$iname, potential isoform of $gene, skipped as possible bad ORF prediction\n";
			next
		}
		my $iso = makeCDSHash(@matchdeets) ;
		#print Dumper $iso ; exit 0;
		
		$excount = scalar(keys %$iso) -1 ;
		next if ($excount == 1) ;
		my $ovl ;
		if ($strict){
			if(! checkContained($gmod,$iso,1,$strand)){
				my $num = scalar(@isoforms)+1  ;
				$matchdeets[3] = $gene."_x$matchdeets[3]" ;
				push @isoforms, join("\t",@matchdeets) ;
			}
			else{
				next
			}				
		}
		else{
			if(! checkContained($gmod,$iso,0,$strand)) {
				my $num = scalar(@isoforms)+1  ;
				$matchdeets[3] = $gene."_x$matchdeets[3]" ;
				push @isoforms, join("\t",@matchdeets) ;
			}
			else{
				next
			}
		}		
		## Compare all isoforms to each other
		for (my $i=0;$i<scalar(@isoforms);$i++){
			my $first = makeCDSHash(split(/\t/,$isoforms[$i])) ;
			my $id = 0;
			if (scalar(@good_indexes)==0){
				push @good_indexes, $i
			}
			else{
				for (my $j=0;$j<scalar(@good_indexes);$j++){
					my $second = makeCDSHash(split(/\t/,$isoforms[$good_indexes[$j]])) ;
					#$id  = checkIfIdentical($first,$second) ;
					if (checkIfIdentical($first,$second)){
						if (($first->{name} =~ /pb/) && ($second->{name} !~ /pb/)){
							splice @good_indexes, $j, 1 ;
							push @good_indexes, $i ;
						}
						$id = 1 ;
						last
					}
				}
				if ($id == 0){
					push @good_indexes, $i ;
				}
			}
		}
		
		

		
	}
	#print Dumper @good_indexes ;
	#print "\n" ;
	#print Dumper @bad_indexes  ;
	foreach my $c (@good_indexes){
		push @final_isoforms, $isoforms[$c] ;
	}

}

print $fh join ("\n",@final_isoforms)."\n" ;
close $fh unless ($outfile eq "") ;
