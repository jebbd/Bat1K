#!/usr/bin/env perl
use warnings;
use strict ;
use File::Slurp ;
use Getopt::Long ;
use Data::Dumper ;
use Bat1K::Utils ;

my $genes='';
my $updates='' ;
my $outfile='' ;
my $help ;
my $tmpfile=`mktemp` ;
my $strict ;

GetOptions (
	"genes=s" => \$genes,
	"updates=s" => \$updates,
	"intersect=s" => \$tmpfile,
	"out=s" => \$outfile,
	"help" => \$help,
	"strict" => \$strict 
) ;

my $usage="USAGE: $0 --genes(-g) [BED file of predicted genes] --updates(-u) [BED file of potential mdoels to incorporate]

other optional arguments
	--intersect(-i) | File of overlapping genes and transcripts created with `bedtools intersect -wao -a GENES -b TRANSCRIPTS
	--out(-o)       | Output file name, default is standard out

	--help(-h)      | Print this kinda helpful message		   
" ;

if ($help){
	die $usage
}

print STDERR "No genes and/or isoforms provided, will use intersected file $tmpfile\n\n" if (($genes eq "")||($updates eq "")) ;

my $fh ;
if ($outfile eq ""){
	$fh = *STDOUT
}
else{
	open $fh, ">$outfile"
}

chomp $tmpfile ;
my $bedtools = "bedtools intersect -s -wao -a $genes -b $updates > $tmpfile" ;
system("$bedtools")==0 || die "could not intersect files, $genes and $updates\n\n" ;
my $getnames = `cut -f4 $tmpfile | awk -F"_x" '{print \$1}' |sort -u` ;
die "Could not retrieve gene identifiers from input bed, $genes\n\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
my @names = split(/\n/,$getnames) ;
my @final_updates ;

foreach my $n (@names){
	#print "collecting all isoforms of gene, $n\n" ;
	my $pattern = 'grep -P \'\t'.$n.'[_\t]\' '.$tmpfile ;
	#print "$pattern\n" ;
	my $hits = `$pattern` ;
	die "ERROR could not perform grep on $tmpfile\n" if ($? != 0 || ${^CHILD_ERROR_NATIVE} != 0) ;
	next if (! length($hits)) ;
	my @matches=split(/\n/,$hits) ;
	#print Dumper \@matches ;
	
	### make two array of the "already" in place isoforms and the potential "updates"
	my @already ;
	my @a_ids ;
	my @updates ;
	my @u_ids ;
	foreach my $line (@matches){
		my @deets = split(/\t/, $line) ;
		next if ($deets[5] ne $deets[17]) ;
		my @adeets = (@deets)[0,1,2,3,4,5,6,7,8,9,10,11] ;
		if (!grep (/\b$deets[3]\b/, @a_ids)){
			push @a_ids, $deets[3] ;
			push @already, join("\t",@adeets) ;	
		}
		#print "checking ".join("\t",@adeets)." is new gene model\n" ;
		my @udeets = (@deets)[12,13,14,15,16,17,18,19,20,21,22,23] unless ($deets[12] eq ".");
		if (checkNMDTarget(@udeets)){
			print STDERR $deets[15].", potential isoform of ".$deets[3].", skipped as likely NMD target\n";
			next
		}
		if ((countFivePrimeUtrExons(@udeets)>2)&&($strict)){
			print STDERR $deets[15].", potential isoform of ".$deets[3].", skipped as likely bad ORF prediction\n" ;
			next
		}
		#print "checking ".join("\t",@udeets)." is new potential update\n" ;

		
		if (!grep (/$deets[15]/, @u_ids)){
			push @u_ids, $deets[15] ;
			push @updates, join("\t",@udeets) ;
		}		
		
	}
	
	foreach my $u (@updates){
		my $already_in = 0 ;
		#my $ovl ;
		my $uHash = makeCDSHash(split(/\t/,$u)) ;
		my $uname = (split(/\t/,$u))[3] ;
		#print "Checking if $uname is an acceptable isoform\n" ;
		my $strand = (split(/\t/,$u))[5] ;
		foreach my $a (@already){
			my $aHash = makeCDSHash(split(/\t/,$a)) ;
			#$ovl = checkContained($aHash,$uHash,1,$strand) ;
			if ($strict){
				if (checkContained($aHash,$uHash,1,$strand)){
					$already_in = 1 ;
					last
				}
			}
			elsif(checkIfIdentical($aHash,$uHash)){
					$already_in = 1 ;
					last
			}
			
		}
		if ($already_in < 1){
			$u =~ s/($uname)/$n\_x$1/ ;
			push @already, $u 
		}
		else{
			print STDERR "$uname, a potential isoform of $n, is already included\n"
		}
	}
	push @final_updates, @already ;
}

print $fh join("\n",@final_updates)."\n" ;
