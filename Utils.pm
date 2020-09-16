package Bat1K::Utils ;
use strict ;
use warnings ;
use Carp ;
use Capture::Tiny ;

use Exporter qw (import) ;

our @EXPORT = qw(makeCDSHash checkContained checkIfIdentical printCDS checkNMDTarget checkPseudo countFivePrimeUtrExons) ;

sub makeCDSHash {
	my @details = @_ ;
	my ($istart,$iname,$istrand,$icdstart,$icdend,$excount,$isize,$iblock)=(@details)[1,3,5,6,7,9,10,11] ;
	my @isizes=split(/,/,$isize) ;
	my @iblocks=split(/,/,$iblock) ;
	my %iso ;
	my $cdsCount = 0 ;
	$iso{name}=$iname ;
	for (my $i=0;$i<scalar(@isizes);$i++){
		my $s = $istart + $iblocks[$i] ;
		my $e = $istart + $iblocks[$i] + $isizes[$i] ;
		if (($e<$icdstart)||($s>$icdend)){
			next
		}
		elsif (($s<=$icdstart)&&($e>=$icdend)){
			$iso{$cdsCount}{left}=$icdstart ;
			$iso{$cdsCount}{right}=$icdend ;
			$cdsCount++ ;			
		}
		elsif(($s<$icdstart)&&($e<$icdend)){
			$iso{$cdsCount}{left}=$icdstart ;
			$iso{$cdsCount}{right}=$e ;
			$cdsCount++ ;				
		}
		elsif(($s>$icdstart)&&($e>$icdend)){
			$iso{$cdsCount}{left}=$s ;
			$iso{$cdsCount}{right}=$icdend ;
			$cdsCount++ ;				
		}			
		else{
			$iso{$cdsCount}{left}=$s ;
			$iso{$cdsCount}{right}=$e ;
			$cdsCount++ ;
			}
	}
	return \%iso ;
}

sub checkContained{
	my $g=shift ;
	my $x=shift ;
	my $str=shift ;
	my $strand=shift ;
	my $size = scalar(keys %$x) ;
	my $gs = scalar(keys %$g) - 1 ;
	my $last = $size -1 ;
	my $check = 0 ;
	my $fivedeg = 0 ;
	for (my $i=0;$i<$last;$i++){
		my $l = $x->{$i}{left} ;
		my $r = $x->{$i}{right} ;
		if ((! $l)||(! $r)){
		my ($stdout, $stderr) = capture {printGene($x)} ;
		confess "Gene CDS has some undefined boundaries\n$stdout"
		}
		foreach (my $k=0;$k<$gs;$k++){
			my $gl = $g -> {$k}{left} ;
			my $gr = $g -> {$k}{right} ;
			if ((! $gl)||(! $gr)){
				my ($stdout, $stderr) = capture {printGene($g)} ;
				confess "Gene CDS has some undefined boundaries\n$stdout"
			}
			if (($i==0)&&($gr eq $r)){
				if ($str>0){ $check++}
				elsif ($gl eq $l){$check++}
				elsif (($gl < $l)&&($strand eq "+")){$check++}
				
				$fivedeg++ if (($strand eq "+")&&($k > 0))
				
			}
			elsif (($i==($last)-1)&&($gl eq $l)){
				if ($str>0){ $check++}
				elsif ($gr eq $r){$check++}
				elsif (($gr > $r)&&($strand eq "-")){$check++}
				
				$fivedeg++ if (($strand eq "-")&&($k < $gs))
			}
			elsif (($gl eq $l) && ($gr eq $r)){
				$check++
			}
		}
	}
	if (($check == $last) && ($check == $gs)){
		return 1 
	}
	elsif (($check == $last)&&($str>0)&&($fivedeg>0)){
		return 1
	}
	elsif ($check == 0){
		return 1
	}
	else {
		return
	}		
}

sub checkIfIdentical{
	my $g=shift ;
	my $x=shift ;
	my $size = scalar(keys %$x) -1 ;
	my $gs = scalar(keys %$g) -1 ;
	my $perfect = 0 ;
	
		for (my $i=0;$i<$size;$i++){
			my $xl = $x->{$i}{left} ;
			my $xr = $x->{$i}{right} ;
			my $gl = $g->{$i}{left} ;
			my $gr = $g->{$i}{right} ;		
			if (($xl == $gl)&&($xr == $gr)){
				$perfect ++
			}
			else {
				last
			}
		}
		if ($perfect == $size){
			return 1
		}
		else {
			return
		}	
}

sub printCDS {
	#no warnings 'uninitialized' ;
	my $g = shift ;
	my $s=scalar(keys %$g) - 1 ;
	print $g->{name}."\n" ;
	for (my $i=0;$i<$s;$i++){
		my $l = $g->{$i}{left} ;
		my $r = $g->{$i}{right} ;
		print "exon$i -> left $l -> right $r\n" ;
	}
	print "\n" 
}

sub checkNMDTarget {
	my @isoform = @_ ;
	my $strand = $isoform[5] ;
	my $chromStart = $isoform[1] ;
	my $introns = 0 ;
	my $gtfifty=  0 ;
	my $utrLen = 0 ;
	if ($strand eq "+"){
		my $stop = $isoform[7] ;
		my @exonStarts = split(/,/,$isoform[11]) ;
		my @exonLens = split(/,/,$isoform[10]) ;
		foreach (my $i=0;$i<scalar(@exonStarts);$i++){
			if (($chromStart+$exonStarts[$i]+$exonLens[$i]) > $stop){
				$introns++ ;
				$gtfifty = 1 if ($utrLen > 50);
				if (($chromStart+$exonStarts[$i])<$stop){
					$utrLen += (($chromStart+$exonStarts[$i]+$exonLens[$i]) - $stop)
				}
				else{
					$utrLen += $exonLens[$i]
				}
			}	
		}
		$introns--
	}
	elsif ($strand eq "-"){
		my $stop = $isoform[6] ;
		my @exonStarts = split(/,/,$isoform[11]) ;
		my @exonLens = split(/,/,$isoform[10]) ;
		foreach (my $i=0;$i<scalar(@exonStarts);$i++){
			if (($chromStart+$exonStarts[$i]) < $stop){
				$introns++ ;
				$gtfifty = 1 if ($utrLen > 50);
				if (($chromStart+$exonStarts[$i]+$exonLens[$i])>$stop){
					$utrLen += ($stop - ($chromStart+$exonStarts[$i]))
				}
				else{
					$utrLen += $exonLens[$i]
				}
			}	
		}
		$introns--
	}
	else{
		confess "$strand not recognised!\n\n" ;
	}
	
	if ($introns > 2){
		return 1
	}
	elsif (($introns ==2 )&&($gtfifty ==1)){
		return 1
	}
	else{
		return ;
	}
}

sub checkPseudo {
	my $seq = shift ;
	#my $pseudo=0 ;
	my @stops = qw(TAG TAA TGA) ;
	my $end = substr($seq,-3) ;
	if (grep /$end/,@stops) {
		$seq = substr($seq,0,-3) ;
	}
	my $orf0=0 ;
	my $orf1=0 ;
	my $orf2=0 ;
	for (my $i=0;$i<length($seq);$i+=3){
		next if ($i+3>length($seq)) ;
		my $codon=substr($seq,$i,3) ;
		if (grep /$codon/,@stops){
			$orf0=1 ;
			last
		}
	}
	for (my $i=1;$i<length($seq)-2;$i+=3){
		#next if ($i+3>length($seq)) ;
		my $codon=substr($seq,$i,3) ;
		if (grep /$codon/,@stops){
			$orf1=1 ;
			last
		}
	}
	for (my $i=2;$i<length($seq)-1;$i+=3){
		#next if ($i+3>length($seq)) ;
		my $codon=substr($seq,$i,3) ;
		if (grep /$codon/,@stops){
			$orf2=1 ;
			last
		}
	}
	
	my $Ns = () =$seq =~ /[Nn]/g ;
	$Ns = $Ns/length($seq) ;
	#return $pseudo ;
	if (($orf0+$orf1+$orf2)==3){
		return 1
	}
	elsif ($Ns > 0.2){
		return 1
	}
	else {return 0}
	
}

sub countFivePrimeUtrExons{
	my @isoform=@_ ;
	my $strand = $isoform[5] ;
	my $chromStart = $isoform[1] ;
	my @exonStarts = split(/,/,$isoform[11]) ;
	my @exonLens = split(/,/,$isoform[10]) ;
	my $count=0 ;
	if ($strand eq "+"){
		my $start=$isoform[6] ;
		foreach (my $i=0;$i<scalar(@exonStarts);$i++){
			$count++ if (($chromStart+$exonStarts[$i]+$exonLens[$i])<$start) ;		
		}
	}
	if ($strand eq "-"){
		my $start=$isoform[7] ;
		foreach (my $i=0;$i<scalar(@exonStarts);$i++){
			$count++ if (($chromStart+$exonStarts[$i])>$start) ;		
		}
	}
	return $count 
}


1 ;
