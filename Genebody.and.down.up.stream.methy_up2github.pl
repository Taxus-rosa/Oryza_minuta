use strict;
use List::Util qw(sum);
use List::Util qw(min);
use POSIX;

my $fn_gene = shift;
my $fn_me = shift;
my $me_type = shift;
my $fn_ave= shift;
my $fn_sgl = shift;

# gene bed list
my $fn_gene="XXX.txt";
my $ini=0; my $ll=6;

## chromosome location and length file
my $fn_gff="BtCt.bed";
my $fn_len="BtCt.len";

## methy.txt
my $fn_me="/.../D1.methylcytosine.txt";
my $me_type="CG";

## output file
my $fn_ave=">XXX.txt";
my $fn_sgl=">XXX_D1_CG.single.txt";
## choose the length of upstream & dowunstream; window size; step size
my $ud=2000; my $window=600; my $step=0.1;

## len file read
my %hs_len; open(LEN,$fn_len) or die $!;
while(<LEN>){
	chomp(); my @a=split(); $hs_len{$a[0]}=$a[1];
} close(LEN); print "$fn_len file done\n";

## gff file read
my %hs_gff; open(GFF,$fn_gff) or die $!;
while(<GFF>){
	chomp(); my @a=split();
	my $up=get_max($a[2]-$ud,1);
	my $dw=get_min($a[3]+$ud,$hs_len{$a[0]});
	$hs_gff{$a[1]}={'c'=>$a[0],'s'=>$a[2],'e'=>$a[3],'u'=>$up,'d'=>$dw};
} close(GFF); print "$fn_gff file done\n";

## main function
open(AVE,$fn_ave) or die $!; open(SG,$fn_sgl) or die $!;
my %hs_gene; my %hs_stat;
for(my $m=$ini;$m<$ll;$m++){
	## bias gene file read
	undef %hs_gene; open(GENE,$fn_gene) or die $!;
	while(<GENE>){
		if($_!~/Om/){next;} # filter title
		chomp(); my @a=split();
		if($a[$m]!~/Om_/){last;} # filter middle line
		$hs_gene{$a[$m+1]}{$a[$m]}=[];
		for(my $i=$hs_gff{$a[$m]}{'s'}-$ud;$i<=$hs_gff{$a[$m]}{'e'}+$ud;$i++){
			if($i<1 or $i>$hs_gff{$a[$m]}{'d'}){push(@{$hs_gene{$a[$m+1]}{$a[$m]}},-1);}
			else{push(@{$hs_gene{$a[$m+1]}{$a[$m]}},0);}
		}
		# print "$a[$m]\t$hs_gff{$a[$m]}{'s'}\t$hs_gff{$a[$m]}{'e'}\t$ud\n@{$hs_gene{$a[$m+1]}{$a[$m]}}\n"; last;
	} close(GENE);
	if(!%hs_gene){next;} # filter middle line
	print "bias gene file $m done\n";

	## Methy file read and set value when overlap
	open(ME,$fn_me) or die $!;
	while(<ME>){
		chomp(); my @a=split();
		my ($chr,$pos,$rnm,$rnn,$type,$p)=($a[0],$a[1],$a[3],$a[4],$a[5],$a[7]);
		if($me_type!~/$type/){next;} 
		foreach my $mk (sort{$a cmp $b} keys %{$hs_gene{$chr}}){
			if($hs_gff{$mk}{'u'}<= $pos and $hs_gff{$mk}{'d'}>=$pos){
				my $index=$pos-$hs_gff{$mk}{'u'};
				## To check whether this site was methylated
				my $x=0; if($p<0.05){$x=$rnm;} 
				my $y=$rnm+$rnn; 
				$hs_gene{$chr}{$mk}[$index]=$x."chu".$y;
			}
		}
	}

	## calculater windows density and print into file
	my $ln=0; undef %hs_stat;
	foreach my $kc (sort{$a cmp $b} keys %hs_gene){
		foreach my $mk (sort{$a cmp $b} keys %{$hs_gene{$kc}}){
			my $s1=$hs_gff{$mk}{'s'}-$ud;
			my $s2=$hs_gff{$mk}{'u'};
			my $s3=$hs_gff{$mk}{'s'};
			my $e3=$hs_gff{$mk}{'e'};
			my $e2=$hs_gff{$mk}{'d'};
			my $e1=$hs_gff{$mk}{'e'}+$ud;
			# print "$mk\t$s1\t$s2\t$s3\t$e3\t$e2\t$e1";
			## To calculate the methylation of upstream
			sl_wd($mk,0,$ud-1,$window,$window*$step,0);
			## To calculate the methylation level of gene regions, first calculate the number of bin based on the fixed 2kb of upstream and downstream regions and the set window step size, and then apply the number of bin to each individual gene body region, because the length of each gene body is not consistent. Therefore, each gene body separately calculates the length of each window according to the number of bin. If the calculated length is a decimal, the corresponding step is an integer added by 1. Therefore, upstream, gene region and downstream can be segmented into the same number of Windows. Take the average after each gene is computed.
			my $st_g=ceil(($e3-$s3+1)/($ud/($window*$step)));
			if($st_g*(($ud-$window)/$window/$step)>=$e3-$s3+1){$st_g-=1;}
			my $wd_g=($e3-$s3+1)-$st_g*($ud-$window)/($window*$step);
			sl_wd($mk,$ud,$e3-$s1,$wd_g,$st_g,5000);
			## To calculate the methylation of downstream
			sl_wd($mk,$e3-$s1+1,$e1-$s1,$window,$window*$step,10000);
			# print "\n";
			$ln++; print "$ln gene done\n";
		}
	}

	foreach my $mk(sort{$a <=> $b}keys %hs_stat){
		my $ave=$hs_stat{$mk}{'s'}/$hs_stat{$mk}{'n'};
		print AVE $ave."\t";
	} print AVE "\n";
} close(OUT); close(AVE);

sub sl_wd{
	my ($id,$is,$ie,$wd,$st,$fg)=@_;
	my $c=$hs_gff{$id}{'c'};
	my @a=@{$hs_gene{$c}{$id}};
	print "$id,$is,$ie,$wd,$st\n";
	print SG "$id";
	for(my $i=$is;$i+$wd-1<$ie+$st;$i+=$st){ 
		my $ine=$i+$wd-1; if($ine>$ie){$ine=$ie;}
		# print "$i,$ine\t@a[$i..$ine]\n";
		my $den=0;
		if(min(@a[$i..$ine])==-1){$den=-1;print "Warning. there is a -1 in $i\n";}
		else{$den=calculater(@a[$i..$ine]);}
		print SG "\t$den";
		my $order=1+($i-$is)/$st;
		if($den>=0){
			$hs_stat{$fg+$order}{'n'}++;
			$hs_stat{$fg+$order}{'s'}+=$den;}
	}
	print SG "\n";
}
sub calculater{
	my @a=@_; my $x=0; my $s=0;
	foreach my $item(@a){
		if($item=~/chu/){
			my @b=split(/chu/,$item); 
			$x+=$b[0]; $s+=$b[1];
		}
	}
	if($x==0 and $s==0){return 0;}
	else{return $x/$s;}
}
sub get_min{
	my $min = 1e15;
	foreach my $item (@_){
		if($item < 0){next;}
		if($item < $min){$min = $item;}
	}
	return $min;
}
sub get_max{
	my $max = -1;
	foreach my $item (@_){
		if($item < 0){next;}
		if($item > $max){$max = $item;}
	}
	return $max;
}
