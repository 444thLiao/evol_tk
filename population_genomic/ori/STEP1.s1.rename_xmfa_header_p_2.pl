#!/usr/bin/perl -w

open IN,"<","parsnp.xmfa.renamed" or die $!;
open OUT,">","mauve_out.xmfa.renamed";
while(<IN>)
{
	chomp;
	if(/^>/)
	{
		print OUT "$_\n";
	}
	else
	{		
		$seq=$_;
		$seq=~tr/atgc/ATGC/;
		print OUT "$seq\n";
	}
}
close IN;
close OUT;
