####################################################
# Script extracts gene families according to user-defined groups
# Strain groups are provided in a strain list file
# 
# Input: 1. Strain list file 2. OrthoFinder Output
# Output: Distribution of gene families in core & pangenome.
####################################################
use warnings;

use strict;


my $usage = "\nInsufficient command line arguments provided.\nUSAGE: $0 \n 1. strain list\n 2.  mclOutput\n\n";

my $strainlist = shift(@ARGV) or die ($usage);
my $mcloutput = shift(@ARGV) or die ($usage);

my $Maribacter_621_count = 0;
my $Maribacter_count = 0;
my $Zobellia_count = 0;

my $panAll = 0;
my $coreAll = 0;

my $panMaribacter_621 = 0;
my $coreMaribacter_621 = 0;

my $panMaribacter = 0;
my $coreMaribacter = 0;

my $panZobellia = 0;
my $coreZobellia = 0;

my $panMaribacterAll = 0;
my $coreMaribacterAll = 0;

my $panMaribacterZobellia = 0;
my $coreMaribacterZobellia = 0;

my $panMaribacter621Zobellia = 0;
my $coreMaribacter621Zobellia = 0;


my %Maribacter_621; 
my %Maribacter;
my %Zobellia;
my %paralogCheck;
my %homologs;
my %singletons;

my @line_array;
my @homologs_array;



open(LIST,'<'.$strainlist);

#store all genes in a hash with their corresponding gene family assignment


while (<LIST>){
	
	
		
	chomp($_);
	if ($_ =~ /^Strain/){next};
	
	@line_array = split('\t',$_);
	
	if ($line_array[5] eq "1"){

		$Maribacter_621{$line_array[1]} = "yep";
	}
	
	if ($line_array[5] eq "2"){

		$Maribacter{$line_array[1]} = "yep";
	}
	
	if ($line_array[5] eq "3"){

		$Zobellia{$line_array[1]} = "yep";
	}
}

open(MCL,'<'.$mcloutput);

open(panAll,'>All_pan_genes.txt');
open(coreAll,'>All_core_genes.txt');

open(panMaribacter_621,'>Maribacter_621_pan_genes.txt');
open(coreMaribacter_621,'>Maribacter_621_core_genes.txt');

open(panMaribacter,'>Maribacter_pan_genes.txt');
open(coreMaribacter,'>Maribacter_core_genes.txt');

open(panZobellia,'>Zobellia_pan_genes.txt');
open(coreZobellia,'>Zobellia_core_genes.txt');

open(panMaribacterAll,'>MaribacterAll_pan_genes.txt');
open(coreMaribacterAll,'>MaribacterAll_core_genes.txt');

open(panMaribacterZobellia,'>MaribacterZobellia_pan_genes.txt');
open(coreMaribacterZobellia,'>MaribacterZobellia_core_genes.txt');

open(panMaribacter621Zobellia,'>Maribacter621Zobellia_pan_genes.txt');
open(coreMaribacter621Zobellia,'>Maribacter621Zobellia_core_genes.txt');


#store all genes in a hash with their corresponding gene family assignment


while (<MCL>){
	
	
		
	chomp($_);
	@line_array = split('\s',$_);
	@homologs_array =@line_array[1 .. $#line_array];
	
	foreach (@homologs_array) { 

    	$homologs{$_} = "yes, has homologs";
    	$_=~ /(\S+)_\S+/;
    	my $strain = $1;
    	
    	if (exists $Maribacter_621{$strain}){

    		unless (exists $paralogCheck{$strain}){
  				if ($line_array[0] =~ /OG0001717/)
		    	{
		    	print "Maribacter_621: $strain\n";
				}
        		++$Maribacter_621_count;
        		$paralogCheck{$strain} = "yep";
        	}
         }
    	elsif (exists $Maribacter{$strain}){
    		unless (exists $paralogCheck{$strain}){
				if ($line_array[0] =~ /OG0001717/)
		    	{
		    	print "Maribacter: $strain\n";
				}
        		++$Maribacter_count;
        		$paralogCheck{$strain} = "yep";
        	}
    	}
    	elsif (exists $Zobellia{$strain}){

    		unless (exists $paralogCheck{$strain}){
		    	
		    	if ($line_array[0] =~ /OG0001717/)
		    	{
		    	print "Zobellia: $strain\n";
				}
        		++$Zobellia_count;
        		$paralogCheck{$strain} = "yep";
        	}
		}
	}

	#print "\n";

	#All core or pan genome gene?
	if ($Maribacter_621_count > 0 && $Maribacter_count > 0 && $Zobellia_count > 0){

		print panAll join("\t", @line_array), "\n";
		++$panAll;
		
		if ($Maribacter_621_count == 6 && $Maribacter_count == 9 && $Zobellia_count == 8){
				
			print coreAll join("\t", @line_array), "\n";
			++$coreAll;
		}
	
	}
	
	#Marib621 core or pan genome gene?
	elsif ($Maribacter_621_count > 0 && $Maribacter_count == 0 && $Zobellia_count == 0){
	
		print panMaribacter_621 join("\t", @line_array), "\n";
		++$panMaribacter_621;
		
		if ($Maribacter_621_count == 6 && $Maribacter_count > 0 && $Zobellia_count == 0){
		
			print coreMaribacter_621 join("\t", @line_array), "\n";
			++$coreMaribacter_621;
		}
	
	}

	#Maribacter core or pan genome gene?
	elsif ($Maribacter_621_count == 0 && $Maribacter_count > 0 && $Zobellia_count == 0){
	
		print panMaribacter join("\t", @line_array), "\n";
		++$panMaribacter;
		
		if ($Maribacter_621_count > 0 && $Maribacter_count == 9 && $Zobellia_count == 0){
		
			print coreMaribacter join("\t", @line_array), "\n";
			++$coreMaribacter;
		}
	
	}	


	#Zobellia core or pan genome gene?
	elsif ($Maribacter_621_count == 0 && $Maribacter_count == 0 && $Zobellia_count > 0){
	
		print panZobellia join("\t", @line_array), "\n";
		++$panZobellia;
		
		if ($Maribacter_621_count == 0 && $Maribacter_count == 0 && $Zobellia_count == 8){
		
			print coreZobellia join("\t", @line_array), "\n";
			++$coreZobellia;
		}
	
	}	

	#Maribacter core or pan genome gene?
	elsif ($Maribacter_621_count > 0 && $Maribacter_count > 0  && $Zobellia_count == 0 ){
	
		print panMaribacterAll join("\t", @line_array), "\n";
		++$panMaribacterAll;
		
		if ($Maribacter_621_count == 6 && $Maribacter_count == 9  && $Zobellia_count == 0 ){
		
			print coreMaribacterAll join("\t", @line_array), "\n";
			++$coreMaribacterAll;
		}
	
	}	

	#MaribZobellia core or pan genome gene?
	elsif ($Maribacter_621_count == 0 && $Maribacter_count > 0  && $Zobellia_count > 0 ){
	
		print panMaribacterZobellia join("\t", @line_array), "\n";
		++$panMaribacterZobellia;
		
		if ($Maribacter_621_count == 0 && $Maribacter_count == 9  && $Zobellia_count == 8 ){
		
			print coreMaribacterZobellia join("\t", @line_array), "\n";
			++$coreMaribacterZobellia;
		}
	
	}	
 
	#Marib621Zobellia core or pan genome gene?
	elsif ($Maribacter_621_count > 0 && $Maribacter_count == 0  && $Zobellia_count > 0 ){
	
		print panMaribacter621Zobellia join("\t", @line_array), "\n";
		++$panMaribacter621Zobellia;
		
		if ($Maribacter_621_count == 6 && $Maribacter_count == 0  && $Zobellia_count == 8 ){
		
			print coreMaribacter621Zobellia join("\t", @line_array), "\n";
			++$coreMaribacter621Zobellia;
		}
	
	}	
 
	
	$Maribacter_621_count = 0;
	$Maribacter_count = 0;
	$Zobellia_count = 0;
	%paralogCheck = ();
	 
}

print "

All core genes: $coreAll
Maribacter_621 core genes: $coreMaribacter_621
Maribacter core genes: $coreMaribacter
Zobellia core genes: $coreZobellia
MaribacterAll core genes: $coreMaribacterAll
MaribacterZobellia core genes: $coreMaribacterZobellia
Maribacter621Zobellia core genes: $coreMaribacter621Zobellia

All pan genes: $panAll
Maribacter_621 pan genes: $panMaribacter_621
Maribacter pan genes: $panMaribacter
Zobellia pan genes: $panZobellia
MaribacterAll pan genes: $panMaribacterAll
MaribacterZobellia pan genes: $panMaribacterZobellia
Maribacter621Zobellia pan genes: $panMaribacter621Zobellia
#############################################################################


To merge the outcome _genes.txt with the Kegg annotation list derived from KAAS output:

for i in $(ls -1 *_genes.txt)
do
awk '$0=$2 "\t" $1' ${i} > ${i}_LocusTag.txt
done


rename _genes.txt_LocusTag.txt _LocusTag.txt *

tr -d '\r' < ~/Maribacter_Matti/AllGenes_KEGG_Anno.txt > ~/Maribacter_Matti/AllGenes_KEGG_Anno2.txt

for i in $(ls -1 *_LocusTag.txt)
do
join -a 1 -1 1 -2 1 -t $'\t' <(sort ${i}) <(sort ~/Maribacter/AllGenes_KEGG_Anno2.txt) > ${i}_Anno.txt
done


rename _LocusTag.txt_Anno.txt _Anno.txt *

