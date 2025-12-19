#!/usr/bin/perl

use warnings;
use strict;

my $gene_census_file = 'gene_presence_absence_transposed.csv';
my $itol_nodes_file = 'itol-tree-node-ids.txt';

open(CENSUS_FILE, "<$gene_census_file") or die $!;

my $header_line = <CENSUS_FILE>;
chomp $header_line;
print "$header_line\n";



my %name2data;

while (<CENSUS_FILE>) {
    chomp;
    my @fields = split("\t");
    my $genome_name = shift @fields;
    $genome_name =~ s/dokuwiki-2023-10-03.versus.//;
     $genome_name =~ s/TIR-proteins.versus.//;
    $genome_name =~ s/.contig.tblastn/_contig/;
    $genome_name =~ s/ATCC_33913_contig/ATCC_33913/;
    $genome_name =~ s/-/_/g;
 
    $name2data{$genome_name} = \@fields;
    warn "Genome name: $genome_name";;
}
close CENSUS_FILE; 

open(ITOL_FILE, "<$itol_nodes_file") or die $!;

while (<ITOL_FILE>) {
    chomp;
    if( defined $name2data{$_} ) {
	print "$_";
	foreach my $field (@{$name2data{$_}}) {
	    print "\t$field";
	}
	print "\n";
    } else {
	warn "Could not find data for $_";
	print "$_\n";
    }
}



close ITOL_FILE;

