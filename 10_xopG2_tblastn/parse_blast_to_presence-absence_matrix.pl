#!/usr/bin/perl

use strict ;
use warnings ;
use Bio::SearchIO;
use Bio::SeqIO;


my %file2hits;
my %queries;
my %files;

my $coverage_threshold_min = 30; # example $coverage_threshold_min = 90
my $coverage_threshold_max = 100; # example $coverage_threshold_max = 100
my $pc_identical_threshold = 90; # example $pc_identical_threshold = 90

while (my $file = shift) {

    #warn "$file\n";

    $files{$file}++;
    
    my $parser = new Bio::SearchIO(-format => 'blast', -file => $file) ; 
    
    while(my $result = $parser->next_result) {
	my $query_acc = $result->query_accession;
	my $query_desc = $result->query_description;
	my $query_length = $result->query_length;

	$queries{"$query_acc $query_desc"}++;
	
	while(my $hit = $result->next_hit) {
		my $hit_acc = $hit->accession;
	    
	    while(my $hsp = $hit->next_hsp) {
		
		my $query_start = $hsp->start('query');
		my $query_end = $hsp->end('query');
		my $pc_identical = 100 * $hsp->frac_identical; 

		my $hit_string = $hsp->hit_string;		
		if ($hit_string =~ m/\*/) {
			warn "$file $query_acc versus $hit_acc contains STOP codon\n";
		}
		
		my $coverage = 100 * (($query_end - $query_start + 1) / $query_length);

		#warn "$coverage\t$pc_identical\n";
		
		if ($coverage >= $coverage_threshold_min and
                   $coverage <= $coverage_threshold_max and
		    $pc_identical >= $pc_identical_threshold) {
		    
		    $file2hits{$file}{"$query_acc $query_desc"}{$coverage}=$pc_identical;
		   
		}
		
	    }
	}
    }
}

print "File";
print "\t";
print "Coverage (%)";
print "\t";
print "Sequence identity (%)";;
print "\n";


foreach my $file (sort keys %files) {
    foreach my $query (sort keys %queries) {
	if (defined $file2hits{$file}{$query}) {
	    foreach my $coverage (sort {$b=$a} keys %{$file2hits{$file}{$query}}) {
		my $pc_identical = $file2hits{$file}{$query}{$coverage};
		print "$file";
		print "\t";
		print int($coverage);
		print "\t";
		print int($pc_identical);
		print "\n";
	    }
	} else {
	    print "$file";
	    print "\t";
	    print "";
	    print "\t";
	    print "";
	    print "\n";
	}
    }
}
