#!/usr/bin/perl

use strict ;
use warnings ;
use Bio::SearchIO;
use Bio::SeqIO;
#use Bio::DB::EUtilities;
use Bio::DB::GenBank;

my %file2hits;
my %files;
my %coverages;
my %stop_codons;

my $upstream_flank = 2000;
my $downstream_flank = 2000;

my $coverage_threshold_min = 0.40; # example $coverage_threshold_min = 0.90 
my $coverage_threshold_max = 1.00; # example $coverage_threshold_max = 1.00

#my $coverage_threshold_min = 0.90; # example $coverage_threshold_min = 0.90
#my $coverage_threshold_max = 1.00; # example $coverage_threshold_max = 1.00

my $pc_identical_threshold = 95; # example $pc_identical_threshold = 90

while (my $file = shift) {
    
    #warn "$file\n";
    
    $files{$file}++;
    
    my $parser = new Bio::SearchIO(-format => 'blast', -file => $file) ; 
    
    while(my $result = $parser->next_result) {
	my $query_acc = $result->query_accession;
	my $query_desc = $result->query_description;
	my $query_length = $result->query_length;
	
	while(my $hit = $result->next_hit) {
	    my $hit_acc = $hit->accession;
	    
	    while(my $hsp = $hit->next_hsp) {
		
		my $query_start = $hsp->start('query');
		my $query_end = $hsp->end('query');
		my $hit_start = $hsp->start('hit');
		my $hit_end = $hsp->end('hit');
		
		my $pc_identical = 100 * $hsp->frac_identical; 

		### Keep a tally of stop codons
		my $hit_string = $hsp->hit_string;		
		if ($hit_string =~ m/\*/) {
		    warn "$file $query_acc versus $hit_acc contains STOP codon\n";
		    $stop_codons{$file}{$hit_acc}++;
		}

		### Calculate coverage over the query sequence
		my $coverage = ($query_end - $query_start + 1) / $query_length;
		my $pc_coverage = int($coverage * 100);
		$coverages{$file}{$hit_acc}{$pc_coverage}++;
		
		if ($coverage >= $coverage_threshold_min and
		    $coverage <= $coverage_threshold_max and
		    $pc_identical >= $pc_identical_threshold) {
		    
		    $file2hits{$file}{$hit_acc} = [$hit_start, $hit_end];
		}
	    }
	}
    }
}


### List genomes that have stop codons
foreach my $file (sort keys %stop_codons) {
    warn "$file has STOP codon\n";
}


### List genomes that show complete coverage
my %coverage_to_printline;
print "Coverages:\n";
foreach my $file (sort keys %coverages){
    my $fasta_filename;
    if ($file =~ m/\.versus\.(.*.fasta)\.tblastn/) {
	$fasta_filename = $1;
    } else {
	die "Could not parse FASTA filename from TBLASTN file '$file'\n";
    }
    foreach my $acc (keys %{$coverages{$file}}) {
	my @coverages = sort {$b<=>$a} keys %{$coverages{$file}{$acc}};
	my $highest_coverage = shift @coverages;
	my $printline = "";
	$printline .= "$fasta_filename";
	$printline .= "\t";
	$printline .= "$acc";
	$printline .= "\t";
	$printline .= "$highest_coverage %";
	foreach my $i (@coverages) {
	    $printline .= ", $i %";
	}
	$printline .= "\t";
	if (defined $stop_codons{$file}{$acc}) {
	    $printline .= "STOP codon ";
	} else {
	    $printline .= "";
	}
	push @{$coverage_to_printline{$highest_coverage}}, $printline;
    }
}

foreach my $coverage (sort {$b<=>$a} keys %coverage_to_printline) {
    foreach my $readline (sort @{$coverage_to_printline{$coverage}}) {
	print "$readline\n";
	warn "$readline\n";
    }
}

### Iterate over each BLAST hit
foreach my $file (sort keys %file2hits) {
    my $fasta_filename;  
    if ($file =~ m/\.versus\.(.*.fasta)\.tblastn/) {
	$fasta_filename = $1;
    } else {
	die "Could not parse FASTA filename from TBLASTN file '$file'\n";
    }

    warn "Examining $fasta_filename\n";
    
    foreach my $hit_acc (sort keys %{$file2hits{$file}}) {
	if (defined $file2hits{$file}{$hit_acc}) {
	    my ($hit_start, $hit_end) = @{$file2hits{$file}{$hit_acc}};

	    ### Print header line for this genome
	    print "\n$hit_acc: $hit_start .. $hit_end in genome $fasta_filename\n";
	    if (defined $stop_codons{$file}{$hit_acc}) {
		print "STOP codon in $hit_acc: $hit_start .. $hit_end\n";
		warn "STOP codon in $hit_acc: $hit_start .. $hit_end\n";
	    }

	    ### Print list of coverages
	    my @coverages = @{$coverages{$file}{$hit_acc}};
	    foreach my $coverage (sort {$b<=>$a} @coverages) {
		if ($coverage == 100) {
		} else {
		    print "!!! ";
		}
		print "Hit covers $coverage %\n";
	    }
	    
	    ### Print a blank line
	    print "\n"; 

	    
	    ### Extract features in the flanking sequence
	    my $X = $hit_start - $upstream_flank;
	    my $Y = $hit_end + $downstream_flank;
	    	    
	    my $gb = Bio::DB::GenBank->new(
		-email => 'd.j.studholme@exeter.ac.uk'
		);
	    
	    my $seq = $gb->get_Seq_by_acc($hit_acc);
	    
	    for my $feature ($seq->get_SeqFeatures) {
		
		my $start = $feature->location->start;
		my $end   = $feature->location->end;
		
		# overlap test
		next if $end   < $X;
		next if $start > $Y;

		my $type = $feature->primary_tag;
		next unless
		    $type eq 'CDS' ||
		    is_mobile_element($feature);
		
		my $name = cds_name($feature);   # from earlier, works for CDS
		my $loc  = $feature->location->to_FTstring;
		
		my $product = get_qual($feature, 'product') // 'NA';
		my $note    = get_qual($feature, 'note')    // 'NA';
		my $mobtype = get_qual($feature, 'mobile_element_type') // 'NA';
		my $highlight = "";
		$highlight = "***" if is_mobile_element($feature);
		
		print join("\t",
			   $type,
			   $loc,
			   $name,
			   $product,
			   $mobtype,
			   $highlight,
		    ), "\n";
	    }
	    
	    ### Throttle requests
	    sleep 0.5;
	    
	}
    }
}
    

sub get_qual {
    my ($feat, $tag) = @_;
    return undef unless $feat->has_tag($tag);
    my @vals = $feat->get_tag_values($tag);
    return join(',', @vals);
}

sub cds_name {
    my ($feat) = @_;

    # Priority order: most specific → most generic
    for my $tag (qw(gene locus_tag protein_id product)) {
	if ($feat->has_tag($tag)) {
	    my @vals = $feat->get_tag_values($tag);
	    return join(',', @vals);
	}
    }
    return 'unnamed_CDS';
}


sub is_mobile_element {
    my ($feat) = @_;

    # 1. Feature type
    return 1 if $feat->primary_tag eq 'mobile_element';

    # 2. Known repeat container types
    return 1 if $feat->primary_tag eq 'repeat_region';

    # 3. Look for mobility-related qualifiers
    for my $tag ($feat->get_all_tags) {
	for my $val ($feat->get_tag_values($tag)) {
	                return 1 if $val =~ /\b(
                transposon|
                insertion\s+sequence|
                integrase|
                transposase|
                IS\d+|
                retrotransposon|
                prophage
            )\b/ix;
	}
    }

    return 0;
}
