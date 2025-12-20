### Create directory containing genome sequences
mkdir race-typed_Xcc
cd race-typed_Xcc/
ln -s ../../02_genomes_race-typed/*.contig .

### Create BLAST databases
for i in *.contig; do echo $i; formatdb -pF -i $i ; done

### Perform the TBLASTN searches
for i in *.contig; do echo $i; tblastn -seg no -query ../XopG2.faa -db $i -evalue 0.0000000001 -out XopG2.versus.$i.tblastn; done

### Come back out of the new directory
cd ..

### Parse TBLASTN files
./extract_flanking_regions.pl ./race-typed_Xcc/*.tblastn > xopG2_genomic_contexts_race-typed_Xcc_tblastn.txt
perl parse_blast_to_presence-absence_matrix.pl ./race-typed_Xcc/*.tblastn > gene_presence_absence_race-typed_Xcc_tblastn.tsv
