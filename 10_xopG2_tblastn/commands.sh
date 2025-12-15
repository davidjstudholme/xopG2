### Create directory containing genome sequences
mkdir all_Xcc
cd all_Xcc/
ln -s ../../01_genomes_all_Xcc/*.fasta .

### Create BLAST databases
for i in *.fasta; do echo $i; formatdb -pF -i $i ; done

### Perform the TBLASTN searches
for i in *.fasta; do echo $i; tblastn -seg no -query ../XopG2.faa -db $i -evalue 0.0000000001 -out XopG2.versus.$i.tblastn; done

### Come back out of the new directory                                                                                                                                             
cd ..

### Parse TBLASTN files
./extract_flanking_regions.pl ./all_Xcc/*.tblastn > xopG2_genomic_contexts_all-Xcc-tblastn.txt
perl parse_blast_to_presence-absence_matrix.pl ./all_Xcc/*.tblastn > gene_presence_absence_all-Xcc-tblastn.tsv








