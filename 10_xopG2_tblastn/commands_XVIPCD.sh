### Create directory containing genome sequences
mkdir tblastn
cd tblastn
ln -s ../../02_genomes_race-typed/*.contig .

### Create BLAST databases
for i in *.contig; do echo $i; formatdb -pF -i $i ; done

### Perform the TBLASTN searches
for i in *.contig; do echo $i; tblastn -seg no -query ../XVIPCD.faa -db $i -evalue 0.0000000001 -out XVIPCD.versus.$i.tblastn; done

### Come back out of the new directory
cd ..

### Parse TBLASTN files
./extract_flanking_regions.pl ./tblastn/XVIPCD.versus.*.tblastn > XVIPCD_genomic_contexts.txt
perl parse_blast_to_presence-absence_matrix.pl ./tblastn/XVIPCD.versus.*.tblastn > XVIPCD_presence_absence.tsv
