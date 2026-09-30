# GoldMethylChicken
Determining the best course of action for comparative species bisulfite methylation profiles of published files for a Golden Eagle's and chicken's heart.  

## Galaxy Pipeline (with coding commands):
The Galaxy Australia pipeline all used default settings.

The SRR/SRA files were imported to Galaxy Australia via sra-tools using the commands below: 
```
mkdir -p ~/.ncbi && cp '/mnt/scratch/job_working_directory/016/039/16039010/configs/tmpjt5h4iyd' ~/.ncbi/user-settings.mkfg &&   echo 'SRR17850309' | sed -r 's/(\,|\;|__cn__)/\n/g' > accessions && for acc in $(cat ./accessions); do ( echo "Downloading accession: $acc..." &&    prefetch -X 200000000 "$acc" && sam-dump --log-level fatal --disable-multithreading  --header --unaligned  "$acc"  | samtools view -Sb - 2> /dev/null > "$acc.bam"  ); done; echo "Done with all accessions."
``` 

The imported BAM files were converted to fastq files through FASTQ Groomer (Galaxy Australia):
``` 
gx-fastq-groomer '/mnt/user-data-volD/data41/a/0/c/dataset_a0c3e7f9-55a3-4a44-8f1e-fe7725e4e860.dat' sanger '/mnt/scratch/job_working_directory/016/039/16039187/outputs/dataset_b3433f64-cca0-48e9-beda-98f2134491d3.dat' sanger ascii summarize_input --fix-id
```

The fastq files were trimmed and quality checked through Trim Galore! (Galaxy Australia):
```
ln -s '/mnt/user-data-volD/data41/b/3/4/dataset_b3433f64-cca0-48e9-beda-98f2134491d3.dat' input_1.fastq &&  trim_galore  --cores ${GALAXY_SLOTS:-4}  --phred33    --output_dir ./      input_1.fastq  --dont_gzip   && if [ -f input_1_trimmed.fq.gz ] ; then mv input_1_trimmed.fq.gz input_1_trimmed.fq ; fi && if [ -f input_1_val_1.fq.gz ] ; then mv input_1_val_1.fq.gz input_1_val_1.fq ; fi && if [ -f input_2_val_2.fq.gz ] ; then mv input_2_val_2.fq.gz input_2_val_2.fq ; fi && if [ -f input_1_unpaired_1.fq.gz ] ; then mv input_1_unpaired_1.fq.gz input_1_unpaired_1.fq ; fi && if [ -f input_2_unpaired_2.fq.gz ] ; then mv input_2_unpaired_2.fq.gz input_2_unpaired_2.fq ; fi && if [ -f input_1.clock_UMI.R1.fq.gz ] ; then mv input_1.clock_UMI.R1.fq.gz input_1.clock_UMI.R1.fq ; fi && if [ -f input_2.clock_UMI.R2.fq.gz ] ; then mv input_2.clock_UMI.R2.fq.gz input_2.clock_UMI.R2.fq ; fi   && ls -lah
```

The 'trusted' reads were mapped to the reference genomes (aquChr2 and galGal6) for bisulfite sequenced reads - via Bismark Mapper (Galaxy Australia):
```
ln -s '/mnt/user-data-volD/data41/1/3/d/dataset_13d4f3e6-5ead-40a6-9863-8b37899a829d.dat' 'input_1.fq' &&  python '/mnt/tools/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/bismark/8c191acde702/bismark/bismark_wrapper.py'  --num-threads "${GALAXY_SLOTS:-4}"   --own-file '/mnt/user-data-volD/data41/7/1/b/dataset_71bfd7db-c60e-4e73-8323-59b1ce4f9b21.dat'   --single-paired 'input_1.fq'  --fastq   --chunkmbs 512    --seed-len 20 --seed-mismatches 0 --seed-extention-attempts 15 --max-reseed 2           --stdout '/mnt/scratch/job_working_directory/016/039/16039223/outputs/dataset_e703c14f-2b85-46f1-a4e3-5f4f3b2d45a9.dat'  --output-report-file '/mnt/scratch/job_working_directory/016/039/16039223/outputs/dataset_b13821b2-180e-4ff1-9f6f-8381c43239eb.dat'   --output '/mnt/scratch/job_working_directory/016/039/16039223/outputs/dataset_3886f30c-38be-4155-9eb6-64b35a664a84.dat'
``` 

The mapped reads were extracted through Bismark Meth. Extractor (Galaxy Australia):
```
python '/mnt/tools/shed_tools/toolshed.g2.bx.psu.edu/repos/bgruening/bismark/8c191acde702/bismark/bismark_methylation_extractor.py'  --multicore "${GALAXY_SLOTS:-4}"  --infile '/mnt/user-data-volD/data41/3/8/8/dataset_3886f30c-38be-4155-9eb6-64b35a664a84.dat'  --single-end    --splitting_report '/mnt/scratch/job_working_directory/016/039/16039711/outputs/dataset_34d49cc3-e2d3-440f-a368-ccd9b3c5c74f.dat'  --mbias_report '/mnt/scratch/job_working_directory/016/039/16039711/outputs/dataset_ab884494-5dcd-436c-9d89-4578f125cdb3.dat'  --cytosine_report '/mnt/scratch/job_working_directory/016/039/16039711/outputs/dataset_e3763842-4e46-4718-ad3a-49596d4111ca.dat' --coverage_file '/mnt/scratch/job_working_directory/016/039/16039711/outputs/dataset_4873b06f-8fc3-4999-8b85-035c7d7e22a4.dat' --genome_file '/mnt/user-data-volD/data41/7/1/b/dataset_71bfd7db-c60e-4e73-8323-59b1ce4f9b21.dat' --cx_context  --comprehensive   --compress '/mnt/scratch/job_working_directory/016/039/16039711/outputs/dataset_f3e8c4a7-615a-449e-9e41-a18c20dd4c5a.dat'
```

Finally, the orthologous genes of the two species' genomes were discovered through protein fasta files processed through OrthoFinder (Galaxy Australia):
```
ln -s '/mnt/user-data-volA/data40/1/d/8/dataset_1d806400-9853-455e-b247-1460583aebfb.dat' 'aquChr2_protein_faa uncompressed.fasta' && ln -s '/mnt/user-data-volA/data40/f/e/3/dataset_fe38ef64-3def-4d13-bc28-15497d18210d.dat' 'galGal6 uncompressed.fasta' &&  orthofinder -f . -S diamond   -I 1.5  -M 'dendroblast'  -t ${GALAXY_SLOTS:-1} -a ${GALAXY_SLOTS:-1} &&  mv OrthoFinder/Results_* results
```

## Methods section 2.2 (Gene methylation profile building):
A python script (csv_maker.py) for the data-manipulation of the methylation files from Galaxy Australia into a comprehensive report with coverage, methylation type, genome position, strand, chromosome, count methylated and count of unmethylated. 

The script from the poly-CpG GitHub by the same GitHub Author to this one was modified for the extraction of gene context for the methylation from the GTF file via a sliding window approach - available in the gene_meth.py file within the github. 

## Methods section 2.3 (Chicken genome cross species comparison):
The methylation profile of the chicken reads and Golden Eagle reads were compared against the galGal6 genome of chickens. The section is divided into multiple scripts to best reflect the repetitive use of scripts throughout the study. 

The script pres_vs_absence.py examines the rates of shared methylation sites among the chicken genome for the 2 species' reads. 

Extraction of the nucleotide kmers was conducted through the kmer_extract.py script by applying a sliding window approach to the chicken genome to isolate the methylation sites and the upstream and downstream 15bp flanking regions to the methylation sites.

Distance matrix construction for the 31-kmer fragments are constructed within Python for the script phylo.py - this python script also includes the hierarchical clustering/ phylogenetic tree construction script. 

The partial-distance Redundancy Based Analysis (dbRDA) based analyses within R are presented within the dbRDA.R script. This script provides a statistical framework for the comparison of kmers between species. 

## Methods section 2.4 (Orthologous based methylation analysis - across reference genomes):
Binning of the genes into shared between species v unique genes (and their associated methylation profiles) was carried out using the OrthoBinning.py script.

The downstream analyses are applied for further investigation: pres_vs_absence.py, kmer_extract.py, phylo.py and dbRDA.R

