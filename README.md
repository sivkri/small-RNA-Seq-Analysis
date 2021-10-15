# miRNA-genome mapping using mirdeep2 and mirpro

#!/bin/bash -l

#SBATCH -J micro_map

#SBATCH -N 1 

#SBATCH --ntasks-per-node 4

#SBATCH --mem 4gb

#SBATCH --time=01:00:00

#SBATCH --mail-type=ALL 

module load fastx-toolkit/0.0.14

module load htslib/1.2.1

module load gcc/4.9.4

module load bcftools/1.2

##########################################################

/home/siva/mirPRo.1.1.4/bin/mirpro  -i /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq  -m /home/siva/Documents/microRNA_data/mature_T.fa  -p /home/siva/Documents/microRNA_data/hairpin_T.fa  -d ./mirPro_results -s null -a null -q 0  --gtf /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf  --novel 1 --other null  -g /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  --index /home/siva/Documents/microRNA_data/ath_tair10.idx &> ./mirPro_results/wt_Control_1.txt

#################################################################

bin/mirpro_feature_pro -i /home/siva/Documents/microRNA_data/TRIMMED/wt_Control_1.sorted.out.bam  -t /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf  -o ./miRNA

#################################################################

nohup /home/siva/mirPRo.1.1.4/bin/mirpro  -i /home/siva/Documents/microRNA_data/TRIMMED/wt_Control_1.fastq  -m /home/siva/Documents/microRNA_data/mature_T.fa  -p /home/siva/Documents/microRNA_data/hairpin_T.fa -d ./mirPro_results -s null -a null -q 0  --gtf /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf --novel 1 --other null  -g /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  --index /home/siva/Documents/microRNA_data/ath_tair10.idx >mirpro_scipt.txt


nohup /home/siva/mirPRo.1.1.4/bin/mirpro  -i /home/siva/Documents/microRNA_data/TRIMMED/wt_Control_1.fastq  -m /home/siva/Documents/microRNA_data/mature_T.fa  -p /home/siva/Documents/microRNA_data/hairpin_T.fa -d ./mirPro_results_1 -s null -a null -q 0  --seed 1 --map-detail 1 -v 1 --arm 1 -c 1 --map-len 15 --map-score 60 --5-upstream 3 --3-upstream 3 --5-downstream 3  --3-downstream 3 -n 0 -r 1 --gtf /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf  --novel 1 --other null -g /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  --index /home/siva/Documents/microRNA_data/ath_tair10.idx &>mirpro_scipt_1.txt 


#################################################################

manatee 
-i /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq -o /home/siva/Documents/microRNA_data/Manatee -index /home/siva/Documents/microRNA_data/bowtie_index/arabidopsis_genome -genome /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa  -annotation /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf

################################################################

`exporting oath to bashrc file`


echo 'export PATH=$PATH:/home/siva/perl5/perlbrew/bin' >> ~/.bashrc

echo 'export PATH=$PATH:/home/siva/perl5/perlbrew/perls/perl-5.28.0/bin' >> ~/.bashrc

echo 'export PATH=$PATH:/home/siva/anaconda3/bin' >> ~/.bashrc

echo 'export PATH=$PATH:/home/siva/anaconda3/condabin' >> ~/.bashrc

echo 'export PATH=$PATH:/home/siva/ViennaRNA-2.4.16/bin' >> ~/.bashrc

echo 'export PATH=$PATH:/home/siva/randfold_src' >> ~/.bashrc

echo 'export PATH=$PATH:/home/siva/bowtie' >> ~/.bashrc

###################################################################

**#Step 1: build an index of the genome**

bowtie-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa bowtie_index/arabidopsis_genome

`#step (optional): converting fastq to fasta format`

seqtk seq -a /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq > /home/siva/Documents/microRNA_data/miRNA/FASTA/wt-Control-1.fa

cat /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/siva/Documents/microRNA_data/miRNA/FASTA/wt-Control-1_1.fa

sed -n '1~4s/^@/>/p;2~4p' /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1-3.fastq > /home/siva/Documents/microRNA_data/miRNA/FASTA/wt-Control-1_2.fa

**#Step 2: process reads and map them to the genome**

mapper.pl reads.fa -c -j -k TCGTATGCCGTCTTCTGCTTGT  -l 18 -m -p cel_cluster -s reads_collapsed.fa -t reads_collapsed_vs_genome.arf -v

mapper.pl TRIMMED/wt-Control-1.fastq -e -h -j -l 18 -m -p bowtie_index/arabidopsis_genome -s Mirdeep_miRNA/wt-Control-1_collapsed.fa -t Mirdeep_miRNA/wt-Control-1_collapsed_vs_genome.arf -v 2> Mirdeep_log/mapper_wt-Control-1.txt

**#Step 3: fast quantitation of reads mapping to known miRBase precursors.**

quantifier.pl -p precursors_ref_this_species.fa -m mature_ref_this_species.fa \ -r reads_collapsed.fa -t cel -y 16_19

quantifier.pl -p precursors.fa -m mature.fa -r reads.fa -s star.fa -y now -t cel

quantifier.pl -p hairpin_T.fa -m mature_T.fa -r Mirdeep_miRNA/wt-Control-1_collapsed.fa -y now -d

**#Step 4: identification of known and novel miRNAs in the deep sequencing data:**

miRDeep2.pl reads_collapsed.fa cel_cluster.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa -t C.elegans 2> report.log

$ dos2unix 

miRDeep2.pl Mirdeep_miRNA/wt-Control-1_collapsed.fa Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Mirdeep_miRNA/wt-ABA-3_collapsed_vs_genome.arf -d none none none 2> wt-Control-1_report.log

**#Step 5: browse the results.**

results.html

################################################################################

`Some other useful links to be considered`

https://fatmab-dincaslan.medium.com/mirdeep2-mirna-sequencing-analysis-example-run-by-using-ubuntu-terminal-7922595bb375

https://www.biostars.org/p/213088/

https://www.mdc-berlin.de/content/mirdeep2-documentation

###############################################




