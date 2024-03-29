# mirRNA-Seq-Analysis
# miRNA-genome mapping using mirdeep2 and mirpro

This repository contains a sample script for performing microRNA analysis using the software package mirPRo and miRDeep2. The script is specifically designed for the dataset available at [https://doi.org/10.3390/ijms22137153](https://doi.org/10.3390/ijms22137153).

### Usage

To run the analysis, follow these steps:

1. Ensure that the required modules are loaded by executing the following commands:
   ```bash
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
   ```

2. Prepare the necessary input files and directories:
   - Place the trimmed FASTQ file at `/home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq`.
   - Prepare the reference files: `mature_T.fa`, `hairpin_T.fa`, `Arabidopsis_thaliana.TAIR10.42.gtf`, and `Arabidopsis_thaliana.TAIR10.dna.toplevel.fa` at their respective locations.

3. Execute the following commands to run the analysis:
```bash
# Step 1: Run mirpro
/home/siva/mirPRo.1.1.4/bin/mirpro -i /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq -m /home/siva/Documents/microRNA_data/mature_T.fa -p /home/siva/Documents/microRNA_data/hairpin_T.fa -d ./mirPro_results -s null -a null -q 0 --gtf /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf --novel 1 --other null -g /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --index /home/siva/Documents/microRNA_data/ath_tair10.idx &> ./mirPro_results/wt_Control_1.txt

# Step 2: Run mirpro_feature_pro
bin/mirpro_feature_pro -i /home/siva/Documents/microRNA_data/TRIMMED/wt_Control_1.sorted.out.bam -t /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf -o ./miRNA

# Step 3: Run mirpro with additional options
nohup /home/siva/mirPRo.1.1.4/bin/mirpro -i /home/siva/Documents/microRNA_data/TRIMMED/wt_Control_1.fastq -m /home/siva/Documents/microRNA_data/mature_T.fa -p /home/siva/Documents/microRNA_data/hairpin_T.fa -d ./mirPro_results_1 -s null -a null -q 0 --seed 1 --map-detail 1 -v 1 --arm 1 -c 1 --map-len 15 --map-score 60 --5-upstream 3 --3-upstream 3 --5-downstream 3 --3-downstream 3 -n 0 -r 1 --gtf /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf --novel 1 --other null -g /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa --index /home/siva/Documents/microRNA_data/ath_tair10.idx &>mirpro_scipt_1.txt

# Step 4: Run manatee
manatee -i /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq -o /home/siva/Documents/microRNA_data/Manatee -index /home/siva/Documents/microRNA_data/bowtie_index/arabidopsis_genome -genome /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -annotation /home/siva/Documents/microRNA_data/Arabidopsis_thaliana.TAIR10.42.gtf
 ```

4. Additional steps and useful links:
   - Exporting paths to `bashrc` file:
     ```bash
     echo 'export PATH=$PATH:/home/siva/perl5/perlbrew/bin' >> ~/.bashrc
     echo 'export PATH=$PATH:/home/siva/perl5/perlbrew/perls/perl-5.28.0/bin' >> ~/.bashrc
     echo 'export PATH=$PATH:/home/siva/anaconda3/bin' >> ~/.bashrc
     echo 'export PATH=$PATH:/home/siva/anaconda3/condabin' >> ~/.bashrc
     echo 'export PATH=$PATH:/home/siva/ViennaRNA-2.4.16/bin' >> ~/.bashrc
     echo 'export PATH=$PATH:/home/siva/randfold_src' >> ~/.bashrc
     echo 'export PATH=$PATH:/home/siva/bowtie' >> ~/.bashrc
     ```

   - Building the genome index:
     ```bash
     bowtie-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa bowtie_index/arabidopsis_genome
     ```

   - Converting FASTQ to FASTA format (optional):
     ```bash
     seqtk seq -a /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq > /home/siva/Documents/microRNA_data/miRNA/FASTA/wt-Control-1.fa
     cat /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > /home/siva/Documents/microRNA_data/miRNA/FASTA/wt-Control-1_1.fa
     sed -n '14s/^@/>/p;24p' /home/siva/Documents/microRNA_data/TRIMMED/wt-Control-1-3.fastq > /home/siva/Documents/microRNA_data/miRNA/FASTA/wt-Control-1_2.fa
     ```

   - Additional mirDeep2 commands:
     ```bash
     miRDeep2.pl reads_collapsed.fa cel_cluster.fa reads_collapsed_vs_genome.arf mature_ref_this_species.fa mature_ref_other_species.fa precursors_ref_this_species.fa -t C.elegans 2> report.log
     $ dos2unix 
     miRDeep2.pl Mirdeep_miRNA/wt-Control-1_collapsed.fa Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Mirdeep_miRNA/wt-ABA-3_collapsed_vs_genome.arf -d none none none 2> wt-Control-1_report.log
     ```

   - Browsing the results: Open `results.html` in a web browser.

### Additional Resources
`Some other useful links to be considered`

For more information on the usage of mirPRo and miRDeep2, please refer to the following resources:

- [Mirdeep2 miRNA Sequencing Analysis Example](https://fatmab-dincaslan.medium.com/mirdeep2-mirna-sequencing-analysis-example-run-by-using-ubuntu-terminal-7922595bb375)
- [miRDeep2 Documentation](https://www.mdc-berlin.de/content/mirdeep2-documentation)
- [BioStars: Mirdeep2 Tutorial](https://www.biostars.org/p/213088/)

Please note that some output files have been excluded from this repository due to GitHub's 25 MB limit. If you require those files, please refer to the original source mentioned in the DOI link above.

