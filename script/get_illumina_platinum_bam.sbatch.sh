#! /bin/bash -l
#SBATCH -A b2013064
#SBATCH -n 1 -p core
#SBATCH -t 8:00:00
#SBATCH -J get_NA12878_platinum_PCRfree
#SBATCH -e /proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation/LOGS/get_illumina_platinum_bam.sbatch.stderr.txt
#SBATCH -o /proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation/LOGS/get_illumina_platinum_bam.sbatch.stdout.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=erik.borgstrom@scilifelab.se

echo "$(date) Running on: $(hostname)"

# from http://www.ebi.ac.uk/ena/data/view/ERS189474
# We have generated paired-end sequence data covering the genome of an CEPH/UTAH individual to a sequence depth of more than 30-fold using the Illumina HiSeq 2000. This individual is a member of the CEPH/UTAH pedigree 1463 (abbreviation: CEPH). The DNA identifier for this individual is NA12878. We obtained the DNA sample NA12878 from The Coriell Institute for Medical Research. Starting with 1ug of DNA, and following random fragmentation, we generated a PCR-Free sequencing library with a median insert size of ~300 bp. 100 base sequence reads were generated from both ends of these templates using the Illumina HiSeq 2000. We carried out purity-filtering (PF) to remove mixed reads, where two or more different template molecules are close enough on the surface of the flow-cell to form a mixed or overlapping cluster. No other filtering of the data has been carried out prior to submission.

path=/proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation/data/NA12878_illumina_PCRfree

mkdir $path -p
cd $path

wget "ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12878_S1.bam"
#     ftp://ftp.sra.ebi.ac.uk/vol1/ERA207/ERA207860/bam/NA12878_S1.bam is alo around though not downloaded currently

echo "$(date) AllDone"


