import multiprocessing

WORK_PATH = '/proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation'

SAMPLES = {
    'GM12878_A':{'bam':WORK_PATH+'/data/GM12878.A/P5357_1006/outs/phased_possorted_bam.bam','name':'GM12878_A'},
    'GM12878_B':{'bam':WORK_PATH+'/data/GM12878.B/P5357_1007/outs/phased_possorted_bam.bam','name':'GM12878_B'},
    'NA12878_10x':{'bam':'/proj/b2013064/private/erik/active/160526_10xGenomics_test_runs/results/wgs_data_test_run/NA12878/outs/phased_possorted_bam.bam','name':'NA12878_10x'},
    'NA12878_PCRfree':{'bam':WORK_PATH+'/data/NA12878_illumina_PCRfree/NA12878_S1.bam','name':'NA12878_PCRfree'}
}
COMPARISON_GROUPS = {
    0:{
        'id':0,
        'samples':{
            'GM12878_A':SAMPLES['GM12878_A'],
            'GM12878_B':SAMPLES['GM12878_B'],
            'NA12878_10x':SAMPLES['NA12878_10x']
        }
    },
    1:{
        'id':1,
        'samples':{
            'NA12878_PCRfree':SAMPLES['NA12878_PCRfree']
        }
    }
}

#DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.DELETEME.db'
#DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.db'
#DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.COMPLETE.db'
DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.1kb.db'
# DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.10kb.db'
#DATABASE_PATH = '/scratch/database.sqlite3.1kb.db'
#DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.with_N.db'

REFERENCE_FASTA_FILE = '/proj/b2014005/SEAseq/ucsc.hg19.fasta'
#REFERENCE_FASTA_FILE = '/sw/data/uppnex/reference/Homo_sapiens/hg19/program_files/samtools/concat.fasta'
#REFERENCE_FASTA_FILE = '/sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta'
#shutil.copy2(REFERENCE_FASTA_FILE,'/scratch/reference.fasta')
#shutil.copy2(REFERENCE_FASTA_FILE+'.fai','/scratch/reference.fasta.fai')
#REFERENCE_FASTA_FILE = '/scratch/reference.fasta.fai'

WINDOW_SIZE = int(1e3)
NUMBER_OF_PARALLEL_PROCESSES = max(len(SAMPLES)*3,multiprocessing.cpu_count()-1) # +1 for the gc calculation
CHUNKS_PER_PROCESS = 10