#! /usr/bin/env python

import sqlite3
import sys
import pysam
from misc import Progress
import multiprocessing
import time
import os
import shutil
#
# set workpath and other variables
#
WORK_PATH = '/proj/b2013064/private/ACTIVE/erik/GA_16_14.10X_Chromium_evaluation'
SAMPLES = {
    'GM12878_A':{'bam':WORK_PATH+'/data/GM12878.A/P5357_1006/outs/phased_possorted_bam.bam','name':'GM12878_A'},
    'GM12878_B':{'bam':WORK_PATH+'/data/GM12878.B/P5357_1007/outs/phased_possorted_bam.bam','name':'GM12878_B'},
    'NA12878_10x':{'bam':'/proj/b2013064/private/erik/active/160526_10xGenomics_test_runs/results/wgs_data_test_run/NA12878/outs/phased_possorted_bam.bam','name':'NA12878_10x'},
    'NA12878_PCRfree':{'bam':WORK_PATH+'/data/NA12878_illumina_PCRfree/NA12878_S1.bam','name':'NA12878_PCRfree'}
}
#DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.DELETEME.db'
#DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.db'
DATABASE_PATH = '/scratch/database.sqlite3.db'
#REFERENCE_FASTA_FILE = '/sw/data/uppnex/reference/Homo_sapiens/hg19/program_files/samtools/concat.fasta'
#REFERENCE_FASTA_FILE = '/sw/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta'
REFERENCE_FASTA_FILE = '/proj/b2014005/SEAseq/ucsc.hg19.fasta'
WINDOW_SIZE = int(1e4)
NUMBER_OF_PARALLEL_PROCESSES = max(len(SAMPLES)*3,multiprocessing.cpu_count()-1) # +1 for the gc calculation
CHUNKS_PER_PROCESS = 10
#shutil.copy2(REFERENCE_FASTA_FILE,'/scratch/reference.fasta')
#shutil.copy2(REFERENCE_FASTA_FILE+'.fai','/scratch/reference.fasta.fai')
#REFERENCE_FASTA_FILE = '/scratch/reference.fasta.fai'

def main():
    
    sys.stderr.write('#INFO: starting script.\n')
    sys.stderr.write('#INFO: going to work with '+str(NUMBER_OF_PARALLEL_PROCESSES+1)+' processes in parallel.\n')


    create_regions_table()
    create_chromosomes_table()
    fill_chromosomes_table()
    create_regions()
    
    gc_calculator_process = multiprocessing.Process(target=calculate_region_gc)
    gc_calculator_process.start()
    import time
    time.sleep(20)
    create_per_sample_tables()
    calculate_average_coverage_of_regions()
    import shutil
    shutil.copy2(DATABASE_PATH,WORK_PATH+'/results/database.sqlite3.COMPLETE.db')

def calculate_average_coverage_of_regions():    
    
    sys.stderr.write('#INFO: calculating per sample coverage.\n')
    
    pool = multiprocessing.Pool(processes=NUMBER_OF_PARALLEL_PROCESSES)
    results = pool.imap(imap_func_region_coverage, sample_region_generator(),chunksize=1)
    
    updatevalues = {sample:[] for sample in SAMPLES}
    
    database = sqlite3.connect(DATABASE_PATH)
    tmp_region_count = len(database.cursor().execute('SELECT region_id FROM regions').fetchall())
    database.commit()
    database.close()    
    p = Progress(tmp_region_count*len(SAMPLES), printint=1, unit='windows_added_to_db')
    
    for return_list in results:
        for (sample, region_id, average_read_depth) in return_list:
            
            p.update()
            updatevalues[sample].append( (region_id, average_read_depth) )

            if len(updatevalues[sample]) >= len(return_list):
                updated = False
                while not updated:
                    try:
                        database = sqlite3.connect(DATABASE_PATH)
                        database.cursor().executemany('INSERT INTO '+sample+'_regions_info VALUES (?, ?)',updatevalues[sample])
                        database.commit()
                        database.close()
                        updatevalues[sample] = []
                        updated = True
                    except sqlite3.OperationalError as err:
                        if str(err) == 'database is locked':
                            sys.stderr.write('# WARNING: database was locked trying to write data for sample '+sample+' waiting 1s and trying again...\n');
                            time.sleep(1)
                        else:
                            raise err
    
    for sample in SAMPLES:
        if updatevalues[sample]:
            updated = False
            while not updated:
                try:
                    database = sqlite3.connect(DATABASE_PATH)
                    database.cursor().executemany('INSERT INTO '+sample+'_regions_info VALUES (?, ?)',updatevalues[sample])
                    database.commit()
                    updated = True
                except sqlite3.OperationalError as err:
                    if str(err) == 'database is locked':
                        sys.stderr.write('# WARNING: database was locked trying to write data for sample '+sample+' waiting 1s and trying again...\n');
                        time.sleep(1)
                    else:
                        raise err
    pool.close()
    pool.join()

def create_regions_table():
    """
     Create the table for storing regions that are to be analyzed ie windows over chromosomes
    """
    
    sys.stderr.write('#INFO: creating regions table in database.\n')
    create_regions_table_string = """CREATE TABLE regions (
    region_id int NOT NULL,
    chromosome_id int(3),
    start_position int(9),
    end_position int(9),
    PRIMARY KEY (region_id))"""
    database = sqlite3.connect(DATABASE_PATH)
    try: database.cursor().execute(create_regions_table_string)
    except sqlite3.OperationalError as err:
        if str(err) == 'table regions already exists': sys.stderr.write('#INFO: regions table already in database.\n')
        else:
            raise err
    database.commit()
    database.close()

def create_chromosomes_table():
    """
     Create table to store the chromosome ids names and size
    """
    
    sys.stderr.write('#INFO: creating chromosomes table in database.\n')
    create_chromosomes_table_string = """CREATE TABLE chromosomes (
    chromosome_id int NOT NULL,
    chromosome_name varchar(10),
    size int(9),
    PRIMARY KEY (chromosome_name))"""
    database = sqlite3.connect(DATABASE_PATH)
    try: database.cursor().execute(create_chromosomes_table_string)
    except sqlite3.OperationalError as err:
        if str(err) == 'table chromosomes already exists': sys.stderr.write('#INFO: chromosomes table already in database.\n')
        else: raise err
    database.commit()
    database.close()

def fill_chromosomes_table():
    """
     Fill the chromosomes table get bam file headers for each sample
    """
    database = sqlite3.connect(DATABASE_PATH)
    in_table_ids = database.cursor().execute('SELECT chromosome_id FROM chromosomes').fetchall()
    if not len(in_table_ids):
    
        #
        # look for chromosome names in refernce fasta
        #
        chromosome_names_found = {}
        reference_fasta = pysam.FastaFile(REFERENCE_FASTA_FILE)
        sample = 'referece_fasta'
        for reference_name in reference_fasta.references:
            refernce_length = reference_fasta.get_reference_length(reference_name)
            try:
                chromosome_names_found[reference_name]['samples'].append(sample)
                chromosome_names_found[reference_name]['lengths'].append(refernce_length)
            except KeyError:
                chromosome_names_found[reference_name] = {'samples':[sample], 'lengths':[refernce_length]}

        #
        # look for chromosome names in bam file headders
        #        
        for sample, sample_info in SAMPLES.iteritems():
            bam_file_name = sample_info['bam']
            bamfile = pysam.AlignmentFile(bam_file_name, "rb")
            for reference in bamfile.header['SQ']:
                try:
                    chromosome_names_found[reference['SN']]['samples'].append(sample)
                    chromosome_names_found[reference['SN']]['lengths'].append(reference['LN'])
                except KeyError:
                    chromosome_names_found[reference['SN']] = {'samples':[sample], 'lengths':[reference['LN']]}
            bamfile.close()
        
        #
        # check chromosome names to see if present in all samples and if lengths agree
        #
        for name, info in chromosome_names_found.iteritems():
            if len(info['samples']) < len(SAMPLES)+1:
                sys.stderr.write('#WARNING: skipping reference sequence '+name+' only found in subset of the samples ('+','.join(info['samples'])+').\n')
                chromosome_names_found[name]['pass_filter']=False
            else:
                lengths = set([ str(length) for length in info['lengths'] ])
                if len(lengths) == 1 :
                    sys.stderr.write('#INFO: chromosome '+name+' found in all samples with length '+','.join(lengths)+'.\n')
                    chromosome_names_found[name]['pass_filter']=True
                else:
                    chromosome_names_found[name]['pass_filter']=False
                    sys.stderr.write('#WARNING: chromosome '+name+' found in all samples but has several different lengths '+','.join(lengths)+'.\n')
                #print name,'found in',info['samples'],'with lengths',info['lengths']
        sys.stderr.write('#INFO: '+str(sum([1 for name in chromosome_names_found if chromosome_names_found[name]['pass_filter']]))+' reference sequences out of '+str(sum([1 for name in chromosome_names_found]))+' pass the filters and will be analyzed.\n')
        
        #
        # add the chromosomes that pass filter to the database
        #
        counter = xrange(int(1e9)).__iter__()
        sys.stderr.write('#INFO: filling chromosomes table with data.\n')
        fill_chromosome_table_string = 'INSERT INTO chromosomes VALUES (?, ?, ?)'
        database.cursor().executemany(fill_chromosome_table_string, [    (counter.next(),name, info['lengths'][0]) for name, info in chromosome_names_found.iteritems() if info['pass_filter']    ])
        #database.cursor().executemany(fill_chromosome_table_string, [    (counter.next(),line.rstrip().split('\t')[0], line.rstrip().split('\t')[1]) for line in open(CHR_TSV_FILENAME)          ])
    
    else: sys.stderr.write('#INFO: '+str(len(in_table_ids))+' rows already in chromosomes table.\n')
    database.commit()
    database.close()

def create_regions():
    """
     Create windows of WINDOW_SIZE over the chromosomes and save to regions table
    """
    
    database = sqlite3.connect(DATABASE_PATH)
    if not len(database.cursor().execute('SELECT region_id FROM regions').fetchall()):
        tmp_region_count = 0
        sys.stderr.write('#INFO: creating windows of '+str(WINDOW_SIZE)+' bases.\n')
        counter = xrange(int(1e9)).__iter__()
        for chromosome_id, size in database.cursor().execute('SELECT chromosome_id, size FROM chromosomes'):
            for start in xrange(0,int(size),WINDOW_SIZE):
                if start+WINDOW_SIZE < int(size): end = start+WINDOW_SIZE
                else: end = size
                fill_chromosome_table_string = 'INSERT INTO regions VALUES (?, ?, ?, ?)'
                database.cursor().execute(fill_chromosome_table_string, (counter.next(),chromosome_id,start,end))
                tmp_region_count += 1
        sys.stderr.write('#INFO: '+str(tmp_region_count)+' genomic windows created.\n')
    else:
        tmp_region_count = len(database.cursor().execute('SELECT region_id FROM regions').fetchall())
        sys.stderr.write('#INFO: '+str(len(database.cursor().execute('SELECT region_id FROM regions').fetchall()))+' genomic windows already in database.\n')
    database.commit()
    database.close()
    
def create_per_sample_tables():
    """
     Creating per sample information tables
    """
    
    database = sqlite3.connect(DATABASE_PATH)
    for sample in SAMPLES:
        sys.stderr.write('#INFO: creating '+sample+'_regions_info table.\n')

        create_sample_table_string = """CREATE TABLE """+sample+"""_regions_info (
        region_id int NOT NULL,
        average_read_depth int,
        PRIMARY KEY (region_id))"""

        # try: database.cursor().execute(create_sample_table_string)
        # except sqlite3.OperationalError as err:
        #     if str(err) == 'table '+sample+'_regions_info already exists': sys.stderr.write('#INFO: '+sample+'_regions_info table already in database.\n')
        #     else: raise err
        updated = False
        while not updated:
            try:
                database.cursor().execute(create_sample_table_string)
                updated = True
            except sqlite3.OperationalError as err:
                if str(err) == 'table '+sample+'_regions_info already exists':
                    sys.stderr.write('#INFO: '+sample+'_regions_info table already in database.\n')
                    break
                elif str(err) == 'database is locked':
                    sys.stderr.write('# WARNING: database was locked while trying to write per sample tables, waiting 1s and trying again...\n');
                    time.sleep(1)
                else:
                    raise err

    database.commit()
    database.close()

def sample_region_generator():
    """
     Generator for splitting up the data im smaller chunks to be processed in parallel
    """

    #
    # find regions that are already analyzed
    #
    already_analyzed = {}
    database = sqlite3.connect(DATABASE_PATH)
    for sample in SAMPLES: already_analyzed[sample] = {region_id[0]:True for region_id in database.cursor().execute('SELECT region_id FROM '+sample+'_regions_info').fetchall()}
    database.close()
    already_analyzed_sum = sum([len(tmp) for tmp in already_analyzed.values()])
    
    #
    # calculate the size of chunks to be created
    #
    sys.stderr.write('#INFO: distributing work between processes ...\n')
    import math
    database = sqlite3.connect(DATABASE_PATH)
    tmp_region_count = len(database.cursor().execute('SELECT region_id FROM regions').fetchall())
    database.close()
    tmp_chunk_size = int(math.ceil(  ((tmp_region_count*len(SAMPLES))-already_analyzed_sum) / (CHUNKS_PER_PROCESS*float(NUMBER_OF_PARALLEL_PROCESSES)))  )
    #print 'tmp_chunk_size='+str(tmp_chunk_size)
    sys.stderr.write('#INFO: each process will work with chunks of '+str(tmp_chunk_size)+' windows.\n')
    
    #
    # create chunks of regions that are not analyzed
    #
    sys.stderr.write('#INFO: creating the chunks ...\n')
    chunk = []
    chunks = []
    database = sqlite3.connect(DATABASE_PATH)
    regions = database.cursor().execute('SELECT regions.region_id, chromosomes.chromosome_name, regions.start_position, regions.end_position FROM regions JOIN chromosomes ON chromosomes.chromosome_id = regions.chromosome_id').fetchall()
    database.close()
    for sample_info in SAMPLES.values():
        sys.stderr.write('#INFO: creating chunks for '+str(sample_info['name'])+'...\n')
        sys.stderr.write('#INFO: skipping '+str(len(already_analyzed[sample_info['name']]))+' already analyzed chunks for '+str(sample_info['name'])+'.\n')
        for region_info in regions:
            if region_info[0] not in already_analyzed[sample_info['name']]:
                chunk.append( region_info )
                if len(chunk) == tmp_chunk_size:
                    chunks.append( (sample_info,chunk) )
                    chunk = []
        chunks.append( (sample_info,chunk) )
        chunk = []
    sys.stderr.write('#INFO: '+str(len(chunks))+' chunks created.\n')

    #
    # yield the chunks one by one to the worker processes 
    #
    from random import shuffle
    shuffle(chunks)
    for chunk in chunks: yield chunk

def imap_func_region_coverage(tmp):
    """
     Function to run for each chunk of regions to calculate average coverage of each region
    """

    # imports
    import pysam
    import time
    import os
    
    # extract the info
    sample_info, chunk = tmp
    sample = sample_info['name']
    bam_file_name = sample_info['bam']
    
    if not chunk:
        sys.stderr.write('#WARNING: process '+str(os.getpid())+' got an empty chunk for sample '+sample+' maybe there is no more regions data to analyze.\n')
        return []
    else:
        sys.stderr.write('#INFO: process '+str(os.getpid())+' getting coverage data for sample '+sample+' chunk '+str(chunk[0][0])+'-'+str(chunk[-1][0])+'.\n')
    
    bamfile = pysam.AlignmentFile(bam_file_name, "rb")
    #p_meter = Progress(len(chunk), printint=10, unit=str('windows_analyzed_by_'+str(os.getpid())))
    return_list = []
    for region_info in chunk:
        region_id, chromosome_name, start_position, end_position = region_info
        #average_read_depth = bamfile.count_coverage(reference='chr'+str(chromosome_name), start=end_position-5, end=end_position, read_callback='all')#, quality_threshold=15,)
        average_read_depth = sum([pileup_column.nsegments for pileup_column in bamfile.pileup(str(chromosome_name), start_position, end_position, stepper='all')])/float(WINDOW_SIZE)
        return_list.append((sample, region_id, average_read_depth))
        #p_meter.update()
    bamfile.close()
    
    return return_list 

def calculate_region_gc():
    
    import sys
    from misc import Progress, formatSecods
    import sqlite3
    import os
    import multiprocessing
    import operator
    import time
    
    sys.stderr.write('#INFO: process '+str(os.getpid())+' calculating GC content for regions ...\n')
    
    #
    # check if table column exists and add if needed
    #
    database = sqlite3.connect(DATABASE_PATH)
    columnNames = [col[1] for col in database.cursor().execute('PRAGMA table_info(regions)').fetchall()]
    if 'gc_content' not in columnNames:
        sys.stderr.write('#INFO: process '+str(os.getpid())+' adding column gc_content to regions table in database ...\n')
        database.cursor().execute('alter table regions add column gc_content REAL')
    database.commit()
    database.close()
    
    #
    # get count of regoins to analyze and create progress meter
    #
    database = sqlite3.connect(DATABASE_PATH)
    tmp_region_count_already_analyzed = len(database.cursor().execute('SELECT region_id FROM regions WHERE gc_content IS NOT NULL').fetchall())
    sys.stderr.write('#INFO: '+str(tmp_region_count_already_analyzed)+' regions already analyzed for GC content.\n')
    tmp_region_count = len(database.cursor().execute('SELECT region_id FROM regions WHERE gc_content IS NULL').fetchall())
    p_meter = Progress(tmp_region_count, printint=5, unit='windows_gc_analyzed_by_process_'+str(os.getpid())+'')
    sys.stderr.write('#INFO: process '+str(os.getpid())+' starting to parse '+str(tmp_region_count)+' regions for GC content ...\n')
    regions = database.cursor().execute('SELECT regions.region_id, chromosomes.chromosome_name, regions.start_position, regions.end_position FROM regions JOIN chromosomes ON chromosomes.chromosome_id = regions.chromosome_id WHERE regions.gc_content IS NULL').fetchall()
    database.close()
    sys.stderr.write('#INFO: process '+str(os.getpid())+' splitting regions into chunks for gc calculation.\n')
    region_chunks = [[]]
    index = 0
    region_chunk_size = (len(regions)/(multiprocessing.cpu_count()-len(SAMPLES)*3))/10
    sys.stderr.write('#INFO: process '+str(os.getpid())+' will into use a chunks size of '+str(region_chunk_size)+' regions for gc calculation.\n')
    for region in regions:
        region_chunks[index].append( region )
        if len(region_chunks[index]) >= region_chunk_size:
            region_chunks.append([])
            index += 1
    sys.stderr.write('#INFO: process '+str(os.getpid())+' will work with '+str(len(region_chunks))+' region chunks for gc calculation.\n')
    if len(region_chunks) == 0: return ''
    
    #
    # do work
    #
    sys.stderr.write('#INFO: process '+str(os.getpid())+' starting '+str(multiprocessing.cpu_count()-len(SAMPLES)*3)+' children for region chunks gc calculation.\n')
    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-len(SAMPLES)*3)
    results = pool.imap(imap_func_region_gc, region_chunks, chunksize=1)

    update_values = {}
    last_update_time = time.time()
#    for region_chunk in region_chunks:
#        return_chunk = imap_func_region_gc(region_chunk)
    for return_chunk in results:

        for gc_content, region_id in return_chunk:
            p_meter.update()
            update_values[region_id] = gc_content

            if len(update_values) >= len(return_chunk) and time.time()-last_update_time > 60:#1e3*2:
                sys.stderr.write('#INFO: process '+str(os.getpid())+' '+str(len(update_values))+' regions gc content calculated updating database ....\n')
                update_start_time = time.time()
                update_values = [(gc_content,region_id) for region_id, gc_content in sorted(update_values.iteritems(), key=operator.itemgetter(0))]
                # print update_values[0:100]
                updated = False
                while not updated:
                    try:
                        database = sqlite3.connect(DATABASE_PATH)
                        database.cursor().executemany('UPDATE regions SET gc_content=? WHERE region_id=?',update_values)
                        database.commit()
                        database.close()
                        update_values = {}
                        updated = True
                        last_update_time = time.time()
                    except sqlite3.OperationalError as err:
                        if str(err) == 'database is locked':
                            sys.stderr.write('# WARNING: from process '+str(os.getpid())+' database was locked trying to write data GC content calculation waiting 1s and trying again...\n');
                            time.sleep(1)
                        else:
                            raise err
                sys.stderr.write('#INFO: process '+str(os.getpid())+' database updated in '+formatSecods(int(round(time.time()-update_start_time,0)))+'.\n')

    pool.close()
    pool.join()

    sys.stderr.write('#INFO: process '+str(os.getpid())+' last gc chunk complete, updating database ....\n')
    update_values = [(gc_content,region_id) for region_id, gc_content in sorted(update_values.iteritems(), key=operator.itemgetter(0))]
    updated = False
    while not updated:
        try:
            database = sqlite3.connect(DATABASE_PATH)
            database.cursor().executemany('UPDATE regions SET gc_content=? WHERE region_id=?',update_values)
            database.commit()
            database.close()
            updated = True
        except sqlite3.OperationalError as err:
            if str(err) == 'database is locked':
                sys.stderr.write('# WARNING: from process '+str(os.getpid())+' database was locked trying to write data GC content calculation waiting 1s and trying again...\n');
                time.sleep(1)
            else:
                raise err
    
    sys.stderr.write('#INFO: process '+str(os.getpid())+' GC content calculation complete.\n')

def imap_func_region_gc(region_chunk):
    
    import pysam
    import os
    
    sys.stderr.write('#INFO: process '+str(os.getpid())+' starting to work on gc calculation chunk with regions '+str(region_chunk[0][0])+' to '+str(region_chunk[-1][0])+' ...\n')

    return_chunk = []
    
    reference_fasta = pysam.FastaFile(REFERENCE_FASTA_FILE)
    
    for region_id, chromosome_name, start_position, end_position in  region_chunk:
        
        region_sequence = reference_fasta.fetch(reference=str(chromosome_name),start=start_position,end=end_position)
        gc_content = round(100*(region_sequence.count('G')+region_sequence.count('C')+region_sequence.count('N')*0.5)/float(len(region_sequence)),2)
        return_chunk.append((gc_content,region_id))
        
    reference_fasta.close()
    
    sys.stderr.write('#INFO: process '+str(os.getpid())+' chunk with regions '+str(region_chunk[0][0])+' to '+str(region_chunk[-1][0])+' complete, returing chunk with results to parent.\n')

    return return_chunk
    
#database.crusor().execute('alter table regions add column '+sample+' int')
#database.crusor().execute('alter table regions add column '+sample+' varchar('+str(max([len(sample_name) for sample_name in SAMPLES.keys()]))+')')
#print round(100*float(region_id)/tmp_region_count,2), chromosome_name, start_position, end_position, average_read_depth

if __name__ == '__main__': main()