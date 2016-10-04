#! /usr/bin/env python

import sqlite3
import sys
from misc import Progress,percentage
import time
import os
import multiprocessing
import pysam

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

DATABASE_PATH = WORK_PATH+'/results/database.sqlite3.COMPLETE.db'
REFERENCE_FASTA_FILE = '/proj/b2014005/SEAseq/ucsc.hg19.fasta'
WINDOW_SIZE = int(1e4)
NUMBER_OF_PARALLEL_PROCESSES = max(len(SAMPLES)*3,multiprocessing.cpu_count()-1) # +1 for the gc calculation
CHUNKS_PER_PROCESS = 10
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


def main():

    total_mapped_reads = {}
    for sample, sample_info in SAMPLES.iteritems():
        bam_file_name = sample_info['bam']
        bamfile = pysam.AlignmentFile(bam_file_name, "rb")
        total_mapped_reads[sample] = bamfile.mapped
        bamfile.close()
    interesting_regions = []
    database = sqlite3.connect(DATABASE_PATH)
    select_string = """SELECT
    regions.region_id,
    chromosomes.chromosome_name,
    regions.start_position, regions.end_position, regions.gc_content,
    """+', '.join([sample+'_regions_info.average_read_depth' for sample in SAMPLES])+"""
    FROM regions
        INNER JOIN chromosomes
            ON chromosomes.chromosome_id = regions.chromosome_id """+' '.join(['INNER JOIN '+sample+'_regions_info ON '+sample+'_regions_info.region_id = regions.region_id' for sample in SAMPLES])
    result_format = ['region_id','chromosome_name','start_position','end_position','gc_content']+[sample for sample in SAMPLES]
    regions = database.cursor().execute(select_string)#.fetchall()
    
    formatted_regions = []
    tmp_count =0
    #for i in xrange(10000): regions.next()
    print '\t'.join(result_format)
    for region in regions:
        
        out_str = ''
        norm_avg_rd = {}
        
        region_dict = {result_format[i]:region[i] for i in xrange(len(region))}
        
        
        formatted_regions.append(
           region_dict
        )
        
        for header in ['region_id','chromosome_name','start_position','end_position','gc_content']: out_str += str(region_dict[header]).replace('.',',')+'\t'

        for sample in SAMPLES:
                norm_avg_rd[sample] = region_dict[sample]/(1e-6*1e-2*total_mapped_reads[sample]) # average read depth normalized per 100M mapped reads
                out_str += str(norm_avg_rd[sample]).replace('.',',')+'\t'
                
        group_average = {}
        for group in COMPARISON_GROUPS.values():
            group_average[group['id']] = sum([norm_avg_rd[sample] for sample in group['samples']])/len(group['samples'])
            group['average_normalized_rd'] = group_average[group['id']]
        
        if group_average[1] == group_average[0]: pass
        elif group_average[1] == 0 or group_average[0] == 0:
            out_str += 'one group zero not other'
            print out_str
            interesting_regions.append(region_dict)
        elif group_average[1]/group_average[0] > 2 or group_average[0]/group_average[1] > 2:
            out_str += 'ojojoj'
            print out_str
            interesting_regions.append(region_dict)
        
        out_str += '\n'
        #print '\t'.join([str(i).replace('.',',') for i in region])
        #if tmp_count > 1000: break
        tmp_count +=1
    database.close()
    
    gc_comp = region_gc_histogram()
    gc_comp.add('all',formatted_regions)    
    gc_comp.add('interesting',interesting_regions)
    print len(interesting_regions),'interesting regions found out of',tmp_count,'regions (',str(percentage(len(interesting_regions),tmp_count)),')'
    gc_comp.summarize()

class region_gc_histogram(object):
    
    def __init__(self,bin_count=10):
        
        self.bins = []
        import sys
        tmp_bin_counter = xrange(bin_count).__iter__()
        step = int(round(100.0 / bin_count,0))
        for i in xrange(0,100,int(step)):
            bin_id = tmp_bin_counter.next()
            self.bins.append( {'id':bin_id,'low':i,'high':min(i+step,100),'samples':{}} )
            #print i,min(i+step,100)
            if i+step > 100: sys.stderr.write('#WARNING: last bin smaller than others.\n')
        self.bins_by_id = { binn['id']:binn for binn in self.bins}
        self.bins_by_low = { binn['low']:binn for binn in self.bins}
        self.bins_by_high = { binn['high']:binn for binn in self.bins}
        self.sample_total_regions = {}
    
    def add(self,name, regions):
        
        import sys
        names = []
        for binn in self.bins: names += binn['samples'].keys()
        if name in set( names ):
            sys.stderr.write('#WARNING: sample '+str(name)+'already in histogram over writing!\n')
        
        for binn in self.bins: binn['samples'][name] = 0
        self.sample_total_regions[name] = 0
        
        for region in regions:
            
            added = False
            for binn in self.bins:
                
                if region['gc_content'] >= binn['low'] and region['gc_content'] < binn['high']:
                    
                    binn['samples'][name] += 1
                    self.sample_total_regions[name] += 1
                    added = True
                    break
    
    def summarize(self, ):
        import sys
        sys.stdout.write( 'Bin\t'+'\t'.join([sample for sample in self.sample_total_regions])+'\n' )
        for binn in self.bins:
            sys.stdout.write( str(binn['low'])+'-'+str(binn['high'])+'\t'+'\t'.join([str(binn['samples'][sample]) for sample in self.sample_total_regions])+'\n' )

if __name__ == '__main__': main()