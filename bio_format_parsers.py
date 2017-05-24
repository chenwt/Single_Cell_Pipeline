#!/usr/bin/env python


import sys
import pprint
import math
import matplotlib.pyplot as plt
import pybedtools 
import numpy as np
import itertools




class gene(object):
    
    def __init__(self):
        self.gene_id = None
        self.transcripts = {}
        self.start = None
        self.end = None
        self.strand = None
        self.chromosome = None
        self.attributes = None
        self.score = 0
        self.sequence = None

    def print_gene_info(self):
        print self.gene_id,self.chromosome,self.strand,self.start,self.end,len(self.transcripts)

    def print_transcript_info(self):
        for transcript in self.transcripts:
            print transcript,self.transcripts[transcript].exon_starts,self.transcripts[transcript].exon_ends
    
    def num_transcripts(self):
        return len(self.transcripts)
    
    def num_exons(self):
        return sum( [ len(self.transcripts[transcript].exon_starts) for transcript in self.transcripts ] )
    
    def length(self):
        pass
    
    def get_start_stop(self):
        pass
    
    
    def get_bed_entry(self,transcripts=None,**kwargs):
        
        bed_lines = []
        if transcripts:
            for transcript_id in self.transcripts:
                transcript_obj = self.transcripts[transcript_id]
                transcript_start = transcript_obj.exon_starts[0]
                transcript_end = transcript_obj.exon_ends[-1]
                transcript_block_size = [str(exon_end-exon_start) for exon_start,exon_end in itertools.izip(transcript_obj.exon_starts,transcript_obj.exon_ends)]
                transcript_block_start = [str(exon_start - transcript_start) for exon_start in transcript_obj.exon_starts]
                
                bed12_format_line_transcript_block_size = ','.join(transcript_block_size)
                bed12_format_line_transcript_block_start = ','.join(transcript_block_start)
                
                transcript_name ='%s:%s' % ( self.gene_id,transcript_id)
                
                #bed12 format
                bed_line = [self.chromosome, transcript_start, transcript_end, transcript_name, self.score, self.strand,
                            transcript_start , transcript_start, kwargs.get('rgb',0), len(transcript_obj.exon_ends),
                            bed12_format_line_transcript_block_size, bed12_format_line_transcript_block_start]
                
                #convert all the entries to string inside the list of a transcript
                #and append to the main list for storing entries of transcripts for a gene
                bed_lines.append([str(x) for x in bed_line ])
            
            return bed_lines

        else:
            bed_line = [self.chromosome,self.start,self.end,self.gene_id,self.score,self.strand]
            #convert all the entries to string inside the list
            bed_lines.append( [str(x) for x in bed_line ] ) 
            return bed_lines
        
        
    
        
class transcript(object):
    
        def __init__(self):
            self.exon_starts = []
            self.exon_ends = []
            self.sequence = None
            
        def num_exon(self):
            pass
        
        def debug(self):
            print self.transcript_id,self.exon_start,self.exon_end
            
    
    
    
def summarize_gtf_annotation_file(gtf_file,summarize_by=None,outFile=None):
    '''
    given a gtf file mainly produced by cuffmerge this method can
    summarize it by gene or transcript level
    **PS: assumption is that for each gene the exons are reported, sorted by position**
    variable : summarize_by can take the value of gene,transcript
    '''
    
    genes = []
    
    #last seen gene and transcript
    #start with a raw initialized gene and transcript
    last_gene = gene()
    last_transcript = transcript()
    
    
    for count,gtf_line in enumerate(pybedtools.BedTool(gtf_file)):
        current_gene_id = gtf_line.attrs.get('gene_id','NA')
        current_transcript_id = gtf_line.attrs.get('transcript_id','NA')
        feature_name = gtf_line.fields[2].strip()  #gtf_line.name wont work as it takes the value from the gene_name attrs
        
        if feature_name != "exon":
            print gtf_line.name
            print 'method only works on the exon feature'
            continue
        
        
        #new gene found
        if ( current_gene_id != last_gene.gene_id ):
            #store the last gene information first
            if last_gene.gene_id is not None: 
                genes.append(last_gene)
            
            #updating the information for new gene seen after last gene info is saved
            last_gene = gene() #new gene
            last_gene.gene_id = current_gene_id
            last_gene.start = gtf_line.start
            last_gene.end = gtf_line.end
            last_gene.strand = gtf_line.strand
            last_gene.chromosome = gtf_line.chrom
            last_gene.transcripts[current_transcript_id] = transcript()
            last_gene.attributes = gtf_line.attrs
            
            #new gene will have new transcripts too
            new_transcript = last_gene.transcripts[current_transcript_id]
            new_transcript.exon_starts.append(gtf_line.start)
            new_transcript.exon_ends.append(gtf_line.end)
            
        #else if we see the same gene again
        else:
            last_gene.end   = gtf_line.end
            last_gene.attributes = gtf_line.attrs
            
            #existing transcript of a gene
            if (last_gene.transcripts.get(current_transcript_id,None) is not None):
                current_transcript = last_gene.transcripts.get(current_transcript_id,None)
                current_transcript.exon_starts.append(gtf_line.start)
                current_transcript.exon_ends.append(gtf_line.end)
            
            #new isoform of an existing gene
            else:
                last_gene.transcripts[current_transcript_id] = transcript()
                new_isoform = last_gene.transcripts[current_transcript_id]
                new_isoform.exon_starts.append(gtf_line.start)
                new_isoform.exon_ends.append(gtf_line.end)
            
            
    #appending the last gene which will not be compared in the loop 
    genes.append(last_gene)
    
    return genes

    '''
    print "number of genes ",len(genes)
    
    for g in genes:
        print g.print_gene_info()
        print g.print_transcript_info()
    '''
    
    
   

if __name__ == "__main__":
    
    
    #generate gene/transcript structure from a exon annotation file
    #mainly useful for cufflinks generated gtf but could be further extended and made more general
    genes = summarize_gtf_annotation_file(sys.argv[1])
    
    
    #do something with the assembled genes structure
    total_transcripts = 0 
    total_exons = 0
    exons_per_transcript = []
    gene_lengths         = []
    bed_entries_genes = []
    bed_entries_transcripts = []
    
    for g in genes:
        #create a list of bed entries
        bed_entries_genes.extend(g.get_bed_entry())
        bed_entries_transcripts.extend(g.get_bed_entry(transcripts=True))
        total_transcripts += g.num_transcripts()
        total_exons       += g.num_exons()
        gene_lengths.append(g.end - g.start)

        #num of exon in each transcript of a gene
        exons_per_transcript.extend ( [  len(g.transcripts[trans].exon_starts) for trans in g.transcripts ] )
    
     
        
#    print math.log(exons_per_transcript,2)
    plt.hist( [ math.log(x,2) for x in exons_per_transcript], bins=50)
    plt.grid(True)
    plt.xlabel('log2(exons/transcript)')
    plt.ylabel('frequency')
    plt.legend()
    plt.savefig('test_exons_per_transcript.png',dpi=175)
    
    
    print 'Total number of genes : %d ' % len(genes)
    print 'Total number of transcripts : %d ' % total_transcripts
    print 'Total number of exons : %d ' % total_exons
    print 'Transcript length %0.2f +/- %0.2f' % (np.mean(gene_lengths),np.std(gene_lengths))
    
    
     #this is mainly needed as pybedtools.create_interval_from_list cant create a BedTool object using list of list
    #but it can take an iterator so need to create a generator first in order to convert list of bed entries to bedtool object 
    def helper_make_Bedtool_intervaltype_obj(bed_entries):
        for bed_entry in bed_entries:
            yield pybedtools.create_interval_from_list(bed_entry)
        
    genes_bed = pybedtools.BedTool(helper_make_Bedtool_intervaltype_obj(bed_entries_genes))
    transcripts_bed = pybedtools.BedTool(helper_make_Bedtool_intervaltype_obj(bed_entries_transcripts))
    
    genes_bed.saveas('genes.bed')
    transcripts_bed.saveas('transcripts.bed')

    
    
    
    
def split_bed_by_strand(bedFile):
    """
    given a bed file split the file into two
    by strand +/-
    and return the two bed file paths
    """
    
    features = pybedtools.BedTool(bedFile)
    
    features_sense_strand_file = pybedtools.BedTool([feature for feature in features if feature.strand == '+' ]).saveas().fn
    features_antisense_strand_file = pybedtools.BedTool([feature for feature in features if feature.strand == '-' ]).saveas().fn
    
    return(features_sense_strand_file,features_antisense_strand_file)
    
    
    







##obselete
## using the the pybedtools parser
def parse_gtf_file(gtf_file):
    """
    given a gtf file as input this function will parse and return a single line
    is a generator : yield single parsed gtf entry one at a time
    """
    
    #sample_gtf_line (split in 3 lines for display purposes)
    #chromosome_1    Cufflinks       exon    20356   20687   .       +       .       
    #gene_id "XLOC_000036"; transcript_id "TCONS_00000050"; exon_number "1"; gene_name "g3"; oId 
    #"PAC:26903463"; contained_in "TCONS_00000051"; nearest_ref "PAC:26903463"; class_code "="; tss_id "TSS41"; p_id "P3";

    ##just creating an abtract class/container to hold the values for features of a gtf line
    class gtf_line_container(object):
            
            def __init__(self):
                self.chr=None
                self.start=None
                self.end=None
                self.transcript_id=None
                self.gene_id=None
                self.feature_type=None
                self.attributes=None
        
    for num,line in enumerate(open(gtf_file,'r')):
        
        line_split = line.strip().split('\t')
        
        #instantiate the container
        gtf_line = gtf_line_container()
        
        attributes_split=line_split[8].split(';')
        #sample attribute split line
        #['gene_id "XLOC_000001"', ' transcript_id "TCONS_00000001"', ' exon_number "1"', ' oId "CUFF.31.1"', ' class_code "u"', ' tss_id "TSS1"', '']
        
        #gene_id
        gene_id=attributes_split[0].split(' ')[1].replace('"','').strip()
        
        #transcrpit_id
        transcript_id=attributes_split[1].split(' ')[2].replace('"','').strip()
        
        
        gtf_line.chr=line_split[0].strip()
        gtf_line.start=line_split[3].strip()
        gtf_line.end=line_split[4].strip()
        gtf_line.strand=line_split[6].strip()
        gtf_line.feature_type=line_split[2].strip()
        gtf_line.gene_id = gene_id
        gtf_line.transcript_id = transcript_id
        gtf_line.attributes = line_split[8].strip()
        yield gtf_line



    
    

#custom method
#posted here for FYI purposes

"""


def create_pandas_df_from_gtf_file(merged_gtf_file):
    '''
    given a merged.gtf file produced by cuffmerge this method can
    summarize it by gene or transcript level
    **PS: assumption is that for each gene the exons are reported, sorted by position**
    variable : summarize_by can take the value of gene,transcript
    '''
    
    features = []
    for count,line in enumerate(pybedtools.BedTool(merged_gtf_file)):
        #parse following from the attrs
        attributes = ['gene_id','transcript_id','exon_number','gene_name','class_code','tss_id','p_id']
        attr_value = [ line.attrs.get(attr,'NA').strip() for attr in attributes ]
        feature = [ x.strip() for x in line[0:8]] + attr_value
        features.append(feature)
    
    #creating pandas df
    colnames = ['chr','source','feature','start','stop','score','strand','frame'] + attributes
    out_pandas_df = pandas.DataFrame.from_records(features,columns=colnames)
    
    return out_pandas_df




def create_gtf_entry_from_df(group):
    output = []
    for index,row in group.T.iteritems():
        gtf_line = [row['chr'], row['source'], row['feature'], 
                    str(row['start']), str(row['stop']), row['score'],
                    row['strand'], str(row['frame'])]
        
        gtf_line_attrs = ['gene_id','transcript_id','exon_number','gene_name','tss_id','p_id']
        attr_line = ''
        
        for attr in gtf_line_attrs:
            value=str(row.get(attr,'NA'))
            if value.strip() != 'NA':
                attr_line += ('%s "%s"; ' % (attr,value))
        
        gtf_line.append(attr_line)
        output.append(gtf_line)
        
    return (output)
   
"""
   
  #interesting peice of logic
            #depending on the user entered summarize_by varibale either the transcript_id or gene_id is used
#            this_feature_id=eval('gtf_line.%s' % feature_name)
