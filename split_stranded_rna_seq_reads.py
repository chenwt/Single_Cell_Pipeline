#!/usr/bin/env python



#custom wrapper done for RQC 
#to report sense/antisense reads

import os
import sys
import argparse

##updating the library look up path
dir = os.path.dirname(__file__)
sys.path.append(os.path.join(dir,'lib'));
import bamUtils



def main():
    SUB = 'main'
    description  =  """
                    ######
                    Script to split First Strand RNA-Seq reads and report numbers 
                    and store the split reads into sense/antisense files
                    ######
                    """
    #user command line args processing
    parser = argparse.ArgumentParser(description=description,add_help=True)
    parser.add_argument('-f','-files',dest='input_files',default=None,metavar='',help='A comma separated list of bam/sam files')
    parser.add_argument('-cf' ,'-write_split_files',dest='create_files',default=False,
                        action='store_true',help='')
    user_args = parser.parse_args()
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    
    numbers = [ bamUtils.split_RNASeq_reads_bam_by_strand(mapping_file,library_protocol='first-strand',create_split_files=user_args.create_files) for mapping_file in user_args.input_files.split(',') ]
    

if __name__ == "__main__":
    main()
    
