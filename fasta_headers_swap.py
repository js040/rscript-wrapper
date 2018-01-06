#!/usr/bin/env python3

__author__ = "Sean P. Jungbluth" and "Robert M. Bowers"
__copyright__ = "Copyright 2017"
__license__ = "GPL 3.0"
__maintainer__ = "Sean P. Jungbluth"
__email__ = "jungbluth.sean@gmail.com"

#import os, sys, argparse, shutil, subprocess, shlex, time, re, glob
import os, argparse
from argparse import RawTextHelpFormatter
from os.path import join

from imag.versioncontrol import *


try:
    if sys.version_info.major != 3:
        sys.stderr.write('\nError: python version 3 is required - you have python version %d.\n\n' % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("Unable to detect python version - assuming python 3 installed.\n\n")
    
  
if __name__ == '__main__':
    print("""
Wrapper for rscript fasta_headers_swap
""" % (imag_version))
    parser = argparse.ArgumentParser(prog='fasta_headers_swap',usage='%(prog)s.py -l [pathtolookuptable] -e [bioelement] -b [pathtobins] --version', description="""
    description of program"""
    ,formatter_class=RawTextHelpFormatter)
    pathtoimag='~/bin/iMAG/'
    requiredNamed = parser.add_argument_group('required arguments')
    requiredNamed.add_argument("-l", dest="lookuptable", help="""Indicate path to the bioelement/hmm lookup table. (default: /path/to/iMAG/../doc/bioelement-process-gene-function-lookup-table.txt)""", default=pathtoimag+'doc/bioelement-process-gene-function-lookup-table.txt', required=False)
    parser.add_argument("-e", dest="bioelement", default='All', help="""Indicate the bioelement/biocompounds to be profile. (options: All, Arsenic, C1-compounds, Carbon-fixation, Carbon-monoxide, Halogenated-compounds, Hydrogen, Metals, Methane, Nitriles, Nitrogen, Oxygen, Selenium, Sulfur, Urea) (default: All)""")
    parser.add_argument("-b", dest="pathtobins", default=str(os.getcwd()), help="""Indicate the full path to the bins to be inspected. (default: cwd)""")
    parser.add_argument("-o", dest="outputdir", default="imag-profiler-output", help="""Indicate output directory to be used. (default: imag-profiler-output)""")
    parser.add_argument('--version', action='version', version='%(prog)s v1.0')
    args = parser.parse_args()
    if len(sys.argv) is None:
        parser.print_help()
    elif len(sys.argv) is < 3:
        print('Need more arguements')
        parser.print_help()
   else:
        starttime = time.time()
        workingdirectory=os.getcwd()
        if not os.path.exists(args.outputdir):
            os.makedirs(args.outputdir)
        else:
            print('Error: output directory "'+args.outputdir+'" already exists, select a different output directory!\n')
            print('Exiting...\n')
            sys.exit()
        sys.stdout = Logger(os.path.join(str(workingdirectory)+"/"+args.outputdir+"/"+args.filenameprefix))
        print('''
***************************************************************************
*                        fasta_header_swap start                          *
***************************************************************************
        ''')
        print("Parameters used to run fasta_headers_swap:")
        print("")
        print("    Bioelement(s)/biocompound(s) to profile: " + str(args.bioelement))
        print("    Path to bins: " + str(args.pathtobins))
        print("    Path to bioelement/hmm lookup take: " + str(args.lookuptable))
        print("    Output directory: " + str(args.outputdir))
        print("Software and versions detected by iMAG-profiler:")
        print("")
        print("    %s" % (python_version_found))
        print("    Rscript: v%s" % (rscript_version_found))  ##########
        print("")
#        my_file=Path(workingdirectory+'/'+args.filenameprefix+'_bioelement-table.txt')
#        if my_file.exists():
#            print('Warning: '+args.filenameprefix+'_bioelement-table.txt file already exists! - removing and regenerating')
#            os.remove(''.join(workingdirectory+'/'+args.filenameprefix+'_bioelement-table.txt'))
#            generate_querytable(args.bioelement,args.filenameprefix,workingdirectory,args.lookuptable)
#        else:
#            generate_querytable(args.bioelement,args.filenameprefix,workingdirectory,args.lookuptable)



        print('''
***************************************************************************
*                          fasta_header_swap end                          *
***************************************************************************
        ''')
        print("fasta_header ran in %s seconds" % round((time.time() - starttime),2))
 #       for f in os.listdir(workingdirectory):
 #           if os.path.isdir(f) is True and "hmm-profiler-output" not in f:
 #               shutil.move(f, os.path.join(str(workingdirectory)+"/"+args.outputdir))
 #           elif "bioelement-table" in f:
 #               shutil.move(f, os.path.join(str(workingdirectory)+"/"+args.outputdir))
            
