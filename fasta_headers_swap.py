#!/usr/bin/env python3

__author__ = "Sean P. Jungbluth" and "Robert M. Bowers"
__copyright__ = "Copyright 2018"
__license__ = "GPL 3.0"
__maintainer__ = "Sean P. Jungbluth"
__email__ = "jungbluth.sean@gmail.com"


#import os, sys, argparse, shutil, subprocess, time
import os, sys, argparse, time
from argparse import RawTextHelpFormatter
from os.path import join
from libs.logger import Logger

from libs.versioncontrol import *


try:
    if sys.version_info.major != 3:
        sys.stderr.write('\nError: python version 3 is required - you have python version %d.\n\n' % sys.version_info.major)
        sys.exit(-1)
except Exception:
    sys.stderr.write("Unable to detect python version - assuming python 3 installed.\n\n")
    



def run_fasta_headers_swap(shortnamefasta,keeptaxonomylookup,newlongnamefastafile,outputdir,filenameprefix):
    if not os.path.exists("output-directory"):
        os.makedirs("output-directory")
    logfile = open(str(filenameprefix)+"_program.log", 'w')
    logfile.write("***Start fasta_headers_swap***\n")
    fastaheaderstime = time.time()
    command='fasta_headers_swap.r --shortnamefasta'+str(shortnamefasta)+' --keeptaxonomylookup '+str(keeptaxonomylookup)+' --newlongnamefastafile '+str(newlongnamefastafile)
    print("Full command: \n\n    "+str(command)+"\n")
    process=subprocess.Popen(command, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in process.stdout:
        sys.stdout.write(str(line))
        logfile.write(str(line)[2:-3])
        logfile.write('\n')
    process.wait()
    print("\n\nfasta_headers_swap ran in %s seconds" % round((time.time() - fastaheaderstime),2))
    logfile.write("***End fasta_headers_swap***\n")
    out, err = process.communicate()
    if process.returncode != 0: sys.exit("*Error running fasta_headers_swap*")
    else:
        logfile.close()
#    for f in os.listdir(workingdirectory):
#        if "MaxBin2" in f:
#            shutil.move(f, os.path.join(str(workingdirectory)+"/bin-output_MaxBin2"))


  
if __name__ == '__main__':
    print("""
Wrapper for fasta_headers_swap.r version %s
""" % (fastaheadersswap_version))
    parser = argparse.ArgumentParser(prog='fasta_headers_swap',usage='%(prog)s.py --shortnamefasta [shortnamefasta] --keeptaxonomylookup [keeptaxonomylookup] --newlongnamefastafile [newlongnamefastafile] --outputdir [outputdir] --filenameprefix [filenameprefix] --version', description="""
    description of program"""
    ,formatter_class=RawTextHelpFormatter)
    pathtoimag='~/bin/iMAG/'
    parser.add_argument("--shortnamefasta", dest="shortnamefasta", help="""Description of shortnamefasta""")
    parser.add_argument("--keeptaxonomylookup", dest="keeptaxonomylookup", default='y', help="""Description of keeptaxonomylookup""") #took a guess here that this is a y/n parameter
    parser.add_argument("--newlongnamefastafile", dest="newlongnamefastafile", default="imag-profiler-output", help="""Description of newlongnamefastafile""")
    parser.add_argument("--outputdir", dest="outputdir", default=os.getcwd(), help="""Indicate output directory (default: current working directory)""")
    parser.add_argument("--filenameprefix", dest="filenameprefix", default="fasta-rename_", help="""Indicate output directory (default: current working directory""")
    parser.add_argument('--version', action='version', version='%(prog)s v1.0')
    args = parser.parse_args()
    if len(sys.argv) is None:
        parser.print_help()
    elif len(sys.argv) < 3:
        print('Need more arguements')
        parser.print_help()
    else:
        starttime = time.time()
        workingdirectory=os.getcwd()
        if not os.path.exists(args.outputdir):
            os.makedirs(args.outputdir)
        sys.stdout = Logger(os.path.join(args.outputdir+"/"))
        print('''
***************************************************************************
*                        fasta_header_swap start                          *
***************************************************************************
        ''')
        print("Parameters used to run fasta_headers_swap:")
        print("")
        print("    shortnamefasta: " + str(args.shortnamefasta))
        print("    keeptaxonomylookup: " + str(args.keeptaxonomylookup))
        print("    newlongnamefastafile: " + str(args.newlongnamefastafile))
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



        run_fasta_headers_swap(args.shortnamefasta,args.keeptaxonomylookup,args.newlongnamefastafile,args.outputdir,args.filenameprefix)


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
            
