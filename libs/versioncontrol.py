#!/usr/bin/env python3

__author__ = "Sean P. Jungbluth"
__copyright__ = "Copyright 2017"
__license__ = "GPL 3.0"
__maintainer__ = "Sean P. Jungbluth"
__email__ = "jungbluth.sean@gmail.com"


import subprocess


##lookup python version
proc=subprocess.Popen(['python3 --version &> /dev/stdout'], stdout=subprocess.PIPE, shell=True)
python_version_found=proc.stdout.read()[0:-1]
python_version_found=python_version_found.decode("utf-8")

##lookup rscript version
proc=subprocess.Popen(['Rscript --version &> /dev/stdout | sed \'s/^.*version //\' | sed \'s/ (.*//\''], stdout=subprocess.PIPE, shell=True)
rscript_version_found=proc.stdout.read()[0:-1]
rscript_version_found=rscript_version_found.decode("utf-8")

fastaheadersswap_version='1.0'




