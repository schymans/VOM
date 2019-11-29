#!/bin/sh
# Script to generate html documentation from comments in .f90 files.
# It presumes that the OS meets the criteria set out at:
# http://sourceforge.net/projects/fordocu/ and that the fordocu is 
# installed in the folder you specify when running the script, e.g.
# in a terminal: $./generate_fordocu.sh ~/Programs/Fordocu/Fordocu/

fordocu_dir=$1
cd ../../VOM-code
perl $fordocu_dir/fordocu.pl -1 -free *.f90 \
-O ../VOM-docu/Fordocu/html -H ../VOM-docu/Fordocu/header 
# Deleting files created in VOM-code directory:
rm -r html
rm fordocu.warn
rm comments.xml

