#!/bin/sh
# Script to generate html documentation from comments in .f90 files.
# It presumes that the OS meets the criteria set out at:
# http://sourceforge.net/projects/fordocu/ and that the fordocu code is present
# in ~/Programs/Fordocu/.

cd ../../VOM-code
perl ~/Programs/Fordocu/Fordocu/fordocu.pl -1 -free *.f90 \
-O ../VOM-docu/Fordocu/html -H ../VOM-docu/Fordocu/header 
# Deleting files created in VOM-code directory:
rm -r html
rm fordocu.warn
rm comments.xml

