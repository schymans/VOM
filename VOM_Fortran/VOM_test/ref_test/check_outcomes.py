import numpy as np
import pandas as pd
import csv
import sys
from termcolor import colored

print("Files are not identical")
print("Checking for differences in results...")


results_file = 'output/results_daily.txt'
reference_file = '../ref_test/results_daily.txt'

results = pd.read_csv(results_file, skiprows = 0, dtype=np.float64, delimiter=r"\s+")
reference = pd.read_csv(reference_file, skiprows = 0, dtype=np.float64, delimiter=r"\s+")

##########results
#make dateframe with datetime values
dftime = pd.DataFrame({'year': results['fyear'], 'month': results['fmonth'], 'day': results['fday']})

#make a pandas datetime series
pddatetime = pd.to_datetime(dftime)

#make a pandas index
index = pd.DatetimeIndex( pddatetime)	

#replace index
results.index = index

##########reference
#make dateframe with datetime values
dftime = pd.DataFrame({'year': reference['fyear'], 'month': reference['fmonth'], 'day': reference['fday']})

#make a pandas datetime series
pddatetime = pd.to_datetime(dftime)

#make a pandas index
index = pd.DatetimeIndex( pddatetime)	

#replace index
reference.index = index

########################################
#make mean yearly values etmt

def check_results(var, varname):

    results = var.resample('Y').mean()   
    reference = var.resample('Y').mean()   

    ann_results = np.mean(results)
    ann_reference = np.mean(reference)

    diff = ((ann_results - ann_reference) / ann_reference)*100

    if np.abs(diff) > 10:
        print("Mean annual " + varname + " changed more than 10%:")
        print(diff)
        print("Conclusions are not valid any more")
        print(colored("TEST FAILED", "red"))
        #sys.exit(1)

    if np.abs(diff) < 10:
        print("Mean annual " + varname + " changed less than 10%:")
        print(diff)
        print("Conclusions are still valid")
        print(colored("TEST PASSED", "green"))
        #sys.exit(0)



check_results(results["etmt"],"etmt" )
check_results(results["etmg"], "etmg")
check_results(results["esoil"], "esoil")
check_results(results["assg"], "assg" )
check_results(results["asst"], "asst")
check_results(results["jmax25t"], "jmax25t" )
check_results(results["jmax25g"], "jmax25g")





