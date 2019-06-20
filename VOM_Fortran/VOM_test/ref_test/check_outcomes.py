import numpy as np
import pandas as pd
import csv
import sys
from termcolor import colored


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

def check_results(var, var_ref, varname):

    results = var.resample('Y').mean()   
    reference = var_ref.resample('Y').mean()   

    ann_results = np.mean(results)
    ann_reference = np.mean(reference)

    diff = ((ann_results - ann_reference) / ann_reference)*100

    if np.abs(diff) > 10:
        print("Mean annual " + varname + " changed more than 10%:")
        print(diff)
        print("Conclusions are not valid any more")
        print(colored("TEST FAILED", "red"))
        success = False


    if np.abs(diff) < 10:
        print("Mean annual " + varname + " changed less than 10%:")
        print(diff)
        print("Conclusions are still valid")
        print(colored("TEST PASSED", "green"))
        success = True
    return success

succes = np.ones((7), dtype=bool)
succes[0] = check_results(results["etmt"], reference["etmt"],"etmt", )
succes[1] = check_results(results["etmg"], reference["etmg"], "etmg")
succes[2] = check_results(results["esoil"], reference["esoil"], "esoil")
succes[3] = check_results(results["assg"], reference["assg"], "assg")
succes[4] = check_results(results["asst"], reference["asst"], "asst")
succes[5] = check_results(results["jmax25t"], reference["jmax25t"], "jmax25t")
succes[6] = check_results(results["jmax25g"], reference["jmax25g"], "jmax25g")

print("")
print("===================================================")
print("")
print("Result:")
print(str(np.sum(succes)) + "tests out of 6 passed")

if( np.sum(succes) > 3):
    print("More then 3 conclusions still valid")
    sys.exit(0)
else:
    print("Less then 3 conclusions valid")
    print("Throwing an error, good luck next time!")
    sys.exit(1)



