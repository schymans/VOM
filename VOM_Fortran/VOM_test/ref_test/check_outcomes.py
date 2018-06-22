import numpy as np
import pandas as pd
import csv
import sys


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
#make mean yearly values

etmt_results = results['etmt'].resample('Y').mean()   
etmt_reference = results['etmt'].resample('Y').mean()   

ann_etmt_results = np.mean(etmt_results)
ann_etmt_reference = np.mean(etmt_reference)

d_etmt = ((ann_etmt_results - ann_etmt_reference) / ann_etmt_reference)*100

if np.abs(d_etmt) > 10:
    print("Mean annual etmt changed more than 10%:")
    print(d_etmt)
    print("Conclusions are not valid any more")
    print("TEST FAILED")
    sys.exit(1)

if np.abs(d_etmt) < 10:
    print("Mean annual etmt changed less than 10%:")
    print(d_etmt)
    print("Conclusions are still valid")
    print("TEST PASSED")
    sys.exit(0)


