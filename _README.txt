Info file for "last rew based" tag version of VOM

rew stands for: representative elementary watershed

Documenation still needs to be updated.

Following changes wrt other tagged versions are included:

- sample output added (VOM-output_samples)

- program for extracting weather data from lsh dataset added (VOM-tools)

- code should be compiler independent

- globally used variables are now in modules.f90

- waterbalance is separated (gasexcwatbal.f90 to watbal.f90 and transpmodel.f90)

- some input variables are separated from input- and code-files to weather-files:
    - air pressure values included in dailyweather.prn (Pres, used for calculation of hourly vd)
    - before air pressure was part (as standard value) of an equation in gasexcwatbal.f90
    
    - molar ratio of CO2 in air included in dailyweather.prn (Ca) and hourlyweather.prn (cah)
    - before it was part of input.par (ca)

- some changes in calculation of waterbalance:
    - see svn commit messages: 53-51, 34-26

- updated start scripts for better performance