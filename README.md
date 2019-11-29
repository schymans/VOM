# VOM
The full name of the project is "Coupled Water Balance and Vegetation Optimality Model", in short "VOM". 
The model predicts vegetation water use based on meteorological information, soils and topography only, without the need for prescribing site-specific vegetation properties or calibration against observed fluxes. 
The model is described in the documentation in folder VOM-docu and at:

http://www.bgc-jena.mpg.de/bgc-theory/index.php/Pubs/2009-WRR-SS 

http://theses.library.uwa.edu.au/adt-WU2007.0095/

## Related papers:
Schymanski, S. J., Roderick, M. L., Sivapalan, M., Hutley, L. B. and Beringer, J.: A test of the optimality approach to modelling canopy properties and CO2 uptake by natural vegetation, Plant, Cell & Environment, 30(12), 1586–1598, [doi:10.1111/j.1365-3040.2007.01728.x](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1365-3040.2007.01728.x), 2007.

Schymanski, S. J., Roderick, M. L., Sivapalan, M., Hutley, L. B. and Beringer, J.: A canopy-scale test of the optimal water use hypothesis, Plant Cell & Environment, 31, 97–111, [doi:10.1111/j.1365-3040.2007.01740.x](https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-3040.2007.01740.x), 2008a.

Schymanski, S. J., Sivapalan, M., Roderick, M. L., Beringer, J. and Hutley, L. B.: An optimality-based model of the coupled soil moisture and root dynamics, Hydrology and Earth System Sciences, 12(3), 913–932, [doi:10.5194/hess-12-913-2008](https://www.hydrol-earth-syst-sci.net/12/913/2008/), 2008b.

Schymanski, S. J., Sivapalan, M., Roderick, M. L., Hutley, L. B. and Beringer, J.: An Optimality-Based Model of the Dynamic Feedbacks between Natural Vegetaton and the Water Balance, Water Resources Research, 45, W01412, [doi:10.1029/2008WR006841](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008WR006841), 2009.

Schymanski, S. J., Roderick, M. L. and Sivapalan, M.: Using an optimality model to understand medium and long-term responses of vegetation water use to elevated atmospheric CO2 concentrations, AoB PLANTS, 7, plv060, [doi:10.1093/aobpla/plv060](https://academic.oup.com/aobpla/article/doi/10.1093/aobpla/plv060/201663), 2015.

## Code description

The initial author of the code and documentation is Stan Schymanski, and initial contributors are Andreas Ostrowski and Steffen Richter.

 The main development line (trunk) of the VOM is stored in the following structure:

    VOM_Fortran
        /vom-code
            The model code written in Fortran 90 (file extension: f90) 
        /vom-scripts
            Different unix scripts to run the model (file extension: script) for:
                /single_computer
                /cluster 
            You have to adapt the scripts (see documentation: How to run the model, and also the description in script files) 
        /vom-docu
            The model documentation is available in html-code, including some text and pdf files
        /vom-input
            input files needed by the model and example forcing data
 
## Documentation

The model documentation can be found at:
https://vom.readthedocs.io/en/latest/

An offline version of the documentation can be found in html-format in the folder VOM-docu/build/html/.
Two additional pdf files containing equations can be found in VOM-docu:
- Equations.pdf
- Watbal3.pdf

The directory VOM-docu/Fordocu contains automatically generated documenation
based on commens in the .f90 files. This one is helpful for understanding the
structure of the source code.








