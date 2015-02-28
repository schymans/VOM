# VOM
The full name of the project is "Coupled Water Balance and Vegetation Optimality Model", in short "VOM". 
The model predicts vegetation water use based on meteorological information, soils and topography only, without the need for prescribing site-specific vegetation properties or calibration against observed fluxes. 
The model is described in the documentation in folder VOM-docu and at:

http://www.bgc-jena.mpg.de/bgc-theory/index.php/Pubs/2009-WRR-SS 

http://theses.library.uwa.edu.au/adt-WU2007.0095/

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
