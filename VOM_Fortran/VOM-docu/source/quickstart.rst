Quickstart
===============================

Dependencies
--------------------------------
Before starting, make sure an appropriate fortran compiler is installed. The code has been tested and developed mainly with gfortran.

In addition, make sure the netcdf-fortran library is installed:

.. code-block:: bash 

    sudo apt install netcdf-bin libnetcdf-dev libnetcdff-dev

The VOM will by default search for the netcdf installation in usr/include. 


Compiling and testing the model
--------------------------------

To compile the model on a unix machine with
the gfortran compiler, type :

.. code-block:: bash 

    make

This will create an executable file called "model.x". Now test the model
by typing:

.. code-block:: bash

    make check

When the test is passed, the model is ready to use. The VOM_namelist contains all settings,
and the model should always be run with the directory of this namelist as working directory. Input and output directories can be defined in the namelist with absolute or relative paths from the working directory. 

Running the model
-----------------
The executable "model.x" can be run : 

.. code-block:: bash

    ./model.x

To run successfully, these input-files are needed:

**vom_namelist** 
    Contains all settings to run the VOM.

**pars.txt**
    Contains the (optimized) vegetation parameters (only needed for single run).

**dailyweather.prn**
    Contains the meteorological forcing.

By default, the executable looks for the vom_namelist in the current workdirectory. The default directory for the other files is /input, relative to the workdirectory. 
This can be changed in vom_namelist, or on the command line:

-i Inputpath to directory with dailyweather.prn, and optionally pars.txt. 

-o Outputpath for all outputfiles.

-n The VOM_namelist (filename can be different)



Model modes
-----------------
The model can be run in 4 different modes, defined by VOM_command in the VOM_namelist:

**1** 	Optimize the model with the Shuffled Complex Evolution algorithm.

**2**   Run without optimization, based on the parameters in pars.txt.

**3**   Run without optimization, based on the parameters in pars.txt. Returns only NCP values as output.

**4** 	Run the model with a set of random parameters.


Other options
-----------------
Initially, the VOM schematized the vegetation as two big leaves, for the perennial and seasonal vegetation. In the newest VOM version, leaf area dynamics can be included in a dynamic way as well, by setting the parameter i_lai_function in the VOM_namelist:

**1** 	No LAI dynamics are included.

**2**   LAI is dynamically modelled, but there is no distinction between shaded and sunlit leaves, as well as different radiation components.

**3**   LAI is dynamically modelled. Shaded and sunlit fractions are determined and a distinction between direct and diffuse radiation is made, but leaves are still treated as a big leave with a single photosynthetic capacity. 

**4**   LAI is dynamically modelled. Shaded and sunlit fractions are determined and a distinction between direct and diffuse radiation is made. The model considers shaded and sunlit leaves with different photosynthetic capacities. 










