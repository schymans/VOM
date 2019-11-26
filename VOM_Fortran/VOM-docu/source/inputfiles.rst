Quickstart
===============================


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
The executable "model.x" can be run with these input-files: 

**vom_namelist** 
    Contains all settings to run the VOM.

**pars.txt**
    Contains the (optimized) vegetation parameters (only needed for single run).

**dailyweather.prn**
    Contains the meteorological forcing.


