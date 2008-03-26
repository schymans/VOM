***************************************************************
Coupled Water Balance and Vegetation Optimality Model

Author: 
Stan Schymanski
School for Environmental Systems Engineering
University of Western Australia
Mailstop M015

Now at:
Max-Planck Institute for Biogeochemistry
PO Box 100164
07701 Jena, Germany
sschym@bgc-jena.mpg.de

July 2006

THIS CODE REPRESENTS THE MODEL DESCRIBED IN THE FOLLOWING PH.D. THESIS:

S. J. Schymanski. Transpiration as the Leak in the Carbon Factory: A Model of Self-Optimising Vegetation. PhD thesis, University of Western Australia, Perth, Australia, 2007.

***************************************************************

COMPILED ON COMPAQ VISUAL FORTRAN 6 FOR WINDOWS
PROGRAM CAN BE STOPPED AT ANY TIME BY PRESSING CTRL+c

The folder 'executables' contains pre-compiled executables. To run the model, copy the executables into a new folder on your computer, together with the input files contained in the folder 'input'.

For graphical display of the optimisation results on a Windows machine, you will need the library aview160.dll. This DLL is not offered under the GNU license and can therefore not be offered as part of the model code. Please Obtain the DLL from Compaq or Intel, or type 'no' when asked by the program if you wish to display the optimisation results.

Model parameters can be set in input6.par and shuffle.par. 

To compute the optimal vegetation properties, make sure the first line in the file 'shuffle.par' says '-optimise', then start the executable from a command line. The optimisation progress will be printed on the screen and saved into the file 'progress.txt'. If you need to interrupt the optimisation at any time, just press CTRL+c. To continue an interrupted optimisation, change the first line in 'shuffle.par' to '-continue' and restart the executable from the command line. The optimisation can take one to several days, depending on your computer, the length of the time series, climate data, soil depth and properties and depending on how quickly the optimisation converges. Usually, the model converges after 3,000 to 6,000 runs. After each run, a new line is displayed on the screen and added to 'progress.txt'. If no new lines appear for a long time, the model has probably hung up.

IF YOU RUN THE EXECUTABLE WITH '-optimise' IN THE FILE 'shuffle.par', ALL PRIOR RESULTS CONTAINED IN THE SAME FOLDER WILL BE OVER-WRITTEN!

Once the optimisation is finished, change the first line in 'shuffle.par' into '-compute' and re-run the executable to compute and write the daily and diurnal fluxes to a file.


---

Copyright (C) 2006  Stan Schymanski

This program is free software: you can redistribute it and/or modify    it under the terms of the GNU General Public License as published by    the Free Software Foundation, either version 3 of the License, or    (at your option) any later version.    This program is distributed in the hope that it will be useful,    but WITHOUT ANY WARRANTY; without even the implied warranty of    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    GNU General Public License for more details.    You should have received a copy of the GNU General Public License    along with this program (e.g. gpl.txt).  If not, see <http://www.gnu.org/licenses/>.





