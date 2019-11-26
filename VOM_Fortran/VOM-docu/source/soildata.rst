
Soil data
===============================


Soil layers 
-------------------------------

Soil layer thickness is globally defined with the parameter i_delz, or in the file soilprofile.par for each soil layer separately. The sum of the soil layers need to be equal to the topographical parameter i_cz (maximum elevation). The VOM will raise a warning and correct i_cz if this is not the case. The river bed level i_zr needs to align with the soil layers. Also here, the VOM will give a warning and correct i_zr if needed. 

Soil properties 
-------------------------------
Soil water retention and hydraulic conductivity parameters for 12 major soil textural groups according to [1]_. Some units of the parameters (alpha and Ksat) are converted for the use in the VOM.


+-----------------+---------+------------+--------------+------+--------------+
|Texture          | theta-r |	theta-s  | alpha (1/m)  | n    | K-sat (m/s)  |
+-----------------+---------+------------+--------------+------+--------------+
|Sand             | 0.045   |	0.43 	 | 14.5         | 2.68 | 8.25*10-5    |
+-----------------+---------+------------+--------------+------+--------------+
|Loamy Sand       | 0.057   |   0.41     | 12.4         | 2.28 | 4.053*10-5   |
+-----------------+---------+------------+--------------+------+--------------+
|Sandy Loam       | 0.065   |   0.41     | 7.5          | 1.89 | 1.228*10-5   |
+-----------------+---------+------------+--------------+------+--------------+
|Loam             | 0.078   |   0.43     | 3.6          | 1.56 | 2.889*10-6   |
+-----------------+---------+------------+--------------+------+--------------+
|Silt             | 0.034   |   0.46     | 1.6          | 1.37 | 6.944*10-7   |
+-----------------+---------+------------+--------------+------+--------------+
|Silt Loam        | 0.067   |   0.45     | 2.0          | 1.41 | 1.25*10-6    |
+-----------------+---------+------------+--------------+------+--------------+
|Sandy Clay Loam  | 0.100   |   0.39     | 5.9          | 1.48 | 3.639*10-6   |
+-----------------+---------+------------+--------------+------+--------------+
|Clay Loam        | 0.095   |   0.41     | 1.9          | 1.31 | 7.222*10-7   |
+-----------------+---------+------------+--------------+------+--------------+
|Silty Clay Loam  | 0.089   |   0.43     | 1.0          | 1.23 | 1.944*10-7   |
+-----------------+---------+------------+--------------+------+--------------+
|Sandy Clay       | 0.100   |   0.38     | 2.7 	        | 1.23 | 3.333*10-7   |
+-----------------+---------+------------+--------------+------+--------------+
|Silty Clay       | 0.070   |   0.36     | 0.5          | 1.09 | 5.555*10-8   |
+-----------------+---------+------------+--------------+------+--------------+
|Clay             | 0.068   |   0.38     | 0.8          | 1.09 | 5.555*10-7   |
+-----------------+---------+------------+--------------+------+--------------+





.. [1] Carsel, R.F. & R.S.Parrish (1988): Developing joint probability distributions of soil water retention characteristics.-Water Resource Research 24:755-769.


