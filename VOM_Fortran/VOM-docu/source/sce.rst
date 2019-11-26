
Shuffled Complex Evolutionary algorithm
===========================================
The Shuffled Complex Evolutionary algorithm of Duan et al. (1994) is used to optimize the VOM. The main settings can be found in the VOM_namelist:

+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_ncomp          | MAXIMUM NUMBER OF COMPLEXES (p)                                                                    | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_ncompmin       | MINIMUM NUMBER OF COMPLEXES (pmin)                                                                 | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_resolution     | RESOLUTION OF OPTIMISATION (% OF MAX VARIATION WHEN OPTIMISATION STOPS)                            | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_patience       | NUMBER OF LOOPS WITHOUT INCREASE IN OF BEFORE OPTIMISATION STOPS                                   | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_nsimp          | NUMBER OF OPTIMISATIONS PER COMPLEX AND RUN                                                        | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_focus          | IF <1.0, THE SPREAD OF THE RANDOM SEED AROUND THE INITIAL VALUES IS LIMITED                        | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|i_iter           | Maximum iterations in case of random runs                                                          | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|vom_npar         | number of parameters in shuffle2par used for optimization in SCE                                   | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|n_thread         | number of threads to be used in parallel (one complex per thread)                                  | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|sce_restart      | restart SCE from previous run, TRUE or FALSE                                                       | \-       |
+-----------------+----------------------------------------------------------------------------------------------------+----------+
|runtime_limit    | time in minutes before sce stops at it's earliest possibility                                      | minutes  |
+-----------------+----------------------------------------------------------------------------------------------------+----------+

The SCE-algorithm (mode 1) can be run in parallel by setting the number of threads. When specified, the SCE-algorithm runs over the 
different complexes in parallel (one complex per thread). The maximum runtime (in minutes) can be set, after which the algorithm 
tries to stop at the earliest possibility. Afterwards, the SCE-algorithm can be restarted when sce_restart is set to TRUE and the 
files of the previous round are available.


Outputs
-------------------------------

+-----------------+---------------------------------------------------------------------------------------------------------------+
|Filename         | Description                                                                                                   |
+-----------------+---------------------------------------------------------------------------------------------------------------+
|sce_progress.txt |  Gets progressively filled with messages as the model runs.                                                   |
+-----------------+---------------------------------------------------------------------------------------------------------------+
|sce_out.txt      | Gets progressively filled with an experimental parameter set and the respective value of the                  |
|                 | objective function. Contains a line for every parameter set explored, composed of the parameter values        | 
|                 | followed by the value of the objective function.                                                              |
+-----------------+---------------------------------------------------------------------------------------------------------------+
|sce_lastloop.txt | Written at the end of each optimisation loop and contains all information needed to continue with the next    |
|                 | loop, i.e. the number of complexes the number of previous loops, the number of runs performed already, the    |
|                 | number of runs since the best objective function was achieved, followed by all the parameter sets explored    |                                       
|                 | in the last completed loop and their respective objective function values (similar as sce_out.txt, but limited|
|                 | to the last loop only).                                                                                       |
+-----------------+---------------------------------------------------------------------------------------------------------------+
|sce_lastbest.txt | Contains the best parameter set and objective function of the most recent loop.                               |
+-----------------+---------------------------------------------------------------------------------------------------------------+
|sce_bestpars.txt | Written whenever a parameter set is found that yields the best objective function value.                      |
+-----------------+---------------------------------------------------------------------------------------------------------------+
|sce_status.txt   | Created when optimisation finished successfully, contains one ascii symbol: "1"                               |
+-----------------+---------------------------------------------------------------------------------------------------------------+

