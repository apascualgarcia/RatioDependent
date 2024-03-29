RatioDependent Project
======================

This project consists of two main scripts. A Fortran code to integrate population dynamics, the proyect ```Beyond-MeanField```, and a Perl
wrapper of this script to run it massively using multiple cores ```Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl```. Both
scripts are mainly designed to compute the structural stability of the system. The former allows for the construction
of the system at a feasible state and its perturbation of the growth rates by a certain amplitude. And the latter allows
for the exploration of multiple amplitude values of the perturbations through multiple realizations of the noise. Both scripts are described in the following.

Beyond-MeanField (FORTRAN)
==========================

## DESCRIPTION:

This script builds or reads from files a bipartite network and integrates the system of equations. It is coded in three different files: ```1Beyond-MeanField_2-3.f```, ```2Parameters.f``` and ```3Integration_Adaptive.f```. You will find some description of the contents along these lines and a tree with the different subroutines describing where are located below. It is provided an executable file ```BeyondMeanField_gfortran```, but compilation is system-dependent, so it is advisable to recompile it even if you share the features of the compilation described below.

 The scripts allows you to build the system either providing metaparameter values (features of the distributions around which random numbers will be generated) or reading the parameters and species abundances from files. You can read some of the parameters/abundances and generate others. 
In addition, the script also allows you to build a system at a fixed point and to perturb its growth rates next, which is the basis of our approach to structural stability (see Pascual-García and Bastolla, _Nat Comm_ 2017). The idea is that you can fix all the parameters and solve the system of equations at a steady state for the growth rates (to obtain a feasible system) and then perturb the growth rates and integrate the ODEs. The perturbation of the growth rates is controlled by a specific (Delta) parameter, which should be typically a value between 0 and 1 representing the amount of perturbation relative to the growth rates (so 1 represents a perturbation of 100% of the growth rates, namely that the growth rates may randomly become zero). 

In the following, we label the two pools as plants (P) and animals (A). These can be two pools interacting mutualistically (pollinators and its flowering plants), in a prey-predator way (like plants and herbivores) or even competitively, depending on the sign of the interaction matrices between both pools chosen (see below). All the parameters have one parameter for plants and another for animals, being in this way possible to control them independently.

Each pool has a parameter for the growth rates (Alpha), which is positive or negative and that, as we said, can be inferred at a steady state. In addition, it has
a competition matrix with respect to the members of the own pool. The script can generate internally a mean-field competition matrix which is controlled by two parameters Beta (diagonal terms, intraspecific competition) and Rho (off-diagonal terms, interspecific competition). These coefficients must be positive. Finally, there is an interaction term between both pools including a ratio-dependent functional response. We define this term for plants as:

```
F^P_ij = Gamma^P_ij N^P_i N^A_j/[f^P_0 + g^A_j sum_k |Gamma^P_kj| N^P_k + h^P_i sum_l |Gamma^P_il| N^A_l]
```

and symmetrically for animals, just interchanging the superscripts ```A``` and ```P```. If the Gamma matrices are generated internally, they will be a mean field matrix whose strength is given by the user, and the sign will determine the kind of interaction between the pools. In addition, we have two handling times ```h``` and ```g```. Note the notation followed for ```g``` and ```h```, the superindex ```P``` (plants) or ```A``` (animals) is in correspondence with the dimension of space (plants or animals) in which it is defined, and hence this notation is not necessarily the most meaningful ecologically.

## INPUT
 
1. If you provide the input files for the parameters or for the species abundances, the files should be named and formatted as follows. You can change these names formatting the file ```1Beyond-MeanField_2-3.f``` in routine  ```Read_Parameters``` and recompiling.
     * __plantsIn.dat__ --> Plants abundances (a column vector)         
     * __animalsIn.dat__ --> Animal abundances
     * __betaInP.dat__ --> Competition for plants, a matrix (space separated with no header)
     * __betaInA.dat__ --> Competition for animals
     * __gammaInP.dat__ --> Matrix of interaction of plants (rows) and animals (columns), (space separated with no header)      
     * __gammaInA.dat__ --> Matrix of interaction of animals (rows) plants (columns)
     * __alphaInP.dat__ --> Growth rates of plants (a column vector)
     * __alphaInA.dat__ --> Growth rates of animals
     * __hhandlingInP.dat__ --> Saturation term (vector h^P_i of Sp plants): ```h^P_i sum_l |Gamma^P_il| N^A_l]```
     * __hhandlingInA.dat__ --> Saturation term (vector h^A_i of Sa plants): ```h^A_i sum_l |Gamma^A_il| N^P_l]```
     * __ghandlingInP.dat__ --> Saturation term (vector g^P_i of Sp plants): ```g^P_j sum_k |Gamma^A_kj| N^A_k]```
     * __ghandlingInA.dat__ --> Saturation term (vector g^A_i of Sa plants): ```g^A_j sum_k |Gamma^P_kj| N^P_k]```    
2.  You will also need a file called ```Beyond-MeanField.in``` describing if you want to generate or read from file the parameters. To control for this, all the parameters have (at least) three metaparameters. The first metaparameter is called ```readX```, with ```X``` the label for the three metaparameters. It will control if the parameter is read from the correspondent file (if it is equal to 1), or if it should be generated internally (=0). Then there are two additional parameters ```midX``` and ```widthX```. The role of these parameters is:
     * If ```readX = 1```, then, each value ```X0``` readed from file will become: ```X=X0*midX*(1+widthX*rnd)```, where ```rnd = [-0.5,0.5]```. You can consider different scenarios:
        * Keep just the readed value: Make ```midX=1``` and ```widthX=0```.
        * Randomize the readed value: Make ```midX=1``` and ```widthX > 0```.
        * Rescale the readed value: Make ```midX neq 0``` and ```widthX = 0```.
        * Rescale and randomize the readed value: Make ```midX neq 0``` and ```widthX > 0```. This scenario and the previous one are useful for instance if you have only the topology of the interactions but you should explore the strengths.
     * If ```readX = 0```, then the script will generate the following values: ```X=midX*(1+widthX*rnd)```, where ```rnd = [-0.5,0.5]```.

* In addition, some objects require some more parameters to be generated. In particular, the competition matrix requires two more values to differentiate between inter and intraspecific competition, and the growth rates have one more parameter to control if it should be estimated at a fixed point. Specifically, the file _Beyond-MeanField.in_ must have the following lines in this order (see folder example1 for an example):
     * __readNp__         ! Control parameter for plants abundances, it is equal to 1 if you read Np from file, 0 if it will be generated internally, in which case:
     * __midNp__          ! It will generate the parameter with a distribution of plant densities ```Np=midNp*(1+ widthNp*rnd)``` where ```rnd = [-0.5,0.5]```
     * __widthNp__        ! Note that these two parameters must be given even if ```readNp=1```
     * __readNa__         ! Same for animal densities 
     * __midNa__          ! 
     * __widthNa__        ! 
     * __readAlphaP__     ! Plants growth rates (read or generate?)
     * __midAlphaP__      ! 
     * __widthAlphaP__    !
     * __readAlphaA__     ! Animals growth rates
     * __midAlphaA__      ! 
     * __widthAlphaA__
     * __readBetaP__      ! Sp x Sp competition between plants, equal to 1 if the matrix is readed from file, otherwise:
     * __midBetaP__       ! First parameter to generate the intraspecific competition for plants (diagonal terms)
     * __widthBetaP__     ! Width of the distribution, again: ```BetaP=midBetaP*(1+ widthBetaP*rnd)```
     * __midRhoP__        ! Second parameter to generate the interespecific competition for animals (off-diagonal terms)
     * __widthRhoP__      ! Width of the distribution
     * __readBetaA__      ! Same for for competition between animals
     * __midBetaA__       ! Intraspecific competition for animals
     * __widthBetaA__     ! Width of the distribution
     * __midRhoA__        ! interspecific competition (animals)
     * __widthRhoA__      ! Width of the distribution
     * __readGammaP__     ! Sp x Sa interaction between the pool of plants and the pool of animals, equal to 1 if the matrix is readed from file, otherwise:
     * __midGammaP__      ! Generate a matrix with cells given by ```GammaP=midGammaP*(1+ widthGammaP*rnd)```
     * __widthGammaP__    !
     * __readGammaA__     ! Same for the interaction between animals and plans
     * __midGammaA__      ! 
     * __widthGammaA__    !
     * __readHp__         ! Sp x 1, vector of handling time for plants. Set to 1 if you read from file, 0 if must be generated internally with
     * __midHp__          ! ```Hp=midHp*(1+ widthHp*rnd) ```
     * __widthHp__        ! 
     * __readHa__         ! Sa x 1 vector of handling time for animals.
     * __midHa__          ! 
     * __widthHa__        !
     * __readGp__         ! Sp x 1, second kind of vector of handling time for plants. Set to 1 if you read from file, 0 if must be generated internally with
     * __midGp__          ! 
     * __widthGp__        ! 
     * __readGa__         ! Sa x 1, vector of handling time for animals
     * __midGa__          ! 
     * __widthGa__        ! 
     * __f0P__            ! Scale of the saturating term (real number, plants)
     * __f0A__            ! Scale of the saturating term (real number, animals)
     * __Sp__             ! Number of species of plants
     * __Sa__             ! of animals
     * __Delta__          ! Amplitude of the growth rates fluctuations at steady state: ```AlphaA=(Comp-Int)*(1+2*Delta*rnd)```, note the 2 factor ensuring that for ```Delta=1```, Alpha may become equal to zero since (rnd= [-0.5,0.5])
     * __fixedPoint__     ! Should the growth rates (Alpha) estimated from the rest of the model at the fixed point? (=1) or not (=0). This option will overide above parameters for Alpha. 
     * __pathSummary__    ! A name file for the output summary
     * __seed__           ! An input seed for random numbers

# OUTPUT

* The script will generate the following files. Note that, except for those related with the abundances that should always change any input file provided will bring as output file the same input file unless you rescale or randomize its values with the parameters midX and widthX. 
     * __plantsOut.dat__	# Plant abundances before the simulation starts, unless you simply read from file will differ from plantsIn.dat, if you randomize the values or generate them from scratch
     * __animalsOut.dat__ 	# Same for animals
     * __betaOutP.dat__ 	# The following files are in correspondence with the inputs, unless you generate
     * __betaOutA.dat__ 	# new parameters, you will get the same matrix
     * __gammaOutP.dat__ 
     * __gammaOutA.dat__ 
     * __alphaOutP.dat__  # For the specific case in which are estimated at fixed point, two additional columns are provided showing the total "competitive" and "interaction" terms
     * __alphaOutA.dat__ # from which the alpha are computed, i.e. ```alpha(i) = comp(i) - int(i)```
     * __hhandlingOutP.dat__ # Four new output  files, I leave the output files for paths the last ones
     * __hhandlingOutA.dat__ 
     * __ghandlingOutP.dat__
     * __ghandlingOutA.dat__
     * __plantsOut-Paths.dat__ 	# The trajectories of the abundances of plant species along the integration, one point is printed every 50 integration steps.
     * __animalsOut-Paths.dat__ 	# For the animals
     * __plantsOut-Final.dat__ # The input an output abundances, together with the alpha values and if the species got extincted. 
     * __animalsOut-Final.dat__ # 
     * __Summary\_$label.dat__ 	# A file providing the following information:
         * The parameters you used to run the simulation (one per row), and the following extra information:
         * __BiomassP__  	#  total biomass for plants  
         * __BiomassA__   	#  total biomass of animals    
         * __ExtinctP__   	# plants extincted     
         * __ExtinctA__   	# animals     
         * __SurvivingP__   # plants species surviving     
         * __SurvivingA__   # animals     
         * __NestednessPin__ 	#  Nestedness of plants matrix (defined as in Bastolla, Nature 2009, i.e. ecological overlap)
         * __NestednessAin__	#   Nestedness of animals     
         * __NestednessTin__	#   Total nestedness     
         * __NestednessPout__	#   Nestedness of plants species after the simulation, it changes just if there were extinctions     
         * __NestednessAout__	#   animals     
         * __NestednessTout__	#   total     
         * __ConnectanceIn__ 	#   Starting connectance 
         * __ConnectanceOut__ 	#   Output connectance
         * __SinglesIn__   	# Starting number of species with a single connection    
         * __SinglesOut__ 	#   Final number    
         * __SingleExtinct__ 	#   Species with a single connection that got extincted    
         * __excludedPin__ 	#   Number of plant species forming fully independent modules (not exhaustively tested)
         * __excludedPout__ 	#  Final
         * __excludedAin__ 	#  Same for animals
         * __excludedAout__           
         * __NegAlphaP__ 	#           Number of negative growth rates for plants (good control when growth rates are computed at the fixed point)
         * __NegAlphaA__ 	#           for animals
         * __AvAlphaP__ 	#   Average of the growth rates of plant species     
         * __VarAlphaP__ 	#  Variance    
         * __AvAlphaA__  	#  Average of the growth rates of animal species          
         * __VarAlphaA__ 	#  Variance   
         * __AvDegreeInP__ 	#   Average of the degree of plant species     
         * __VarDegreeInP__ 	#   Variance    
         * __AvDegreeInA__  	#    Average of the degree of animal species 
         * __VarDegreeInA__ 	#   Variance     
         * __Seed__	#     Starting seed (control it changes across realizations)     
         * __TrueDiversityP__ 	#    Effective biodiversity of plants     
         * __TrueDiversityA__ 	#    Of animals
     
## COMPILATION:
  
Simply run in the directory where the .f files are located. You will need a gfortran compiler, but no further requirements are
expected since all the libraries are included in the code file.

  ```$ make    ```

  Compilation has been tested in x86-64 procs under Linux.

## EXECUTION:

Once the files are compiled, just execute the file, ensuring that the necessary files are in the same folder. If the output files already exist the script will not overwrite them and it aborts the execution.

## EXAMPLE:

At the folder _example1_ there is the input and output of a mutualistic system with all parameters generated from scracth, leading
to mean-field matrices, and growth rates estimated at a fixed point.

## TREE of the code structure 
From root (Main) we show from top o down the sequential order in which the code is executed and how the subroutines are nested one within the other. Subroutines with less "dots" are closer to the leafs and are called by the immediately upper subroutine containing one dot more.

     ..... Read_Parameters .... Header                 ... When  .. idate                                   ##(At 1Beyond-MeanField_2-3.f)
                                                           
                                                                 .. itime
                        .... Warning

     ..... Pseudo_Main     .... Read_Generate_Inputs                                                        ##(At 2Parameters.f)
                                                       ... ReadGenerateVec       
                                                       ... ReadGenerateTrixComp
                                                       ... ReadGenerateTrixInt
                                                       ... Topology       
                                                       ... Consistent_Alpha .. ran2 

                           .... Integration            ... Derivatives                                      ##(At 3Integration.f)     
                                                       ... bsstep           .. mmid      . Derivatives
                                                                            .. pzextr
                           .... Output_Globals         ... Topology                                         ##(back to 1Beyond-MeanField_2-3.f)

     ..... Normal_Exit


Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl
=================================================

## DESCRIPTION

This script (located in the folder Perl) is a wrapper to call multiple times Beyond-MeanField.f with a fork manager, to explore
different regimes of parameters. It is designed to explore the structural stability of certain system, and
hence it allows for the exploration of a range of Delta parameters (i.e. perturbations of the growth rates), and
it will run for each Delta value N randomizations, and it will return the mean and standard deviation of all
the quantities included in the Summary file of Beyond-MeanField.f. In addition, it allows to explore
different values of one parameter X of your choice among Beta, Rho, Gamma, H and G parameters. Therefore, the
script will run N times a combination of (X, Delta) parameters (keeping all the remaining parameters fixed, except
for any eventual randomization you may want to include each time you run a realization), and
will return a list with means and stdv for each combination. 

## REQUIREMENTS

To run this script you will need to have installed Perl and, in particular, the package
ForkManager to run multiple times the scripts using different processors in the same computer.

* First install [CPAN:](https://www.cpan.org/modules/INSTALL.html)
* Then open a CPAN shell typing in your Linux shell:
     * ```$> sudo perl -MCPAN -e shell```
* A new prompt (e.g. $CPAN>) will appear indicating that you are at CPAN, then type:
     * ```$CPAN> install Parallel::ForkManager```
     * ```$CPAN> exit```
* Remember to give execution permission to the script:
     * ```$> chmod 755 Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl```

## USAGE

To use the script you should first edit in the code a list of parameters. These parameters are:

* The parameters needed to run Beyond-MeanField.f, i.e. the parameters contained in Beyond-MeanField.in. You should fix there the different choices, if the script should read or generate parameters, etc. Since the idea behind this script is to explore the structural stability, the parameter fixedPoint should be always equal to one, in which case the system will be always prepared at a feasible starting point.
* The range of Delta values to simulate.
* Select the additional variable X you want to explore and the range of values (select a single value for any of the possible choices if you are not interested in exploring any metaparameter)
* The number of randomizations for each pair (X, Delta) values.
* The number of processors to use (variable Nproc). If it is not fixed to certain value it will get the maximum number available in the computer.
* Labels for the output files.

Then, you should run the script in the directory where the binary file BeyondMeanField_gfortran and your .in files (if needed) are located:

```
$ chmod 755 Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl
$ ./Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl
```


## OUTPUT

The script will return three files:

* __Summary-Ensemble\_Beyond-MeanField\_$labels__ This file returns, for each pair of metaparameters X and perturbation Delta, the mean and standard deviation across the N randomizations of all the output variables of the FORTRAN script.
* __ProbExtinctDegrees\_Beyond-MeanField\_$labels__ This file quantifies, for each pair of metaparameters X and perturbation Delta, the probability that an extinction occurs for a species with starting degree d, where the degree is the number of interactions in the interaction matrix Gamma. Each row contains, the two metaparameters, the degree d (from 1 to the number of species) and the probability P(extinction | d).
* __DevfromRefSimulation\_Beyond-MeanField\_$labels__ Similar file than the _Summary_ file but in this case it computes the difference between the means of the simulations for the first pair of parameters X and Delta, and the remaining ones. The idea behind this file is to test the increase (decrease) in the different values with respect to some reference simulation, typically the ones in which the metaparameter X is equal to zero, in order to test the effect of including the metaparameter in the model. For instance, if Gamma is a mutualistic matrix, starting with midGamma=0 we will simulate two purely competitive systems, and we can test the difference between having pure competition and mutialism.

In folder Perl/example1 you will find the output of a mean field mutualistic simulation with four values of Delta and three values of midGamma.

_Example_ Critical perturbation:_ To compute the critical perturbation of a system for some metaparameters given, you should fix the values of the metaparameters
and explore a particular region of the Delta parameter. If you do not have a clue about the specific region to explore, you may want to start exploring a a wide range of Delta with few randomizations and then narrow your search around the region where you observe few extinctions, further increasing the number of randomizations.

The relevant file for this exploration is the file  __Summary-Ensemble\_Beyond-MeanField\_$labels__. In particular, column 1 returns the value of the metaparameter you are dynamically modifying, column 2 the value of Delta, column 7 (9) the mean number of extinctions for the pool P (A) of organisms, and column 8 (10) the standard error for the pool P (A). For instance, considering these columns in the example located at folder ```example1``` within the folder ```Perl```, these columns read:

```
#1midGamma 2Delta, 7ExtinctP, 8StdvExtP, 9ExtinctA, 10StdvExtA,
0.02 0.1 0 0 0 0
0.02 0.25 0 0 0 0
0.02 0.5 6.39 4.86 5.83 4.12
0.02 0.75 13.59 5.76 13.09 5.58
0.04 0.1 0 0 0 0
0.04 0.25 0 0 0 0
0.04 0.5 0.096 0.15 0.036 0.035
0.04 0.75 7.578 4.66 6.91 5.073
0.06 0.1 0 0 0 0
0.06 0.25 0 0 0 0
0.06 0.5 0 0 0 0
0.06 0.75 0.02 0.020 0 0
```

which suggests that, for the metaparameter Gamma=0.02 the Delta critical should be explored for Delta<0.5, for Gamma=0.04, between Delta=0.5 and Delta=0.75, and
for Gamma=0.06, the values should be Delta>0.75.


## PERFORMANCE

For the moment this wrappper simply runs multiple times the Fortran code in independent directories and then retrieves and post-process all the information, what implies reading and writing files and increases the time of computation. In the (I hope near) future I will wrap it into R to avoid read/write operations, but I'm not familiar with fork management in R yet. At this moment, the example presented at Perl/example1 consists of 6000 population dynamics simulations and it takes 5 minutes running in a 12 core AMD computer with two threads per core, and it required around 4GB of memory. 





