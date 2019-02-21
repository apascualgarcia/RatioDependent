Beyond-MeanField.f
==================

# DESCRIPTION:
  * This script builds or reads from files a bipartite network and integrates the system of equations. The scripts allows you to build the system either providing metaparameter values (features of the distributions around which random numbers will be generated) or reading the parameters and species abundances from files. You can read some of the parameters/abundances and generate others. You can additionally perturb the growth rates, controlled in the .in file with the Delta parameter. This should be typically a value between 0 and 1 representing the amount of perturbation relative to the growth rates (so 1 represents a perturbation of 100% of the growth rates, namely that the growth rates may randomly become zero). The script also allows you to build a system at a fix point, to perturb its growth rates next, which is the basis of our approach to structural stability (see Pascual-García and Bastolla, Nat Comm 2017). The idea is that you can fix all the parameters and solve the system of equations at a steady state for the growth rates free, and hence you build a feasible system. 


# INPUT
 
* If you provide the input files for the parameters or for the species abundances, the files should be named and formatted as follows. You can change these names formatting the file 1Beyond-MeanField_2-3.f in routine  Read_Parameters and recompiling.
 * _plantsIn.dat_ Plants abundances (a column vector)         
 * _animalsIn.dat_ Animal abundances
 * _betaInP.dat_ Competition for plants, a matrix (space separated with no header)
 * _betaInA.dat_ Competition for animals
 * _gammaInP.dat_ Mutualistic benefit of plants due to animals, a matrix (space separated with no header)      
 * _gammaInA.dat_ Mutualistic benefit of animals due to plants
 * _alphaInP.dat_ Growth rates of plants (a column vector)
 * _alphaInA.dat_ Growth rates of animals
 * _hhandlingInP.dat_ Handling time for plants due to animals abundances, a vector parameter hP (see documentation to clarify this)
 * _hhandlingInA.dat_ Handling time for animals due to plants abundances, hA
 * _ghandlingInP.dat_ Handling time for animals due to plants abundances, gP
 * _ghandlingInA.dat_ Handling time for plants due to plants abundances, gA
    
* You will also need a file called _Beyond-MeanField.in_ describing if you want to generate or read from file the parameters. To control for this, all the parameters have (at least) three metaparameters. The first metaparameter is called readX, with X the label for the three metaparameters. It will control if the parameter is read from the correspondent file (if it is equal to 1), or if it should be generated internally (=0). Then there are two additional parameters midX and widthX. The role of these parameters is:
 * If readX = 1, then, each value X0 readed from file will become: X=X0*midX*(1+widthX*rnd), where rnd = [-0.5,0.5]. You can consider different scenarios:
  * Keep just the readed value: Make midX=1 and widthX=0.
  * Randomize the readed value: Make midX=1 and widthX > 0.
  * Rescale the readed value: Make midX neq 0 and widthX eq 0.
  * Rescale and randomize the readed value: Make midX neq 0 and widthX > 0. This scenario and the previous one are useful for instance if you have only the topology of the interactions but you should explore the strengths.
 * If readX = 0, then the script will generate the following values: X=midX*(1+widthX*rnd), where rnd = [-0.5,0.5].

* In addition, some metaparameters require some more parameters to be generated. In particular, the competition matrix requires two more values to differentiate between inter and intraspecific competition, and the growth rates have one more parameter to control if it should be estimated at a fixed point. Specifically, the file _Beyond-MeanField.in_ must have the following lines in this order (see folder example1 for an example):
 * readNp         ! Control parameter for plants abundances, it is equal to 1 if you read Np from file, 0 if it will be generated internally, in which case:
 * midNp          ! It will generate the parameter with a distribution of plant densities Np=midNp*(1+ widthNp*rnd) where rnd = [-0.5,0.5]
 * widthNp        ! Note that these two parameters must be given even if readNp=1
 * readNa         ! Same for animal densities 
 * midNa          ! 
 * widthNa        ! 
 * readAlphaP     ! Plants growth rates (read or generate?)
 * midAlphaP      ! 
 * widthAlphaP    !
 * readAlphaA     ! Animals growth rates
 * midAlphaA      ! 
 * widthAlphaA
 * readBetaP      ! Sp x Sp competition between plants, equal to 1 if the matrix is readed from file, otherwise:
 * midBetaP       ! First parameter to generate the intraspecific competition for plants (diagonal terms)
 * widthBetaP     ! Width of the distribution, again: BetaP=midBetaP*(1+ widthBetaP*rnd)
 * midRhoP        ! Second parameter to generate the interespecific competition for animals (off-diagonal terms)
 * widthRhoP      ! Width of the distribution
 * readBetaA      ! Same for for competition between animals
 * midBetaA       ! Intraspecific competition for animals
 * widthBetaA     ! Width of the distribution
 * midRhoA        ! interspecific competition (animals)
 * widthRhoA      ! Width of the distribution
 * readGammaP     ! Sp x Sa interaction between the pool of plants and the pool of animals, equal to 1 if the matrix is readed from file, otherwise:
 * midGammaP      ! Generate a matrix with cells given by GammaP=midGammaP*(1+ widthGammaP*rnd)
 * widthGammaP    !
 * readGammaA     ! Same for the interaction between animals and plans
 * midGammaA      ! 
 * widthGammaA    !
 * readHp         ! Sp x 1, vector of handling time for plants. Set to 1 if you read from file, 0 if must be generated internally with
 * midHp          ! Hp=midHp*(1+ widthHp*rnd) 
 * widthHp        ! 
 * readHa         ! Sa x 1 vector of handling time for animals.
 * midHa          ! 
 * widthHa        !
 * readGp         ! Sp x 1, second kind of vector of handling time for plants. Set to 1 if you read from file, 0 if must be generated internally with
 * midGp          ! 
 * widthGp        ! 
 * readGa         ! Sa x 1, vector of handling time for animals
 * midGa          ! 
 * widthGa        ! 
 * f0P            ! Scale of the saturating term (real number, plants)
 * f0A            ! Scale of the saturating term (real number, animals)
 * Sp             ! Number of species of plants
 * Sa             ! of animals
 * Delta          ! Amplitude of the growth rates fluctuations at steady state: AlphaA=(Comp-Int)*(1+2*Delta*rnd), note the 2 factor ensuring that for Delta=1 Alpha may reach zero (rnd= [-0.5,0.5])
 * fixedPoint     ! Should the growth rates (Alpha) estimated from the rest of the model at the fixed point? (=1) or not (=0). This option will overide above parameters for Alpha. 
 * pathSummary    ! A name file for the output summary
 * seed           ! An input seed for random numbers

# OUTPUT

* The script will generate the following files. Note that, except for those related with the abundances that should always change any input file provided will bring as output file the same input file unless you rescale or randomize its values with the parameters midX and widthX. 
 * _plantsOut.dat_	# Plant abundances before the simulation starts, unless you simply read from file will differ from plantsIn.dat, if you randomize the values or generate them from scratch
 * _animalsOut.dat_ 	# Same for animals
 * _betaOutP.dat_ 	# The following files are in correspondence with the inputs, unless you generate
 * _betaOutA.dat_ 	# new parameters, you will get the same matrix
 * _gammaOutP.dat_
 * _gammaOutA.dat_
 * _alphaOutP.dat_  # For the specific case in which are estimated at fixed point, two additional columns are provided showing the total "competitive" and "interaction" terms
 * _alphaOutA.dat_ # from which the alpha are computed, i.e. alpha(i) = comp(i) - int(i)q
 * _hhandlingOutP.dat_ # Four new output  files, I leave the output files for paths the last ones
 * _hhandlingOutA.dat_
 * _ghandlingOutP.dat_
 * _ghandlingOutA.dat_
 * _plantsOut-Paths.dat_ 	# The trajectories of the abundances of plant species along the integration, one point is printed every 50 integration steps.
 * _animalsOut-Paths.dat_ 	# For the animals
 * _plantsOut-Final.dat_ # The input an output abundances, together with the alpha values and if the species got extincted. 
 * _animalsOut-Final.dat_ # 
 * Summary_$label.dat 	# A file providing the following information:
  * The parameters you used to run the simulation (one per row), and the following extra information
  * BiomassP  	#  total biomass for plants  
  * BiomassA   	#  total biomass of animals    
  * ExtinctP   	# plants extincted     
  * ExtinctA   	# animals     
  * SurvivingP   # plants species surviving     
  * SurvivingA   # animals     
  * NestednessPin 	#  Nestedness of plants matrix (defined as in Bastolla, Nature 2009, i.e. ecological overlap)
  * NestednessAin	#   Nestedness of animals     
  * NestednessTin	#   Total nestedness     
  * NestednessPout	#   Nestedness of plants species after the simulation, it changes just if there were extinctions     
  * NestednessAout	#   animals     
  * NestednessTout	#   total     
  * ConnectanceIn 	#   Starting connectance 
  * ConnectanceOut 	#   Output connectance
  * SinglesIn   	# Starting number of species with a single connection    
  * SinglesOut 	#   Final number    
  * SingleExtinct 	#   Species with a single connection that got extincted    
  * excludedPin 	#   Number of plant species forming fully independent modules (not exhaustively tested)
  * excludedPout 	#  Final
  * excludedAin 	#  Same for animals
  * excludedAout           
  * NegAlphaP 	#           Number of negative growth rates for plants (good control when growth rates are computed at the fixed point)
  * NegAlphaA 	#           for animals
  * AvAlphaP 	#   Average of the growth rates of plant species     
  * VarAlphaP 	#  Variance    
  * AvAlphaA  	#  Average of the growth rates of animal species          
  * VarAlphaA 	#  Variance   
  * AvDegreeInP 	#   Average of the degree of plant species     
  * VarDegreeInP 	#   Variance    
  * AvDegreeInA  	#    Average of the degree of animal species 
  * VarDegreeInA 	#   Variance     
  * Seed	#     Starting seed (control it changes along realizations)     
  * TrueDiversityP 	#    Effective biodiversity of plants     
  * TrueDiversityA 	#    Of animals
     
# COMPILATION:
  Simply run in the directory where the .f files are located. You will need a gfortran compiler, but no further requirements are
expected since all the libraries are included in the code file.
  ```$ make    ```
  Compilation has been tested in x86-64 procs under Linux.
# EXECUTION:
    Once the files are compiled, just execute the file, ensuring that the necessary files are in the same folder. If the output files already exist the script will not overwrite them and it aborts the execution.

# EXAMPLE:

At the folder _example1_ there is the input and output of a mutualistic system with all parameters generated from scracth, leading
to mean-field matrices, and growth rates estimated at a fixed point.

# TREE 
From root (Main):

*     ..... Read_Parameters .... Header                 ... When  .. idate                                   ##(At 1Beyond-MeanField_2-3.f)
*                                                                 .. itime
*                           .... Warning
*
*     ..... Pseudo_Main     .... Read_Generate_Inputs                                                        ##(At 2Parameters.f)
*                                                       ... ReadGenerateVec       
*                                                       ... ReadGenerateTrixComp
*                                                       ... ReadGenerateTrixInt
*                                                       ... Topology       
*                                                       ... Consistent_Alpha .. ran2 
*
*                           .... Integration            ... Derivatives                                      ##(At 3Integration.f)     
*                                                       ... bsstep           .. mmid      . Derivatives
*                                                                            .. pzextr
*                           .... Output_Globals         ... Topology                                         ##(back to 1Beyond-MeanField_2-3.f)
*
*     ..... Normal_Exit
*

Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl
=================================================

# DESCRIPTION

This script is a wrapper to call multiple times Beyond-MeanField.f with a fork manager, to explore
different regimes of parameters. It is designed to explore the structural stability of certain system, and
hence it allows for the exploration of a range of Delta parameters (i.e. perturbations of the growth rates), and
it will run for each Delta value N randomizations, and it will return the mean and standard deviation of all
the quantities included in the Summary file of Beyond-MeanField.f. In addition, it allows to explore
different values of one parameter X of your choice among Beta, Rho, Gamma, H and G parameters. Therefore, the
script will run N times a combination of X, Delta parameters (keeping all the remaining parameters fixed, except
for any eventual randomization you may want to include as in the Fortran code), and
will return a list with means and stdv for each combination. 

# REQUIREMENTS

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

# USAGE

To use the script you should first edit in the code a list of parameters. These parameters are:

* The parameters needed to run Beyond-MeanField.f, i.e. the parameters contained in Beyond-MeanField.in. You should fix there the different choices, if the script should read or generate parameters, etc. Since the idea behind this script is to explore the structural stability, the parameter fixedPoint should be always equal to one, in which case the system will be always prepared at a feasible starting point.
* The range of Delta values to simulate.
* Which is the additional variable X you want to explore and the range of values.
* The number of randomizations for each pair X, Delta values.
* The number of processors to use (variable Nproc). If it is not fixed to certain value it will get the maximum number available in the computer.
* Labels for the output files.

Then, you should run the script in the directory where the binary file BeyondMeanField_gfortran and your .in files (if needed) are located:

```$ Shuttle4_Beyond-MeanField-Adaptive_ForkManager.pl```

To explore a particular region of the Delta parameter (for instance the one in which the probability of observing at least one extinction is 0.5, our original definition) you can explore a wide range of Delta with few randomizations and then narrow your search around the region where you observe few extinctions, further increasing the number of randomizations.

# OUTPUT

The script will return three files:

* _Summary-Ensemble_Beyond-MeanField_$labels_ This file returns, for each pair of metaparameters X and perturbation Delta, the mean and standard deviation across the N randomizations of all the output variables of the FORTRAN script.
* _ProbExtinctDegrees_Beyond-MeanField_$labels_ This file quantifies, for each pair of metaparameters X and perturbation Delta, the probability that an extinction occurs for a species with starting degree d, where the degree is the number of interactions in the interaction matrix Gamma. Each row contains, the two metaparameters, the degree d (from 1 to the number of species) and the probability P(extinction | d).
* _DevfromRefSimulation_Beyond-MeanField_$labels_ Similar file than the _Summary_ file but in this case it computes the difference between the means of the simulations for the first pair of parameters X and Delta, and the remaining ones. The idea behind this file is to test the increase (decrease) in the different values with respect to some reference simulation, typically the ones in which the metaparameter X is equal to zero, in order to test the effect of including the metaparameter in the model. For instance, if Gamma is a mutualistic matrix, starting with midGamma=0 we will simulate two purely competitive systems, and we can test the difference between having pure competition and mutialism.

In folder Perl/example1 you will find the output of a mean field mutualistic simulation with four values of Delta and three values of midGamma

# PERFORMANCE

For the moment this wrappper simply runs multiple times the Fortran code in independent directories and then retrieves and post-process all the information, what implies reading and writing files and increases the time of computation. In the (I hope near) future I will wrap it into R to avoid read/write operations, but I'm not familiar with fork management in R yet. At this moment, the example presented at Perl/example1 consists of 6000 population dynamics simulations and it takes 5 minutes running in a 12 core AMD computer with two threads per core, and it required around 4GB of memory.





