Beyond-MeanField.f
==================

# DESCRIPTION:
  * This script builds or reads from files a bipartite network and integrates the system of equations. The scripts allows you to build the system either providing metaparameter values (features of the distributions around which random numbers will be generated) or reading the parameters and species abundances from files. You can read some of the parameters/abundances and generate others. 
 * The script will integrate the system you provide with files or that you build internally. You can additionally perturb the growth rates, controlled in the .in file with the Delta parameter. This should be typically a value between 0 and 1 representing the amount of perturbation relative to the growth rates (so 1 represents a perturbation of 100% of the growth rates, and you will get extinctions). The script also allows you to build a system at a fix point, to perturb its growth rates next, which is the basis of our approach to structural stability (see Pascual-García and Bastolla, Nat Comm 2017). The idea is that you can fix all the parameters and solve the system of equations at the fix point leaving the growth rates free, that will be fixed with this procedure. 


# INPUT
 
* If you provide the input files for the parameters or for the species abundances, the files should be named and formatted as follows. You can change this names formatting the file 1Beyond-MeanField_2-3.f and recompiling.
 * _plantsIn.dat_ Plants abundances (a column vector)         
 * _animalsIn.dat_ Animal abundances
 * _betaInP.dat_ Competition for plants, a matrix (space separated with no header)
 * _betaInA.dat_ Competition for animals
 * _gammaInP.dat_ Mutualistic benefit of plants due to animals, a matrix (space separated with no header)      
 * _gammaInA.dat_ Mutualistic benefit of animals due to plants
 * _alphaInP.dat_ Growth rates of plants (a column vector)
 * _alphaInA.dat_ Growth rates of animals
 * _hhandlingInP.dat_ Handling time for plants due to animals abundances, parameter hP (see documentation to clarify this)
 * _hhandlingInA.dat_ Handling time for animals due to plants abundances, hA
 * _ghandlingInP.dat_ Handling time for animals due to plants abundances, gP
 * _ghandlingInA.dat_ Handling time for plants due to plants abundances, gA
    
* You will also need a file called _Beyond-MeanField.in_ describing if you want to generate or read from file the parameters. To control for this, all the parameters have (at least) three metaparameters. The first metaparameter is called readX, with X the label for the three metaparameters. It will control if the parameter is read from the correspondent file above (if it is equal to 1), or if it should be generated internally (=0). Then there are two additional parameters midX and widthX. The role of these parameters is:
 * If readX = 1, then, each value X0 readed from file will be randomized if widthX is not zero, as: X=X0*(1+widthX*rnd), where rnd = [-0.5,0.5]. (i.e. the value readed will take a new value of a distribution centered at the value readed with some dispersion controlled by widthX).
 * If readX = 0, then the script will generate the following values: X=midX*(1+widthX*rnd), where rnd = [-0.5,0.5].

* In addition, some metaparameters require some more parameters to be generated. In particular, the competition matrix requires two more values to differentiate between inter and intraspecific competition, and the growth rates have one more parameter to control if should be estimated at a fixed point. Specifically, the file _Beyond-MeanField.in_ must have the following lines in this order:
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
 * Delta          ! Amplitude of the growth rates fluctuations
 * fixedPoint     ! Should the growth rates (Alpha) estimated from the rest of the model at the fixed point? (=1) or not (=0). This option will overide above parameters for Alpha. 

# Output
* The script will generate the following files. Note that, except those related with the abundances that should always change (unless you build a fix point and you do not perturb it, in which case there will not be integration), any input provided (not generated internally) will bring as output file simply the input file. 
 * _plantsOut-Final.out_	#  Plant abundances at the steady state    
 * _animalsOut-Final.out_ 	# Same for animals
 * _betaOutP.out_ 	# The following files are in correspondence with the inputs, unless you generate
 * _betaOutA.out_ 	# new parameters, you will get the same matrix
 * _gammaOutP.out_
 * _gammaOutA.out_
 * _alphaOutP.out_
 * _alphaOutA.out_
 * _hhandlingOutP.dat_ ! Four new output  files, I leave the output files for paths the last ones
 * _hhandlingOutA.dat_
 * _ghandlingOutP.dat_
 * _ghandlingOutA.dat_
 * _plantsOut-Paths.out_ 	# The trajectories of the abundances of plant species along the integration 
 * _animalsOut-Paths.out_ 	# For the animals
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
  * excludedPin 	#   Starting plant species forming an independent module
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
  * EntropyP 	#    Effective biodiversity of plants     
  * EntropyA 	#    Of animals
     
# COMPILATION:
  Simply run in the directory where the .f files are located. You will need a gfortran compiler.
  $ make    
  Compilation has been tested in x86-64 procs under Linux.
# EXECUTION:
    Once the files are compiled, just execute the file, ensuring that the necessary files are in the same folder. If the output files already exist the script will not overwrite them and it aborts the execution.

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


