Beyond-MeanField.f
==================

* DESCRIPTION:
  * This script builds a bipartite mutualistic network and integrates the system of equations. The scripts allows you to build the system either providing metaparameter values (features of the distributions around which numbers  will be generated) or reading the parameters and species abundances from files. You can read some of the parameters/abundances and generate others, the only mandatory input is a network of mutualistic interactions between both pools, typially with plants in rows and animals in columns. If you have other pools keep in mind that, in the following, parameters are distinguished as plants and animals, and one should be consistent.
  * If you provide the input files for the parameters or for the species abundances, the files should be named and formatted follows. You can change this names formatting the file 1Beyond-MeanField_2-3.f.
    * 'plantsIn.dat' Plants abundances (a column vector)         
    * 'animalsIn.dat' Animal abundances
    * 'betaInP.dat' Competition for plants, a matrix (space separated with no header)
    * 'betaInA.dat' Competition for animals
    * 'gammaInP.dat' Mutualistic benefit of plants due to animals, a matrix (space separated with no header)      
    * 'gammaInA.dat' Mutualistic benefit of animals due to plants
    * 'alphaInP.dat' Growth rates of plants (a column vector)
    * 'alphaInA.dat' Growth rates of animals
    
  * You will also need a file called 'Beyond-MeanField.in' with the following lines (We provide some values as examples and after '#' we comment the description):

     * 100	# Plants, Population densities, a uniform distribution centered around this value is generated if >0, if =0 read from file.
     *  5	# Plants, Width of the distribution (+- "value" around the center of the distro. Set it to 0 to get constant values)
     * 	1	# Animals, Population densities, center
     * 	0.05	# Animal, Width
     * 	-1	# Alphas (Growth rates, >0 generate random numbers, =0, read from file, <0 solve at the fixed point the equations after fixing from the other parameters, see below)
     * 	0.05	# Alphas, Width
     * 	1	# Beta0 parameter (intraspecific competition), (>0 random numbers, =0 read from file), it multiplies a distribution centered in one (i=j) or rho (i!=j) with the following width:
     * 	0.15	# Beta0, Width
     * 	0.15	# Rho parameter, interspecific competition (always constant, no noise here)
     * 	1	# Gamma0  (>0 random numbers, <0 read from file), it multiplies a distribution centered in one with the following width:
     * 	0.15	# Gamma0, Width
     * 	1	# Handling time (Plants)
     * 	1	# Handling time (Animals)
     * 	50	# Species number (Plants)
     * 	14	# Species number (Animals)
     * 	0.5	# Delta parameter, amplitude of the perturbation over the growth rates
     * 	../../Trinidad.txt	# Path for the adjacency matrix
     * 	Summary_Trinidad-Observed_Gamma0-1_Delta-0.5_2strong.dat	# Path for the output summary matrix
     * 	0.6890  # This parameter is only available in the "multiple" version. It is needed when you run multiple times the script, and it should change through realizations
	
  * The script will generate the following files. Note that, except those related with the abundances that should always change (unless you build a fix point and you do not perturb it, in which case there will not be integration), any input provided (not generated internally) will bring as output file simply the input file. 
    * 'plantsOut-Final.out'	#  Plant abundances at the steady state    
    * 'animalsOut-Final.out' 	# Same for animals
    * 'betaOutP.out' 	# The following files are in correspondence with the inputs, unless you generate
    * 'betaOutA.out' 	# new parameters, you will get the same matrix
    * 'gammaOutP.out'
    * 'gammaOutA.out'
    * 'alphaOutP.out'
    * 'alphaOutA.out'
    * 'plantsOut-Paths.out' 	# The trajectories of the abundances of plant species along the integration 
    * 'animalsOut-Paths.out' 	# For the animals
    * Summary_$label.dat 	# A file providing the following information:
     > The parameters you used to run the simulation (one per row), and the following extra information
     > BiomassP  	#  total biomass for plants  
     > BiomassA   	#  total biomass of animals    
     >  ExtinctP   	# plants extincted     
     >  ExtinctA   	# animals     
     >  SurvivingP   # plants species surviving     
     >  SurvivingA   # animals     
     >  NestednessPin 	#  Nestedness of plants matrix (defined as in Bastolla, Nature 2009, i.e. ecological overlap)
     >  NestednessAin	#   Nestedness of animals     
     >  NestednessTin	#   Total nestedness     
     >  NestednessPout	#   Nestedness of plants species after the simulation, it changes just if there were extinctions     
     >  NestednessAout	#   animals     
     >  NestednessTout	#   total     
     >  ConnectanceIn 	#   Starting connectance 
     >  ConnectanceOut 	#   Output connectance
     >  SinglesIn   	# Starting number of species with a single connection    
     >  SinglesOut 	#   Final number    
     >  SingleExtinct 	#   Species with a single connection that got extincted    
     >  excludedPin 	#   Starting plant species forming an independent module
     >  excludedPout 	#  Final
     >  excludedAin 	#  Same for animals
     >  excludedAout           
     >  NegAlphaP 	#           Number of negative growth rates for plants (good control when growth rates are computed at the fixed point)
     >  NegAlphaA 	#           for animals
     >  Dissipation	#   Several measures of dissipation of the system (under research)     
     >  AvDissipation 	#  ""     
     >  DissRate 	#  ""     
     >  AvDissRate 	#  ""
     >  AvAlphaP 	#   Average of the growth rates of plant species     
     >  VarAlphaP 	#  Variance    
     >  AvAlphaA  	#  Average of the growth rates of animal species          
     >  VarAlphaA 	#  Variance   
     >  AvDegreeInP 	#   Average of the degree of plant species     
     >  VarDegreeInP 	#   Variance    
     >  AvDegreeInA  	#    Average of the degree of animal species 
     >  VarDegreeInA 	#   Variance     
     >  BiomassRate0  	#  Ratio between average plant and animal abundances (very relevant in obligatory mutualism)     
     >  MaxhP  	#  Maximum handling time for plants (very relevant in obligatory mutualism)     
     >  MaxhA  	#  animals   
     >  MaxRhoP	#     Maximum interspecific competition for plants (another possibility to build obligatory mutualism not explored)
     >  MaxRhoA	#     For animals
     >  Seed	#     Starting seed (control it changes along realizations)     
     >  EntropyP 	#    Effective biodiversity of plants     
     >  EntropyA 	#    Of animals
     
 * The script will integrate your system without perturbing it if you build a system that it is not at a fixed point. You can additionally perturb the growth rates, controlled in the .in file at the Delta parameter. This should be typically a value between 0 and 1 representing the amount of perturbation relative to the growth rates (so 1 represents a perturbation of 100% of the growth rates, and you will get extinctions). The script also allows you to build a system at a fix point, to perturb its growth rates next, which is the basis of our approach to structural stability (see Pascual-García and Bastolla, Nat Comm 2017). The idea is that you can fix all the parameters and solve the system of equations at the fix point leaving the growth rates free, that will be fixed with this procedure. To do this, you should provide a negative value for the growth rates at the .in file.   


* COMPILATION:
  Compiled in x86-64bits under Linux with gnu Fortran 77 with the following options:
	>  g77 -O2 -ffixed-line-length-0 -march=nocona -fcase-upper -o Name.exe *.f (for old compilers)
        >  gfortran -O2 -ffixed-line-length-0 -march=nocona -o Name.exe *.f 


* EXECUTION:
    Once the files are compiled, just execute the file, ensuring that the necessary files are in the same folder. If the output files already exist the script will not overwrite them and it stops.
