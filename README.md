Beyond-MeanField.f
==================

* DESCRIPTION:
  * This script builds a bipartite mutualistic network and integrates the system of equations. The scripts allows you to build the system either providing metaparameter values (features of the distributions around which numbers  will be generated) or reading the parameters and species abundances from files. You can read some of the parameters/abundances and generate others, the only mandatory input is a network of mutualistic interactions between both pools, typially with plants in rows and animals in columns. If you have other pools keep in mind that, in the following, parameters are distinguished as plants and animals, and one should be consistent.
  * If you provide the input files for the parameters or for the species abundances, the files should be named and formatted follows. You can change these names formatting the code in the file 1Beyond-MeanField_2-3.f.
    * 'plantsIn.dat' Plants abundances (a column vector)         
    * 'animalsIn.dat' Animal abundances
    * 'betaInP.dat' Competition for plants, a matrix (space separated with no header)
    * 'betaInA.dat' Competition for animals
    * 'gammaInP.dat' Mutualistic benefit of plants due to animals, a matrix (space separated with no header)      
    * 'gammaInA.dat' Mutualistic benefit of animals due to plants
    * 'alphaInP.dat' Growth rates of plants (a column vector)
    * 'alphaInA.dat' Growth rates of animals
    
  * You will also need a file called 'Beyond-MeanField.in' with the following lines (We provide some values as examples, and after '#' we show the name of the variable in the code and we comment the description):

     * 100	# (midNp) Plants, Population densities, a uniform distribution centered around this value is generated if >0, if =0 read from file.
     *  5	# (widthP) Width of the distribution, computed as:  Na(i)=midNa*(1+rnd*widthP), with rnd=[-0.5,0.5].
     * 	1	# (midNa) Animals, Population densities, center
     * 	0.05	# Animal, Width
     * 	-1	# (Alpha) Growth rates. If >0 generate random numbers, =0, read from file, <0 solve at the fixed point the equations after fixing from the other parameters, see below)
     * 	0.05	# Alpha, Width
     * 	1	# (Beta0) Intraspecific competition. If >0 random numbers, =0 read from file. It multiplies a distribution centered in one (i=j) or rho (i!=j) with the following width:
     * 	0.15	# Beta0, Width
     * 	0.23	# (rho) interspecific competition for plants (always constant, no noise here)
     * 	0.15	# (rho) interspecific competition for animals (always constant, no noise here)
     * 	0.15	# Mutualistic strengh: Gamma0<-100000 read from file, -10000<Gamma0<0 system with no benefit for plants, >0 mutualistic. Generates a Unif(Gamma0,Width).
     * 	0.15	# Gamma0, Width
     * 	0.1	# (hP) Handling time (Plants)
     * 	0.1	# (hA) Handling time (Animals)
     * 	50	# (Sp) Species number (Plants)
     * 	14	# (Sa) Species number (Animals)
     * 	0.5	# (Delta) amplitude of the perturbation over the growth rates
     * 	/home/user/AdjacencyMatrix.txt	# Path for the adjacency matrix (we recommend to use absolute paths)
     * 	/home/user/Summary_NameOutput.dat	# Absolute path for the output summary, with the name of the file you want to use   
     * 	0.6890  # This parameter is only available in the "multiple" version. It is a seed  needed when you run multiple times the script, and it should randomly change through realizations
	
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
     
 * The script will integrate your system without perturbing it if you build a system that it is not at a fixed point. You can additionally perturb the growth rates, controlled in the .in file at the Delta parameter. This should be typically a value between 0 and 1 representing the amount of perturbation relative to the growth rates (so 1 represents a perturbation of 100% of the growth rates, and you will get extinctions). The script also allows you to build a system at a fix point, to perturb its growth rates next, which is the basis of our approach to structural stability (see Pascual-García and Bastolla, Nat Comm 2017). The idea is that you can fix all the parameters and solve the system of equations at the fix point leaving the growth rates free, that will be fixed with this procedure. To do this, you should provide a negative value for the growth rates in the .in file.   


* COMPILATION: 
    To compile just type:
    $ make
 
    Compilation has bee tested in x86-64bits under Linux with older versions of gnu Fortran 
    and with gfortran (current version).The code is structured in three files labeled as 1$Description.f 2$Description.f 3$Description.f. You may find two or more
    versions of each of them, and you should select the version you need to build the executable. For instance, the files 2Parameters.f and 2Parameters_logNorm.f, read
    and generate parameters for the model, but the former generates equilibrium abundances following a uniform distribution and the latter a lognormal distribution, you
    should change the Makefile accordingly.

* EXECUTION:
    Once the files are compiled, just execute the file, ensuring that the necessary files are in the same folder. If the output files already exist the script will not overwrite them and it stops.

* TREE
  Structure of the code/files in which the functions are located:

* ... root (Main):
*
*     ..... Read_Parameters .... Header                 ... When  .. idate                                   ##(At 1Beyond-MeanField_2-3.f)
*                                                                 .. itime
*                           .... Warning
*
*     ..... Pseudo_Main     .... Read_Generate_Inputs   ... ran2                                             ##(At 2Parameters.f)
*                                                       ... Topology       
*                                                       ... Consistent_Alpha .. ran2 
*
*                           .... Integration            ... Derivatives                                      ##(At 3Integration.f)     
*                                                       ... bsstep           .. mmid      . Derivatives
*                                                                            .. pzextr
*                           .... Output_Globals         ... Topology                                         ##(back to 1Beyond-MeanField_2-3.f)

