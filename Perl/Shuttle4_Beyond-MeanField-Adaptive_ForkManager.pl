#!/usr/bin/perl -w
#
# ****************************************
#  Multiple Shuttle for Beyond-MeanField.f
# ****************************************
#
# Multiple Shuttle for the program Beyond-MeanField.f. Our aim is to get values of biomass and
# extinctions for different values of a selected parameter, and
# Delta, which controls the fluctuations over the productivities vector due to environmental
# fluctuations. We will run Nrnd times each set of parameters obtaining the average, standard
# deviation and an error estimation of all outputs provided by the Fortran script in
# a single file.

# -- Additional requirements

use POSIX ":sys_wait_h"; # libraries to create child processes with fork (not needed if you use fork manager)
require("syscall.ph"); # allow system calls
use Parallel::ForkManager; # libraries to manage fork.

####### START EDITING HERE  ##########

# --Number of species and randomizations:

$Sp=50;	# Species number (Plants)
$Sa=50;	# Species number (Animals)
$Nrnd=500;   # Number of realizations (calls to F77), for each set of parameters.

# -- Variable parameters
# .... Set here the variable you want to explore dynamically, and the values for Delta

@Delta=(0.1,0.25,0.5,0.75); # Fix here the values you want to explore for Delta

# .... Possible variables to explore dynamically are "midGamma","midRho","midBeta","midH" and "midG".
# .... Note that the same value will be given for plants and animals.

$selectVar="midGamma"; # A string for the selected variable
@selectValues=(0.02,0.04,0.06); # the values to be explored, will override any value fixed below for the selected variable

# -- Labels for the output files
$Set='Mutualism-Test'; # Name of the system (label for the output files)
$param='Facultative_GammaVariable'; # Label to identify the regime of parameters (label for the output files)

# -- Additional parameters for the .in file of F77 program (Beyond-MeanField.f) A

$readNp=0; #  Control parameter for plants abundances, it is equal to 1 if you read Np from file, 0 if it will be generated internally, in which case:
$midNp =1; # It will generate the parameter with a distribution of plant densities Np=midNp *(1+ widthNp *rnd) where rnd = [-0.5,0.5]
$widthNp =0.15; #  Note that these two parameters must be given even if readNp=1
$readNa =0; # Same for animal densities 
$midNa =1; #  
$widthNa =0.15; # 
$readAlphaP =0; #  Plants growth rates (read or generate?)
$midAlphaP =1; #  
$widthAlphaP =0.15; # 
$readAlphaA =0; #  Animals growth rates
$midAlphaA =1; # 
$widthAlphaA=0.15; #
$readBetaP =0; #  Sp x Sp competition between plants, equal to 1 if the matrix is readed from file, otherwise:
$midBetaP =1; # First parameter to generate the intraspecific competition for plants (diagonal terms)
$widthBetaP =0.15; # Width of the distribution, again: BetaP=midBetaP *(1+ widthBetaP *rnd)
$midRhoP =0.05; # Second parameter to generate the interespecific competition for animals (off-diagonal terms)
$widthRhoP =0.15; #  Width of the distribution
$readBetaA =0; # Same for for competition between animals
$midBetaA =1; #  Intraspecific competition for animals
$widthBetaA =0.15; # Width of the distribution
$midRhoA =0.05; # interspecific competition (animals)
$widthRhoA =0.15; # Width of the distribution
$readGammaP =0; # Sp x Sa interaction between the pool of plants and the pool of animals, equal to 1 if the matrix is readed from file, otherwise:
$midGammaP =0.05; # Generate a matrix with cells given by GammaP=midGammaP *(1+ widthGammaP *rnd)
$widthGammaP =0.15; #
$readGammaA =0; # Same for the interaction between animals and plans
$midGammaA =0.05; # 
$widthGammaA =0.15; #
$readHp =0; # Sp x 1, vector of handling time for plants. Set to 1 if you read from file, 0 if must be generated internally with
$midHp =0.1; # Hp=midHp *(1+ widthHp *rnd) 
$widthHp =0.15; # 
$readHa =0; #  Sa x 1 vector of handling time for animals.
$midHa =0.1; # 
$widthHa =0.15; #
$readGp =0; # Sp x 1, second kind of vector of handling time for plants. Set to 1 if you read from file, 0 if must be generated internally with
$midGp =0; #  
$widthGp =0; # 
$readGa =0; # Sa x 1, vector of handling time for animals
$midGa =0; # 
$widthGa=0; # 
$f0P=1; #  Scale of the saturating term (real number, plants)
$f0A=1; # Scale of the saturating term (real number, animals)
$fixedPoint=1; # Should the growth rates (Alpha) estimated from the rest of the model at the fixed point? (=1) or not (=0). This option will override above parameters for Alpha. 
$outTmp='Summary_Beyond-MeanField.tmp';	# Path for the output summary matrix


# -- Parameters to control processors

$Nproc=`grep -c proc /proc/cpuinfo`; # Get the number of free processors of your computer, set to one if running sequentially (in a cluster with queues for instance)
$threads=$Nproc;   # Number of threads it will use for running FORTRAN. You can fix it to the maximum available ($Nproc) or to a fixed number.

####### STOP EDITING HERE  ##########

# --- Check if you need input files

if(($readNp==1)||($readNa==1)||($readBetaP==1)||($readBetaA==1)||($readGammaP==1)||($readGammaA==1)||($readHp==1)||($readHa==1)||($readGp==1)||($readGa==1)){
    $controlFiles=1;
}else{
    $controlFiles=0;
}

# --- Some reordering, taking sizes and fixing files that will be used
my $max = ($Sp, $Sa)[$Sp < $Sa];
$degreePtmp='plantsOut-Final.dat'; # From these files we will get the degrees, to compute the probability of extinction as a function of the degree
$degreeAtmp='animalsOut-Final.dat';
$length1=@selectValues;
$length2=@Delta;
for($k=0;$k<= $max;$k++){$DegVec[$k]=0;} # Initialize the counter for extincted species as a function of degree
$ThrExt=1e-08; # Biomass lower bound to consider an specie extincted
$Init=1; # Identifier of the first realization

# -- Global input and output files, and handles definition

$fileReadIn = 'Beyond-MeanField.in'; # The input file must be edit at each step
$fileout2 = 'Summary-Ensemble_Beyond-MeanField_'.$Set.'_'.$param.'.dat'; # This is the output file
$fileout3 = 'DevfromRefSimulation_Beyond-MeanField_'.$Set.'_'.$param.'.dat';
$fileout4 = 'ProbExtinctDegrees_Beyond-MeanField_'.$Set.'_'.$param.'.dat';

$Handle[0]=*HANDLE2; $Handle[1]=*HANDLE3; $Handle[2]=*HANDLE4;
open(HANDLE2, ">$fileout2") || die "Couldn't open file $fileout2"; # The output file can already be opened
open(HANDLE3, ">$fileout3") || die "Couldn't open file $fileout3"; # The output file can already be opened
open(HANDLE4, ">$fileout4") || die "Couldn't open file $fileout4"; # The output file can already be opened

# -- Build output headers

@date=localtime(time); # Some functions to construct the date for the header..
@months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
@weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($second, $minute, $hour, $dayOfMonth, $month, $yearOffset, $dayOfWeek, $dayOfYear, $daylightSavings) = localtime();
$year = 1900 + $yearOffset;
$theTime = "$hour:$minute:$second, $weekDays[$dayOfWeek] $months[$month] $dayOfMonth, $year";
print '  ',"\n";
print '****************************************',"\n";
print ' Multiple Shuttle for Beyond-MeanField.f ',"\n";
print '****************************************',"\n";
print '  ',"\n";
print '** Running multiple times Lotka-Volterra mutualistic simulations beyond mean field',"\n";
print '  ',"\n";
for($i=0;$i<=2;$i++) # for each output file
{
    $file=$Handle[$i];
    print $file '# >> Multiple Shuttle for Beyond-MeanField.f <<',"\n";
    print $file '# >> Biomasses and number of extincted species for different parameter values',"\n"; 
    print $file '# >> Running at date: ',$theTime,"\n";
    print $file '# >> Zurich, Theoretical Biology (ETH), apascualgarcia.github.io', "\n";
    print $file '# >> COMMON PARAMETERS:', "\n";
    print $file '# Population densities plants, readed from file?(=1): ',$readNp,'; midNp: ',$midNp,'; widthNp: ',$widthNp,"\n";
    print $file '# Population densities animals, readed from file?(=1): ',$readNa,'; midNa: ',$midNa,'; widthNa: ',$widthNa,"\n";
    print $file '# Growth rates plants, estimated at fixed point?(=1): ',$fixedPoint,'; readed from file?(=1): ',$readAlphaP,'; midNp: ',$midAlphaP,'; widthNp: ',$widthAlphaP,"\n";
    print $file '# Growth rates animals, estimated at fixed point?(=1): ',$fixedPoint,'; readed from file?(=1): ',$readAlphaA,'; midAlphaA: ',$midAlphaA,'; widthAlphaA: ',$widthAlphaA,"\n";
    print $file '# Beta parameter plants, readed from file?(=1): ',$readBetaP,'; midBetaP: ',$midBetaP,'; widthBetaP: ',$widthBetaP,"\n";
    print $file '# Beta parameter animals, readed from file?(=1): ',$readBetaA,'; midBetaA: ',$midBetaA,'; widthBetaA: ',$widthBetaA,"\n";
    print $file '# Rho parameter plants, readed from file?(=1): ',$readBetaP,'; midRhoP: ',$midRhoP,'; widthRhoP: ',$widthRhoP,"\n";
    print $file '# Rho parameter animals, readed from file?(=1): ',$readBetaA,'; midRhoA: ',$midRhoA,'; widthRhoA: ',$widthRhoA,"\n";
    print $file '# Gamma parameter plants, readed from file?(=1): ',$readGammaP,'; midGammaP: ',$midGammaP,'; widthGammaP: ',$widthGammaP,"\n";
    print $file '# Gamma parameter animals, readed from file?(=1): ',$readGammaA,'; midGammaA: ',$midGammaA,'; widthGammaA: ',$widthGammaA,"\n";
    print $file '# H parameter plants, readed from file?(=1): ',$readHp,'; midHp: ',$midHp,'; widthHp: ',$widthHp,"\n";
    print $file '# H parameter animals, readed from file?(=1): ',$readHa,'; midHa: ',$midHa,'; widthHa: ',$widthHa,"\n";
    print $file '# G parameter plants, readed from file?(=1): ',$readGp,'; midGp: ',$midGp,'; widthGp: ',$widthGp,"\n";
    print $file '# G parameter animals, readed from file?(=1): ',$readGa,'; midGa: ',$midGa,'; widthGa: ',$widthGa,"\n";
    print $file '# Species number (Plants): ',$Sp,"\n";
    print $file '# Species number (Animals): ',$Sa,"\n";
    print $file '# Path for the output summary matrix: ',$outTmp,"\n";
    print $file '# Number of realizations for each set of parameters: ', $Nrnd,"\n";
    print $file '# Values for the Delta parameter: ',join("; ",@Delta),"\n";
    print $file '# Values for the ',$selectVar,' parameter: ',join("; ",@selectValues),"\n";

}

# -- MAIN ROUTINE


srand(); # Changes the seed every time you run the script, comment if you want to reproduce results.
for($i=0;$i<=$length1-1;$i++)    
{
    if($selectVar eq "midGamma"){
	$midGammaP=$selectValues[$i];
	$midGammaA=$selectValues[$i];
    }elsif($selectVar eq "midBeta"){
	$midBetaP=$selectValues[$i];
	$midBetaA=$selectValues[$i];
    }elsif($selectVar eq "midRho"){
	$midRhoP=$selectValues[$i];
	$midRhoA=$selectValues[$i];
    }elsif($selectVar eq "midH"){
	$midHp=$selectValues[$i];
	$midHa=$selectValues[$i];
    }elsif($selectVar eq "midG"){
	$midGp=$selectValues[$i];
	$midGa=$selectValues[$i];
    }
    for($j=0;$j<=$length2-1;$j++)
    {
	print ' >> Computing realizations for parameters ',$selectVar,' = ',$selectValues[$i],' Delta = ',$Delta[$j],"\n";
	my $pm = new Parallel::ForkManager($threads); # Set a new environment to start to run different threads
	for($m=$Init;$m<=$Nrnd;$m++) # For each matrix you compute the score
	{
	    my $random = rand();
	    $pm->start and next; # Start the threads
	    # Create input file for F77 program
	    
	    $dir='dir_tmp'.$m; # Create directories to generate different .in files for FORTRAN calls
	    if(($i == 0)&&($j == 0)){`mkdir  $dir`;}
	    if($controlFiles==1){
		system("cp *In*.dat $dir"); # Copy input files if you are gonna use them
	    }
	    chdir $dir; # Call to FORTRAN binary where the .in files are
	    
            # Edit the .in file for FORTRAN
	    open(HANDLE1, ">$fileReadIn") || die "Couldn't open file $fileReadIn";
	    print HANDLE1 $readNp," readNp \n";
	    print HANDLE1 $midNp ," midNp\n";
	    print HANDLE1 $widthNp ," widthNp\n";
	    print HANDLE1 $readNa," readNa\n"; 
	    print HANDLE1 $midNa," midNa\n"; 
	    print HANDLE1 $widthNa," widthNa\n"; 
	    print HANDLE1 $readAlphaP," readAlpha\n"; 
	    print HANDLE1 $midAlphaP," midAlphaP\n"; 
	    print HANDLE1 $widthAlphaP," widthAlphaP\n"; 
	    print HANDLE1 $readAlphaA ," readAlphaA\n";
	    print HANDLE1 $midAlphaA," midAlphaA\n"; 
	    print HANDLE1 $widthAlphaA," widthAlphaA\n";
	    print HANDLE1 $readBetaP," readBeta\n"; 
	    print HANDLE1 $midBetaP ," midBetaP\n";
	    print HANDLE1 $widthBetaP ," widthBetaP\n";
	    print HANDLE1 $midRhoP,"1 $midRhoP\n"; 
	    print HANDLE1 $widthRhoP," widthRhoP\n"; 
	    print HANDLE1 $readBetaA," readBetaA\n"; 
	    print HANDLE1 $midBetaA," midBetaA\n"; 
	    print HANDLE1 $widthBetaA," widthBeta\n"; 
	    print HANDLE1 $midRhoA," midRhoA\n"; 
	    print HANDLE1 $widthRhoA," widthRhoA\n"; 
	    print HANDLE1 $readGammaP," readGammaP\n"; 
	    print HANDLE1 $midGammaP," midGammaP\n"; 
	    print HANDLE1 $widthGammaP," widthGammaP\n"; 
	    print HANDLE1 $readGammaA," readGammaA\n"; 
	    print HANDLE1 $midGammaA," midGammaA\n"; 
	    print HANDLE1 $widthGammaA," widthGammaA\n"; 
	    print HANDLE1 $readHp," readHp\n"; 
	    print HANDLE1 $midHp," midHp\n"; 
	    print HANDLE1 $widthHp," widthHp\n"; 
	    print HANDLE1 $readHa," readHa\n"; 
	    print HANDLE1 $midHa," midHa\n"; 
	    print HANDLE1 $widthHa," widthHa\n"; 
	    print HANDLE1 $readGp," readGp\n"; 
	    print HANDLE1 $midGp," midGp\n"; 
	    print HANDLE1 $widthGp," widthGp\n"; 
	    print HANDLE1 $readGa," readGa\n"; 
	    print HANDLE1 $midGa," midGa,\n"; 
	    print HANDLE1 $widthGa," widthGa\n";
	    print HANDLE1 $f0P," f0P\n";
	    print HANDLE1 $f0A," f0A\n";
	    print HANDLE1 $Sp," Sp\n";
	    print HANDLE1 $Sa," Sa\n";
	    print HANDLE1 $Delta[$j]," Delta[j]\n";
	    print HANDLE1 $fixedPoint," fixedPoint\n";
	    print HANDLE1 $outTmp," outTmp\n";
	    print HANDLE1 $random," random\n";
	    close(HANDLE1);


            ######## CALLING F77 PROGRAM #######  
	    system("../BeyondMeanField_gfortran > /dev/null");                         
	    chdir "..";        
	    print ' << Realization ',$m,' launched!',"\n";

	    $pm->finish; # End of threads  
	}
	$pm->wait_all_children; # Wait until all the randomizations are finished
	#exit; # DEBUG
	for($m=$Init;$m<=$Nrnd;$m++) # For each realization
	{
	    $dir='dir_tmp'.$m; # Go the directory to recover results
	    chdir $dir;
	    # Store values in variables
	    open(HANDLEtmp, $outTmp)|| die "Couldn't open file $outTmp";
	    $k=-1; # Control outputs
	    $l=-1; # Control lines
	    while($line=<HANDLEtmp>)
	    {
		if($line =~ /#/){next;}
		@fields=split(/\s+/,$line);
		$Id=$fields[1];
		$Value=$fields[2];
		$l=$l+1;	      
		#print $m,' m ',$k,' k ', $fields[0],'--',$fields[1],'--',$Value,"\n"; #debug
		#exit; #debug
		#if($l < 17){next;}else{$k=$k+1;} # There are 15 input parameters we are not interested in
		if(($i == 0)&&($j == 0)&&($m == $Init)){
		    if($Id eq "NestednessPin"){
			$NestednessP=$Value;
		    }elsif($Id eq "NestednessAin"){			
			$NestednessA=$Value;
		    }elsif($Id eq "NestednessTin"){			
			$NestednessT=$Value;
		    }elsif($Id eq "ConnectanceIn"){			
			$Connectance=$Value;
		    }elsif($Id eq "SinglesOut"){
			$Singles=$Value;
		    }elsif($Id eq "AvDegreeInP"){			
			$AvDegreeP=$Value;
		    }elsif($Id eq "VarDegreeInP"){			
			$VarDegreeP=$Value;
		    }elsif($Id eq "AvDegreeInA"){			
			$AvDegreeA=$Value;
		    }elsif($Id eq "VarDegreeInA"){			
			$VarDegreeA=$Value;
			# Now that you have information common to all simulations, finish the headers
			for($ii=0;$ii<=2;$ii++){ 
			    $file=$Handle[$ii];
			    print $file '# Input nestedness for plants: ',$NestednessP,"\n";
			    print $file '# Input nestedness for animals: ',$NestednessA,"\n";
			    print $file '# Input total nestedness: ',$NestednessT,"\n";
			    print $file '# Input connectance: ',$Connectance,"\n";
			    print $file '# Average degree for plants: ',$AvDegreeP,"\n";
			    print $file '# Variance degree for plants: ',$VarDegreeP,"\n";
			    print $file '# Average degree for animals: ',$AvDegreeA,"\n";
			    print $file '# Variance degree for animals: ',$VarDegreeA,"\n";
			    print $file '# Number of starting singletons: ',$Singles,"\n";
			    if($ii ==1){ # for the first file $i == 0, we will introduce this information later, because there are some previous measures to compute and include.			     
				print $file '#1',$selectVar,' 2Delta, 3dBioP, 4Err_dBioP, 5dBioA, 6Err_dBioA, 7dSp, 8Err_dSp, 9dSa, 10Err_dSa, 11dSpRel, 12Err_dSpRel, 13dSaRel, 14Err_dSaRel, 15dSrel, 16Err_dSrel, 17dNestOutP, 18Err_dNestOutP, 19dNestOutA, 20Err_dNestOutA, 21dNestOutT, 22Err_dNestOutT, 23dConnectOut, 24Err_dConnectOut, 25dSinglesOut, 26Err_dSinglesOut, 27dSinglesExtinct, 28Err_dSinglesExtinct, 29dExcludeP, 30Err_dExcludeP, 31dExcludeA, 32Err_dExcludeA', "\n";
			    }elsif($ii==2){			     
				print $file '#1',$selectVar,' 2Delta, 3Degree, 4ProbExtinction(degree), 5TotalExtinctions', "\n";
			    }
			}
		    }
		} # end of if(($i == 0)&&($j == 0))&&($m == $Init))
		# Now fill the variables that change through randomizations
		if($Id eq "BiomassP"){		
		    if($m == $Init){ # Initialize		
			$BiomassP[$i][$j]=$Value; $StdvBioP[$i][$j]=$Value**2;
		    }else{
			$BiomassP[$i][$j]+=$Value; $StdvBioP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "BiomassA"){
		    if($m == $Init){ # Initialize							
			$BiomassA[$i][$j]=$Value; $StdvBioA[$i][$j]=$Value**2;
		    }else{
			$BiomassA[$i][$j]+=$Value; $StdvBioA[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "ExtinctP"){
		    if($m == $Init){ # Initialize							
			$ExtinctP[$i][$j]=$Value; $StdvExtP[$i][$j]=$Value**2;
		    }else{
			$ExtinctP[$i][$j]+=$Value; $StdvExtP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "ExtinctA"){
		    if($m == $Init){ # Initialize							
			$ExtinctA[$i][$j]=$Value; $StdvExtA[$i][$j]=$Value**2;
		    }else{
			$ExtinctA[$i][$j]+=$Value; $StdvExtA[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "SurvivingP"){			    
		    if($m == $Init){ # Initialize							
			$SurviveP[$i][$j]=$Value; $StdvSurvP[$i][$j]=$Value**2;
		    }else{
			$SurviveP[$i][$j]+=$Value; $StdvSurvP[$i][$j]+=$Value**2;	
		    }
		}elsif($Id eq "SurvivingA"){
		    if($m == $Init){ # Initialize							
			$SurviveA[$i][$j]=$Value; $StdvSurvA[$i][$j]=$Value**2;		
		    }else{
			$SurviveA[$i][$j]+=$Value; $StdvSurvA[$i][$j]+=$Value**2;			
		    }	
		}elsif($Id eq "NestednessPout"){
		    if($m == $Init){ # Initialize							
			$NestOutP[$i][$j]=$Value; $StdvNestOutP[$i][$j]=$Value**2;
		    }else{
			$NestOutP[$i][$j]+=$Value; $StdvNestOutP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "NestednessAout"){
		    if($m == $Init){ # Initialize										
			$NestOutA[$i][$j]=$Value; $StdvNestOutA[$i][$j]=$Value**2;
		    }else{
			$NestOutA[$i][$j]+=$Value; $StdvNestOutA[$i][$j]+=$Value**2;	
		    }		
		}elsif($Id eq "NestednessTout"){
		    if($m == $Init){ # Initialize										
			$NestOutT[$i][$j]=$Value; $StdvNestOutT[$i][$j]=$Value**2;
		    }else{
			$NestOutT[$i][$j]+=$Value; $StdvNestOutT[$i][$j]+=$Value**2;
		    }		    
		}elsif($Id eq "ConnectanceOut"){
		    if($m == $Init){ # Initialize										
			$ConnectOut[$i][$j]=$Value; $StdvConnectOut[$i][$j]=$Value**2;
		    }else{
			$ConnectOut[$i][$j]+=$Value; $StdvConnectOut[$i][$j]+=$Value**2;	
		    }
		}elsif($Id eq "SinglesOut"){
		    if($m == $Init){ # Initialize										
			$SinglesOut[$i][$j]=$Value; $StdvSinglesOut[$i][$j]=$Value**2;
		    }else{
			$SinglesOut[$i][$j]+=$Value; $StdvSinglesOut[$i][$j]+=$Value**2;	
		    }	       
		}elsif($Id eq "SingleExtinct"){
		    if($m == $Init){ # Initialize										
			$SinglesExtinct[$i][$j]=$Value; $StdvSinglesExtinct[$i][$j]=$Value**2;
		    }else{
			$SinglesExtinct[$i][$j]+=$Value; $StdvSinglesExtinct[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "excludedPin"){ # ix, it is just a control as there may not exist
		    if($m == $Init){ # Initialize										
			$ExcludeP[$i][$j]=$Value; $StdvExcludeP[$i][$j]=$Value**2;
		    }else{
			$ExcludeP[$i][$j]+=$Value; $StdvExcludeP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "excludedAin"){ # ix, it is just a control as there may not exist
		    if($m == $Init){ # Initialize										
			$ExcludeA[$i][$j]=$Value; $StdvExcludeA[$i][$j]=$Value**2;
		    }else{
			$ExcludeA[$i][$j]+=$Value; $StdvExcludeA[$i][$j]+=$Value**2;
		    }		    
		}elsif($Id eq "NegAlphaP"){
		    if($m == $Init){ # Initialize										
			$NegAlphaP[$i][$j]=$Value; $StdvNegAlphaP[$i][$j]=$Value**2;
		    }else{
			$NegAlphaP[$i][$j]+=$Value; $StdvNegAlphaP[$i][$j]+=$Value**2;						
		    }
		}elsif($Id eq "NegAlphaA"){
		    if($m == $Init){ # Initialize										
			$NegAlphaA[$i][$j]=$Value; $StdvNegAlphaA[$i][$j]=$Value**2;
		    }else{
			$NegAlphaA[$i][$j]+=$Value; $StdvNegAlphaA[$i][$j]+=$Value**2;						
		    }
		}elsif($Id eq "AvAlphaP"){
		    if($m == $Init){ # Initialize		    
			$AvAlphaAvP[$i][$j]=$Value; $StdvAlphaAvP[$i][$j]=$Value**2;
		    }else{
			$AvAlphaAvP[$i][$j]+=$Value; $StdvAlphaAvP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "VarAlphaP"){
		    if($m == $Init){ # Initialize		    
			$AvAlphaVarP[$i][$j]=$Value; $StdvAlphaVarP[$i][$j]=$Value**2;
		    }else{
			$AvAlphaVarP[$i][$j]+=$Value; $StdvAlphaVarP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "AvAlphaA"){
		    if($m == $Init){ # Initialize		    
			$AvAlphaAvA[$i][$j]=$Value; $StdvAlphaAvA[$i][$j]=$Value**2;
		    }else{
			$AvAlphaAvA[$i][$j]+=$Value; $StdvAlphaAvA[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "VarAlphaA"){
		    if($m == $Init){ # Initialize		    
			$AvAlphaVarA[$i][$j]=$Value; $StdvAlphaVarA[$i][$j]=$Value**2;
		    }else{
			$AvAlphaVarA[$i][$j]+=$Value; $StdvAlphaVarA[$i][$j]+=$Value**2;
		    }		    
		}elsif($Id eq "TrueDiversityP"){
		    if($m == $Init){ # Initialize		    
			$entropyP[$i][$j]=$Value; $StdvEntropyP[$i][$j]=$Value**2;
		    }else{
			$entropyP[$i][$j]+=$Value; $StdvEntropyP[$i][$j]+=$Value**2;
		    }
		}elsif($Id eq "TrueDiversityA"){
		    if($m == $Init){ # Initialize		    
			$entropyA[$i][$j]=$Value; $StdvEntropyA[$i][$j]=$Value**2;
		    }else{
			$entropyA[$i][$j]+=$Value; $StdvEntropyA[$i][$j]+=$Value**2;
		    }
		}else{next;}
	    } #end while($line=<HANDLEtmp>)
	    close(HANDLEtmp);
	    #exit; #debug

	    # Relate the extinctions with the degrees
	    open(HANDLEtmp, $degreePtmp)|| die "Couldn't open file $degreePtmp";	    
	    while($line=<HANDLEtmp>)
	    {
		if($line =~ /#/){next;}
		@fields=split(/\s+/,$line);
		$BioTmp=$fields[3]; $DegTmp=$fields[5];
		#print join(" ",$fields[0],'0',$fields[1],'1',$fields[2],'2',$fields[3],'3',,$fields[4],'4',$fields[5],'5','DEBUG0'),"\n";
		if($BioTmp < $ThrExt){$DegVec[$DegTmp]+=1;}
	    }
	    #exit; # debug
	    close(HANDLEtmp);
	    open(HANDLEtmp, $degreeAtmp)|| die "Couldn't open file $degreeAtmp";	    
	    while($line=<HANDLEtmp>)
	    {
		if($line =~ /#/){next;}
		@fields=split(/\s+/,$line);
		$BioTmp=$fields[3]; $DegTmp=$fields[5];
		if($BioTmp < $ThrExt){$DegVec[$DegTmp]+=1;}
	    }
	    close(HANDLEtmp);
	    `rm -f  $outTmp *.dat *.in`;
	    chdir "..";
	} # end for($m=$Init;$m<=$Nrnd;$m++)
	#exit; # debug

	# --- Compute Measures from Summary file
	#print join(" ",$i,$j,$BiomassP[$i][$j],$StdvBioP[$i][$j],$BiomassA[$i][$j],$StdvBioA[$i][$j],'DEBUG0'),"\n";
	$BiomassP[$i][$j]=$BiomassP[$i][$j]/$Nrnd; $BiomassA[$i][$j]=$BiomassA[$i][$j]/$Nrnd;
	$StdvBioP[$i][$j]=$StdvBioP[$i][$j]/$Nrnd-$BiomassP[$i][$j]**2; $StdvBioA[$i][$j]=$StdvBioA[$i][$j]/$Nrnd-$BiomassA[$i][$j]**2; 
	#print join(" ",$i,$j,$BiomassP[$i][$j],$StdvBioP[$i][$j],$BiomassA[$i][$j],$StdvBioA[$i][$j],'DEBUG1'),"\n";
	$ExtinctAll=$ExtinctP[$i][$j]+$ExtinctA[$i][$j];  # To normalize degrees extinction later
	$ExtinctP[$i][$j]=$ExtinctP[$i][$j]/$Nrnd; $ExtinctA[$i][$j]=$ExtinctA[$i][$j]/$Nrnd;
	#print join(" ",$i,$j,$ExtinctP[$i][$j],$ExtinctA[$i][$j],$ExtinctAll,'DEBUG0'),"\n";
	#exit; # debug
	$StdvExtP[$i][$j]=$StdvExtP[$i][$j]/$Nrnd-$ExtinctP[$i][$j]**2; if($StdvExtP[$i][$j]<0){$StdvExtP[$i][$j]=0;} # This may happen if there are not extinctions due to numerical precision
	$StdvExtA[$i][$j]=$StdvExtA[$i][$j]/$Nrnd-$ExtinctA[$i][$j]**2; if($StdvExtA[$i][$j]<0){$StdvExtA[$i][$j]=0;}	
	$SurviveP[$i][$j]=$SurviveP[$i][$j]/$Nrnd; $SurviveA[$i][$j]=$SurviveA[$i][$j]/$Nrnd;
	$StdvSurvP[$i][$j]=$StdvSurvP[$i][$j]/$Nrnd-$SurviveP[$i][$j]**2; 
	$StdvSurvA[$i][$j]=$StdvSurvA[$i][$j]/$Nrnd-$SurviveA[$i][$j]**2;
	$NestOutP[$i][$j]=$NestOutP[$i][$j]/$Nrnd; $NestOutA[$i][$j]=$NestOutA[$i][$j]/$Nrnd; $NestOutT[$i][$j]=$NestOutT[$i][$j]/$Nrnd;
	$StdvNestOutP[$i][$j]=$StdvNestOutP[$i][$j]/$Nrnd-$NestOutP[$i][$j]**2; if($StdvNestOutP[$i][$j]<0){$StdvNestOutP[$i][$j]=0;}
	$StdvNestOutA[$i][$j]=$StdvNestOutA[$i][$j]/$Nrnd-$NestOutA[$i][$j]**2; if($StdvNestOutA[$i][$j]<0){$StdvNestOutA[$i][$j]=0;} 
	$StdvNestOutT[$i][$j]=$StdvNestOutT[$i][$j]/$Nrnd-$NestOutT[$i][$j]**2; if($StdvNestOutT[$i][$j]<0){$StdvNestOutT[$i][$j]=0;}
	$ConnectOut[$i][$j]=$ConnectOut[$i][$j]/$Nrnd;
	$StdvConnectOut[$i][$j]=$StdvConnectOut[$i][$j]/$Nrnd-$ConnectOut[$i][$j]**2; if($StdvConnectOut[$i][$j]<0){$StdvConnectOut[$i][$j]=0;}
	$SinglesOut[$i][$j]=$SinglesOut[$i][$j]/$Nrnd;
	$StdvSinglesOut[$i][$j]=$StdvSinglesOut[$i][$j]/$Nrnd-$SinglesOut[$i][$j]**2;
	$SinglesExtinct[$i][$j]=$SinglesExtinct[$i][$j]/$Nrnd;
	$StdvSinglesExtinct[$i][$j]=$StdvSinglesExtinct[$i][$j]/$Nrnd-$SinglesExtinct[$i][$j]**2;
	$ExcludeP[$i][$j]=$ExcludeP[$i][$j]/$Nrnd; $ExcludeA[$i][$j]=$ExcludeA[$i][$j]/$Nrnd;
	$StdvExcludeP[$i][$j]=$StdvExcludeP[$i][$j]/$Nrnd-$ExcludeP[$i][$j]**2; $StdvExcludeA[$i][$j]=$StdvExcludeA[$i][$j]/$Nrnd-$ExcludeA[$i][$j]**2;
	$NegAlphaP[$i][$j]=$NegAlphaP[$i][$j]/$Nrnd; $NegAlphaA[$i][$j]=$NegAlphaA[$i][$j]/$Nrnd;
	$StdvNegAlphaP[$i][$j]=$StdvNegAlphaP[$i][$j]/$Nrnd-$NegAlphaP[$i][$j]**2; $StdvNegAlphaA[$i][$j]=$StdvNegAlphaA[$i][$j]/$Nrnd-$NegAlphaA[$i][$j]**2;
	$AvAlphaAvP[$i][$j]=$AvAlphaAvP[$i][$j]/$Nrnd; $AvAlphaAvA[$i][$j]=$AvAlphaAvA[$i][$j]/$Nrnd; 
	$StdvAlphaAvP[$i][$j]= $StdvAlphaAvP[$i][$j]/$Nrnd-$AvAlphaAvP[$i][$j]**2; $StdvAlphaAvA[$i][$j]= $StdvAlphaAvA[$i][$j]/$Nrnd-$AvAlphaAvA[$i][$j]**2;
	$AvAlphaVarP[$i][$j]=$AvAlphaVarP[$i][$j]/$Nrnd; $AvAlphaVarA[$i][$j]=$AvAlphaVarA[$i][$j]/$Nrnd; 
	$StdvAlphaVarP[$i][$j]= $StdvAlphaVarP[$i][$j]/$Nrnd-$AvAlphaVarP[$i][$j]**2; $StdvAlphaVarA[$i][$j]= $StdvAlphaVarA[$i][$j]/$Nrnd-$AvAlphaVarA[$i][$j]**2;
	$entropyA[$i][$j]=$entropyA[$i][$j]/$Nrnd;
	$StdvEntropyA[$i][$j]=$StdvEntropyA[$i][$j]/$Nrnd-$entropyA[$i][$j]**2;
	$entropyP[$i][$j]=$entropyP[$i][$j]/$Nrnd;
	$StdvEntropyP[$i][$j]=$StdvEntropyP[$i][$j]/$Nrnd-$entropyP[$i][$j]**2;
	if(($i == 0)&&($j == 0))
	{
	    print HANDLE2 '#1',$selectVar,' 2Delta, 3BiomassP, 4StdvBioP, 5BiomassA, 6StdvBioA, 7ExtinctP, 8StdvExtP, 9ExtinctA, 10StdvExtA, 11SurviveP, 12StdvSurvP, 13SurviveA, 14StdvSurvA, 15NestOutP, 16StdvNestOutP, 17NestOutA, 18StdvNestOutA, 19NestOutT, 20StdvNestOutT, 21ConnectOut, 22StdvConnectOut, 23SinglesOut, 24StdvSinglesOut, 25SinglesExtinct, 26StdvSinglesExtinct, 27ExcludeP, 28StdvExcludeP, 29ExcludeA, 30StdvExcludeA, 31NegAlphaP, 32StdvNegAlphaP, 33NegAlphaA, 34StdvNegAlphaA, 35AvAlphaAvP, 36StdvAlphaAvP, 37AvAlphaAvA, 38StdvAlphaAvA, 39AvAlphaVarP, 40StdvAlphaVarP, 41AvAlphaVarA, 42StdvAlphaVarA, 43AvTrueDiversP, 52StdvTrueDiversP, 53AvTrueDiversA, 54StdvTrueDiversA', "\n";
	}		
	print HANDLE2 join(" ",$selectValues[$i],$Delta[$j],$BiomassP[$i][$j],$StdvBioP[$i][$j],$BiomassA[$i][$j],$StdvBioA[$i][$j],$ExtinctP[$i][$j],$StdvExtP[$i][$j],$ExtinctA[$i][$j],$StdvExtA[$i][$j],$SurviveP[$i][$j],$StdvSurvP[$i][$j],$SurviveA[$i][$j],$StdvSurvA[$i][$j],$NestOutP[$i][$j],$StdvNestOutP[$i][$j],$NestOutA[$i][$j],$StdvNestOutA[$i][$j],$NestOutT[$i][$j],$StdvNestOutT[$i][$j],$ConnectOut[$i][$j],$StdvConnectOut[$i][$j],$SinglesOut[$i][$j],$StdvSinglesOut[$i][$j],$SinglesExtinct[$i][$j],$StdvSinglesExtinct[$i][$j],$ExcludeP[$i][$j],$StdvExcludeP[$i][$j],$ExcludeA[$i][$j],$StdvExcludeA[$i][$j],$NegAlphaP[$i][$j],$StdvNegAlphaP[$i][$j],$NegAlphaA[$i][$j],$StdvNegAlphaA[$i][$j],$AvAlphaAvP[$i][$j], $StdvAlphaAvP[$i][$j], $AvAlphaAvA[$i][$j], $StdvAlphaAvA[$i][$j], $AvAlphaVarP[$i][$j], $StdvAlphaVarP[$i][$j], $AvAlphaVarA[$i][$j], $StdvAlphaVarA[$i][$j],$entropyA[$i][$j],$StdvEntropyA[$i][$j],$entropyP[$i][$j],$StdvEntropyP[$i][$j]), "\n";
	for($k=0;$k<=$max;$k++) # Print extincted species as a function of degree
	{
	    if($ExtinctAll != 0){$DegVec[$k]/=$ExtinctAll;}
	    print HANDLE4 join(" ",$selectValues[$i],$Delta[$j],$k,$DegVec[$k],$ExtinctAll), "\n";
	    $DegVec[$k] = 0; # Initialize for the next set of realizations
	}

    }
}
for($j=0;$j<=$length2-1;$j++)
{
    for($i=1;$i<=$length1-1;$i++)
    {
	#print join(" ",$i,$j,$Nrnd,$Np,$Na,$BiomassP[$i][$j],$StdvBioP[$i][$j],$BiomassA[$i][$j],$StdvBioA[$i][$j],'DEBUG2'),"\n";
	$dBioP=($BiomassP[$i][$j]-$BiomassP[0][$j]); $dBioA=($BiomassA[$i][$j]-$BiomassA[0][$j]);
	$Err_dBioP=sqrt(($StdvBioP[$i][$j]+$StdvBioP[0][$j])/$Nrnd); $Err_dBioA=sqrt(($StdvBioA[$i][$j]+$StdvBioA[0][$j])/$Nrnd);

	$dSp=-($ExtinctP[$i][$j]-$ExtinctP[0][$j]); $dSa=-($ExtinctA[$i][$j]-$ExtinctA[0][$j]); $dS=$dSp+$dSa;
	$Err_dSp=sqrt(($StdvExtP[$i][$j]+$StdvExtP[0][$j])/$Nrnd); $Err_dSa=sqrt(($StdvExtA[$i][$j]+$StdvExtA[0][$j])/$Nrnd); $Err_dS=$Err_dSp+$Err_dSa;

	$ErrSurvP0=sqrt($StdvSurvP[$i][$j]/$Nrnd);$ErrSurvA0=sqrt($StdvSurvA[$i][$j]/$Nrnd);
	$Survive0=$SurviveP[0][$j]+$SurviveA[0][$j]; $ErrSurv0=$ErrSurvP0+$ErrSurvA0;

	$dSpRel=$dSp/$SurviveP[0][$j]; $dSaRel=$dSa/$SurviveA[0][$j]; $dSrel=$dS/$Survive0;

	if($dSp == 0){$Err_dSpRel=0;}else{$Err_dSpRel=($dSpRel)*($Err_dSp/$dSp+$ErrSurvP0/$SurviveP[0][$j]);} # If there is no difference the error will diverge 
	if($dSa == 0){$Err_dSaRel=0;}else{$Err_dSaRel=($dSaRel)*($Err_dSa/$dSa+$ErrSurvA0/$SurviveA[0][$j]);}
	if($dS == 0){$Err_dSrel=0;}else{$Err_dSrel=($dSrel)*($Err_dS/$dS+$ErrSurv0/$Survive0);}

	$dNestOutP=($NestOutP[$i][$j]-$NestednessP);
	$Err_dNestOutP=sqrt($StdvNestOutP[$i][$j]/$Nrnd);
	$dNestOutA=($NestOutA[$i][$j]-$NestednessA);
	$Err_dNestOutA=sqrt($StdvNestOutA[$i][$j]/$Nrnd);
	$dNestOutT=($NestOutT[$i][$j]-$NestednessT);
	$Err_dNestOutT=sqrt($StdvNestOutT[$i][$j]/$Nrnd);
	$dConnectOut=($ConnectOut[$i][$j]-$ConnectOut[0][$j]);
	$Err_dConnectOut=sqrt(($StdvConnectOut[$i][$j]+$StdvConnectOut[0][$j])/$Nrnd);

	$dSinglesOut=($SinglesOut[$i][$j]-$SinglesOut[0][$j]);
	$Err_dSinglesOut=sqrt(($StdvSinglesOut[$i][$j]+$StdvSinglesOut[0][$j])/$Nrnd);
	$dSinglesExtinct=($SinglesExtinct[$i][$j]-$SinglesExtinct[0][$j]);
	$Err_dSinglesExtinct=sqrt(($StdvSinglesExtinct[$i][$j]+$StdvSinglesExtinct[0][$j])/$Nrnd);

	$dExcludeP=($ExcludeP[$i][$j]-$ExcludeP[0][$j])/$Sp; $dExcludeA=($ExcludeA[$i][$j]-$ExcludeA[0][$j])/$Sa;
	$Err_dExcludeP=sqrt(($StdvExcludeP[$i][$j]+$StdvExcludeP[0][$j])/($Nrnd*$Sp)); $Err_dExcludeA=sqrt(($StdvExcludeA[$i][$j]+$StdvExcludeA[0][$j])/($Nrnd*$Sa));

	
	print HANDLE3 join(" ",$selectValues[$i],$Delta[$j],$dBioP,$Err_dBioP,$dBioA,$Err_dBioA,$dSp,$Err_dSp,$dSa,$Err_dSa,$dSpRel,$Err_dSpRel,$dSaRel,$Err_dSaRel,$dSrel,$Err_dSrel,$dNestOutP,$Err_dNestOutP,$dNestOutA,$Err_dNestOutA,$dNestOutT,$Err_dNestOutT,$dConnectOut,$Err_dConnectOut,$dSinglesOut,$Err_dSinglesOut,$dSinglesExtinct,$Err_dSinglesExtinct,$dExcludeP,$Err_dExcludeP,$dExcludeA,$Err_dExcludeA), "\n";
    }
}	
`rm -rf  dir*`;
close(HANDLE2); 

print '  ',"\n";
print '****************************',"\n";
print '** Program finished',"\n";
print '** Check your results  ',"\n";
print '** Bye! ',"\n";
print '****************************',"\n";
print '  ',"\n";
print '  ',"\n";
