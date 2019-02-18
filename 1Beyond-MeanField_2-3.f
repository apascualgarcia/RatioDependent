*********************************************************************
*                       Beyond-MeanField.f                          *
*********************************************************************
*
*     Our aim with this program is to explore the parameters regime
*     where we find coexistence for a plant-pollinator network of
*     populations. Following our previous work at Bastolla et al. Nature (2009),
*     we want to know whether the soft mean field assumption can
*     be further relaxed for the mutualistic parameters.
*     We already know that there is a relation between the persistence
*     of the species and the width of the variance of the productivity
*     vectors, as it must be as parallel as possible with respect
*     to the main eigenvector of the normalized effective competition
*     matrix.
*     Finding the correct parameter values by chance out of the
*     mean field approach seems to be rather hard, and we
*     propose here a way to circunvent the numerous choice of
*     parameters. 
*     Given a network with Sa (Sp) species of animals (plants) we are
*     going to fix the range of the distribution for most of the parameters,
*     selecting the abundances in such a way that we
*     are able to get positive bare productivities for plants and
*     negative for animals. In this way we expect to deal
*     with values of the bare productivities that, although they may
*     not be tight, they are at least approximately parallel to
*     the main eigenvector.
*     *****
*                                               Madrid, July 2012
*                                   Alberto Pascual-García (CBMSO) 
*                                            apascual@cbm.uam.es     
*     *****
*     ADAPTIVE VERSION: In this version we include integration routines aimed to 
*     integrate the ODEs with a variable step. We integrate with a Bulirsch-Stoer Method
*     including the Richardson Extrapolation*      
*     *****
*     COMPILATION:
*     Compiled in x86-64bits under Linux with gnu Fortran 77 with the following options:
*     g77 -O2 -ffixed-line-length-0 -march=nocona -fcase-upper -o Name.exe *.f (likely obsolet)
*     gfortran -O2 -ffixed-line-length-0 -march=nocona -o Name.exe *.f 
*     *****
*
*     TREE from root (Main):
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
*
*     ..... Normal_Exit
*
*******************************************************************

*    **********************************************
*     MAIN PROGRAM
*    **********************************************                                                  

      IMPLICIT NONE
      INTEGER*4 today(3),now(3)
      INTEGER unit(30),Nfiles,NfilesIn
      INTEGER sizeA,sizeP,sizeT
      CHARACTER*40 path(30)
      REAL*8 seed,fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     ...COMMONS
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 rhoP,widthRhoP,rhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp
      COMMON/Species/ Sa,Sp

      PRINT *, ' ' 
      PRINT *,'*******************************************'
      PRINT *,'* Mutualism Simulations Beyond Mean Field *'
      PRINT *,'*******************************************'
      PRINT *, ' ' 

      call Read_Parameters(Nfiles,NfilesIn,unit,today,now,path,
     &seed,fixedPoint,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa)
      sizeA=Sa                  ! Attempt to define dynamically the matrices
      sizeP=Sp
      sizeT=Sa+Sp
      call Pseudo_Main(path,Nfiles,NfilesIn,unit,today,now,seed,
     &fixedPoint,sizeA,sizeP,sizeT,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa)
      call Normal_Exit()
      END

*    **********************************************
*     END MAIN PROGRAM
*    **********************************************                                                  
*_______________________________________________________________________
*
*    **********************************************
*     SUBROUTINE Pseudo_Main
*    **********************************************                                                  


      SUBROUTINE Pseudo_Main(path,Nfiles,NfilesIn,unit,today,now,
     &seed,fixedPoint,sizeA,sizeP,sizeT,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa)
      IMPLICIT NONE
      CHARACTER*40 path(30)
      INTEGER sizeA,sizeP,sizeT
      INTEGER unit(30),Nfiles,NfilesIn
      INTEGER*4 today(3), now(3)
      REAL*8 seed,fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     This subroutine is just a trick to define dynamically
*     the matrix sizes, as it is not possible to do it in the main program
*     within the static paradigm of Fortran 77
      REAL AdjA(sizeA,sizeP),AdjP(sizeP,sizeA)
      REAL degreeP(sizeP),degreeA(sizeA)
      INTEGER excludedA,excludedP,NegAlphaP,NegAlphaA
      REAL Connect,single
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP),seedout
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA),NestA,NestP,NestT
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 Dissipation,DissRate,AvDissipation,AvDissRate
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 AvA,AvP,VarA,VarP

      call Read_Generate_Inputs(path,Nfiles,NfilesIn,unit,today,now,seed,fixedPoint,
     &sizeA,sizeP,sizeT,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa,BetaA,BetaP,GammaA,GammaP,Na,Np,u,AlphaA,AlphaP,
     &degreeP,degreeA,NestA,NestP,NestT,AdjA,AdjP,Connect,single,
     &excludedP,excludedA,NegAlphaP,NegAlphaA,hhA,hhP,seedout,AvA,AvP,VarA,VarP)
      call Integration(Nfiles,unit,sizeA,sizeP,sizeT,BetaA,BetaP,
     &GammaA,GammaP,Na,Np,u,AlphaA,AlphaP,hhA,hhP,Dissipation,AvDissipation,
     &DissRate,AvDissRate)
      call Output_Globals(unit,sizeA,sizeP,sizeT,BetaA,BetaP,GammaA,GammaP,Na,Np,
     &u,AlphaA,AlphaP,degreeP,degreeA,NestA,NestP,NestT,AdjA,AdjP,Connect,single,
     &Dissipation,AvDissipation,DissRate,AvDissRate,excludedP,excludedA,
     &NegAlphaP,NegAlphaA,seedout,AvA,AvP,VarA,VarP)

      RETURN
      END

*    **********************************************
*     Subroutine Read_Parameters
*    **********************************************                                                  

      SUBROUTINE Read_Parameters(Nfiles,NfilesIn,unit,today,now,path,
     &seed,fixedPoint,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa)
      IMPLICIT NONE
*     Read parameters and paths from .in file. We do not use "include"
*     statements, although it would be easier to handle, because
*     we want to run the program massively avoiding the
*     recompilation before each realization.
      INTEGER i,Z,ZZ,tmp
      INTEGER iADJUSTL, iTRIM
      INTEGER Nfiles,NfilesIn,NfilesOut,unit(30),unitTmp
      CHARACTER*40 path0,path(30)
      CHARACTER*200 pathSummary
      INTEGER*4 today(3), now(3)
      REAL*8 seed,fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     ...Commons (tragedy of)
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 rhoP,widthRhoP,rhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp
      COMMON/Species/ Sa,Sp
       
      Nfiles=26
      NfilesIn=12
      NfilesOut=Nfiles-NfilesIn
      path0='Beyond-MeanField.in'
c      PRINT *, '  * Reading files from: ',Path0
      OPEN(UNIT=69,STATUS='OLD',ERR=769,FILE=path0)
      READ(69,*) readNp         ! One if you read Np from file, 0 if will be generated internally
      READ(69,*) midNp          ! Distribution of densities centered at midNp.      
      READ(69,*) widthNp
      READ(69,*) readNa         ! One if you read Np from file, 0 if will be generated internally
      READ(69,*) midNa          ! Distribution of densities centered at this value
      READ(69,*) widthNa        ! with this width
      READ(69,*) readAlphaP         ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midAlphaP      ! Distribution of densities centered ..(zero if reading from file)
      READ(69,*) widthAlphaP
      READ(69,*) readAlphaA     ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midAlphaA      ! Distribution of densities centered ..(zero if reading from file)
      READ(69,*) widthAlphaA
      READ(69,*) readBetaP      ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midBetaP         ! Intraspecific competition for plants
      READ(69,*) widthBetaP     ! Width of the distribution
      READ(69,*) readBetaA      ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midBetaA         ! Intraspecific competition for animals
      READ(69,*) widthBetaA     ! Width of the distribution
      READ(69,*) midRhoP           ! Constant for interspecific competition (plants)
      READ(69,*) widthRhoP      ! Width of the distribution
      READ(69,*) midRhoA           ! Constant for interspecific competition (animals)
      READ(69,*) widthRhoA      ! Width of the distribution
      READ(69,*) readGammaP     ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGammaP        ! Coupling interactions with plants in rows and animals in columns
      READ(69,*) widthGammaP
      READ(69,*) readGammaA     ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGammaA        ! Coupling interactions with plants in rows and animals in columns
      READ(69,*) widthGammaA
      READ(69,*) readHp         ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midHp             ! handling time for plants with respect to animals abundances
      READ(69,*) widthHp        ! Width of the distribution
      READ(69,*) readHa         ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midHa             ! handling time for animals
      READ(69,*) widthHa        ! Width of the distribution
      READ(69,*) readGp         ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGp             ! handling time for plants with respect to plants abundances
      READ(69,*) widthGp        ! Width of the distribution
      READ(69,*) readGa         ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGa             ! handling time for animals
      READ(69,*) widthGa        ! Width of the distribution
      READ(69,*) Sp             ! Number of species of plants
      READ(69,*) Sa             ! of animals
      READ(69,*) Delta          ! Amplitude of the growth rates fluctuations
      READ(69,*) fixedPoint     ! Parameter controlling if alphas should be estimated at a fixed point, will override any value given to alpha
      path(1)='plantsIn.dat'        
      path(2)='animalsIn.dat'
      path(3)='betaInP.dat'
      path(4)='betaInA.dat'
      path(5)='gammaInP.dat'
      path(6)='gammaInA.dat'
      path(7)='alphaInP.dat'
      path(8)='alphaInA.dat'
      path(9)='hhandlingInP.dat' ! Four new input  files
      path(10)='hhandlingInA.dat'
      path(11)='ghandlingInP.dat'
      path(12)='ghandlingInA.dat'
      path(13)='plantsOut-Final.out'        
      path(14)='animalsOut-Final.out'
      path(15)='betaOutP.out'
      path(16)='betaOutA.out'
      path(17)='gammaOutP.out'
      path(18)='gammaOutA.out'
      path(19)='alphaOutP.out'
      path(20)='alphaOutA.out'
      path(21)='hhandlingOutP.dat' ! Four new output  files, I leave the output files for paths the last ones
      path(22)='hhandlingOutA.dat'
      path(23)='ghandlingOutP.dat'
      path(24)='ghandlingOutA.dat'
      path(25)='plantsOut-Paths.out'        
      path(26)='animalsOut-Paths.out'
      DO i=1,NfilesIn
         unit(i)=20+i
      ENDDO
      DO i=1,NfilesOut !2 ! NfilesOut ! We just need to open two files
         tmp=NfilesIn+i
         unit(tmp)=50+i
         OPEN(UNIT=unit(tmp),STATUS='NEW',ERR=771,FILE=path(tmp))
c         unit(Nfiles+1-i)=40+i ! I use this line and the following to print the integration paths (I print four files in this way), comment otherwise
c         OPEN(UNIT=unit(Nfiles+1-i),STATUS='NEW',ERR=771,FILE=path(Nfiles+1-i))
c         PRINT *, '  * Unit ',unit(NfilesIn+i),' corresponds to ',path(NfilesIn+i)
      ENDDO
      READ(69,'(a)') pathSummary
c      PRINT *, '  * The path for the output summary is: ',pathSummary
      READ(69,*) seed ! read a random number seed from the perl shuttle
      unitTmp=81
      unit(Nfiles+2)=unitTmp
      Z=iTRIM(pathSummary)
      ZZ=iADJUSTL(pathSummary)
      OPEN(UNIT=unitTmp,STATUS='NEW',ERR=773,FILE=pathSummary(ZZ:Z))
      call Header(pathSummary,unitTmp,today,now,
     &readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa) ! Write the header for the output files      
      RETURN
 769  CONTINUE
      PRINT *,'  * * Err 769: Problems opening file .in: ',path0
      call Warning()
 771  CONTINUE
      PRINT *,'  * *  Err 771: Problems opening unit: ',unit(NfilesIn+i),' >> ',path(i)
      call Warning()
 773  CONTINUE
      PRINT *,'  * *  Err 773: Problems opening unit: ',unitTmp,' >> ',pathSummary
      call Warning()
      RETURN
      END

*    **********************************************
*     Subroutine Header
*    **********************************************                                                  

      SUBROUTINE Header(pathOut,unit,today,now,
     &readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa)
      IMPLICIT NONE
      INTEGER unit
      CHARACTER*200 pathOut
*     Writes the output file header.
      INTEGER*4 today(3), now(3)
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     ...Commons
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 rhoP,widthRhoP,rhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp
      COMMON/Species/ Sa,Sp
 
      call When(Today,Now)      ! Compute the date for the file header
      WRITE(unit,*) '# Beyond-MeanField.f'
      WRITE(unit,111) Today,Now
      WRITE(unit,112) '# Summary output: ',pathOut

      IF(readNp .eq. 1)THEN
         WRITE(unit,*) 'Np was read from file '
      ELSE
         WRITE(unit,*) 'midNp ',midNp ! Distribution of densities centered ..
         WRITE(unit,*) 'widthNp ',widthNp
      ENDIF
      IF(readNa .eq. 1)THEN
         WRITE(unit,*) 'Na was read from file '
      ELSE
         WRITE(unit,*) 'midNa ',midNa ! Distribution of densities centered at this value
         WRITE(unit,*) 'widthNa ',widthNa ! with this width 
      ENDIF
      IF(readAlphaP .eq. 1)THEN
         WRITE(unit,*) 'AlphaP was read from file '
      ELSE
         WRITE(unit,*) 'midAlphaP ', midAlphaP ! Distribution of densities centered ..
         WRITE(unit,*) 'widthAlphaP ',widthAlphaP
      ENDIF
      IF(readAlphaA .eq. 1)THEN
         WRITE(unit,*) 'AlphaA was read from file '
      ELSE
         WRITE(unit,*) 'midAlphaA ', midAlphaA ! Distribution of densities centered ..
         WRITE(unit,*) 'widthAlphaA ',widthAlpha
      ENDIF
      IF(readBetaP .eq. 1)THEN
         WRITE(unit,*) 'BetaP and rhoP were read from file '
      ELSE
         WRITE(unit,*) 'midBetaP ', midBetaP ! Paramater multiplying a distribution centered around one
         WRITE(unit,*) 'widthBetaP ', widthBetaP ! Width of the distribution
         WRITE(unit,*) 'midRhoP ',midRhoP ! Constant for interspecific competition for plants
         WRITE(unit,*) 'widthRhoP ',widthRhoP ! Constant for interspecific competition for plants
      ENDIF
      IF(readBetaA .eq. 1)THEN
         WRITE(unit,*) 'BetaA and rhoA were read from file '
      ELSE
         WRITE(unit,*) 'midBetaA ', midBetaA ! Parameter multiplying a distribution centered around one
         WRITE(unit,*) 'widthBetaA ', widthBetaA ! Width of the distribution
         WRITE(unit,*) 'midRhoA ',midRhoA ! Constant for interspecific competition for plants
         WRITE(unit,*) 'widthRhoA ',widthRhoA ! Constant for interspecific competition for plants
      ENDIF
      IF(readGammaP .eq. 1)THEN
         WRITE(unit,*) 'GammaP was read from file '
      ELSE
         WRITE(unit,*) 'midGammaP ',midGammaP ! Paramater multiplying a distribution centered around one
         WRITE(unit,*) 'widthGammaP ',widthGammaP
      ENDIF
      IF(readGammaA .eq. 1)THEN
         WRITE(unit,*) 'GammaA was read from file '
      ELSE
         WRITE(unit,*) 'midGammaA ',midGammaA ! Paramater multiplying a distribution centered around one
         WRITE(unit,*) 'widthGammaA ',widthGammaA    
      ENDIF
      IF(readHp .eq. 1)THEN
         WRITE(unit,*) 'hP was read from file '
      ELSE
         WRITE(unit,*) 'midHp ',midHp ! handling time for plants
         WRITE(unit,*) 'widthHp ',widthHp ! handling time for plant
      ENDIF
      IF(readHa .eq. 1)THEN
         WRITE(unit,*) 'hA was read from file '
      ELSE
         WRITE(unit,*) 'midHa ',midHa ! handling time for plants
         WRITE(unit,*) 'widthHa ',widthHa ! handling time for plants 
      ENDIF
      IF(readGp .eq. 1)THEN
         WRITE(unit,*) 'gP was read from file '
      ELSE
         WRITE(unit,*) 'midGp ',midGp ! handling time for plants
         WRITE(unit,*) 'widthGp ',widthGp ! handling time for plant
      ENDIF
      IF(readGa .eq. 1)THEN
         WRITE(unit,*) 'gA was read from file '
      ELSE
         WRITE(unit,*) 'midGa ',midGa ! handling time for plants
         WRITE(unit,*) 'widthGa ',widthGa ! handling time for plants 
      ENDIF
      WRITE(unit,*) 'Sp ',Sp    ! of plants
      WRITE(unit,*) 'Sa ',Sa    ! Number of species of animals
      WRITE(unit,*) 'Delta  ',Delta ! Amplitude of the growth rates fluctuations
      
 111  format ( ' # Running at: Date ', i2.2, '/', i2.2, '/', i4.4, 
     &     '; Time ',i2.2, ':', i2.2, ':', i2.2 )
 112  format (A20,' ',A120)
      Return
      END

*    **********************************************
*     Function When
*    **********************************************                                                  

      SUBROUTINE When(today,now)
*     It brings the date and time the program is being running
      integer*4 today(3), now(3)

      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second

      Return
      end

*     *********************************************
*     SUBROUTINE Writer
*     *********************************************                                                  
 
      SUBROUTINE Output_Globals(unit,sizeA,sizeP,sizeT,BetaA,BetaP,GammaA,GammaP,Na,Np,
     &u,AlphaA,AlphaP,degreeP,degreeA,NestA,NestP,NestT,AdjA,AdjP,Connect,single,
     &Dissipation,AvDissipation,DissRate,AvDissRate,excludedP,excludedA,
     &NegAlphaP,NegAlphaA,seedout,AvA,AvP,VarA,VarP)
      IMPLICIT NONE
      REAL Thr
      PARAMETER(Thr=1e-8)
      INTEGER unit(30),i,j,k
      INTEGER sizeA,sizeP,sizeT
      REAL degreeP(sizeP),degreeA(sizeA)
      REAL AdjA(sizeA,sizeP),AdjP(sizeP,sizeA)
      REAL single,SingleExtinct
      INTEGER excludedP,excludedA,NegAlphaP,NegAlphaA
      REAL Connect
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 NestA,NestP,NestT
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 AvA,AvP,VarA,VarP
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 Dissipation,DissRate,AvDissipation,AvDissRate,seedout
*     --- Writes a summary computing total biomass, dissipation and extinctions.
      LOGICAL alive(sizeT),key_i,key_j
      REAL degreePout(sizeP),degreeAout(sizeA)
      REAL*8 BiomassP,BiomassA,ExtinctP,ExtinctA
      REAL*8 NestAout,NestPout,NestTout
      REAL singleOut
      REAL ConnectOut
      INTEGER excludedPout,excludedAout
      REAL*8 f(sizeT),Frac
      REAL*8 ShannonA,ShannonP
      REAL*8 AvDegreeP,AvDegreeA,VarDegreeP,VarDegreeA
*     ...Commons
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 rhoP,widthRhoP,rhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp
      COMMON/Species/ Sa,Sp

      BiomassP=0.0d0
      ExtinctP=0.0d0
      BiomassA=0.0d0
      ExtinctA=0.0d0
      singleOut=0
      singleExtinct=0
      DO i=1,sizeP              ! Test for surviving species
         alive(i)=.true.
         IF(u(i).lt.Thr)THEN
            u(i)=0.0d0
            ExtinctP=ExtinctP+1
            IF(degreeP(i).eq.1) SingleExtinct=SingleExtinct+1
            alive(i)=.false.
         ELSE
            BiomassP=BiomassP+u(i)
         ENDIF
      ENDDO
      DO i=1,sizeA
         alive(i+sizeP)=.true.
         IF(u(i+sizeP).lt.Thr)THEN 
            u(i+sizeP)=0.0d0
            ExtinctA=ExtinctA+1
            IF(degreeA(i).eq.1) singleExtinct=singleExtinct+1
            alive(i+sizeP)=.false.
         ELSE
            BiomassA=BiomassA+u(sizeP+i)
         ENDIF
      ENDDO
c$$$      DO i=1,sizeP
c$$$         PRINT *, i,alive(i),' debug alive inside i'
c$$$      ENDDO
      call Topology(sizeA,sizeP,sizeT,AdjP,AdjA,alive,degreePout,
     &     degreeAout,NestPout,NestAout,NestTout,ConnectOut,
     &     singleOut,excludedPout,excludedAout)
      AvDegreeP=0
      VarDegreeP=0
      ShannonP=0
      ShannonA=0
      WRITE(unit(9),*) '#Final Plant Values'
      WRITE(unit(9),*) '#SpecieId, BiomassIn, BiomassOut, Alpha, DegreeIn, 
     & DegreeOut, Alive?'
      DO i=1,sizeP
         AvDegreeP=AvDegreeP+degreeP(i)
         VarDegreeP=VarDegreeP+degreeP(i)**2
         IF(alive(i).eqv..false.)THEN
            !degreePout(i)=0
         ELSE
            Frac=u(i)/BiomassP
            ShannonP=ShannonP+Frac*LOG(Frac)
         ENDIF
         WRITE(unit(9),*) i,Np(i),u(i),AlphaP(i),
     &degreeP(i),degreePout(i),alive(i)
      ENDDO
      ShannonP=EXP(-ShannonP)
      AvDegreeA=0
      VarDegreeA=0
      WRITE(unit(10),*) '#Final Animals Values'
      WRITE(unit(10),*) '#SpecieId, BiomassIn, BiomassOut, Alpha, DegreeIn, 
     & DegreeOut, Alive?'
      DO i=1,sizeA
         AvDegreeA=AvDegreeA+degreeA(i)
         VarDegreeA=VarDegreeA+degreeA(i)**2
         IF(alive(sizeP+i).eqv..false.)THEN
            !degreeAout(i)=0
         ELSE
            Frac=u(i+sizeP)/BiomassA
            ShannonA=ShannonA+Frac*LOG(Frac)
         ENDIF
         WRITE(unit(10),*) i,Na(i),u(sizeP+i),AlphaA(i),
     &degreeA(i),degreeAout(i),alive(sizeP+i)
      ENDDO
      ShannonA=EXP(-ShannonA)
      AvDegreeP=AvDegreeP/float(sizeP)
      VarDegreeP=VarDegreeP/float(sizeP)-AvDegreeP**2
      AvDegreeA=AvDegreeA/float(sizeA)
      VarDegreeA=VarDegreeA/float(sizeA)-AvDegreeA**2
      IF(single.gt.0)THEN
         SingleExtinct=(single-SingleExtinct)/single
      ENDIF

      WRITE(61,*) 'BiomassP ', BiomassP
      WRITE(61,*) 'BiomassA ',BiomassA
      WRITE(61,*) 'ExtinctP', ExtinctP
      WRITE(61,*) 'ExtinctA', ExtinctA
      WRITE(61,*) 'SurvivingP',sizeP-ExtinctP
      WRITE(61,*) 'SurvivingA',sizeA-ExtinctA
      WRITE(61,*) 'NestednessPin',NestP
      WRITE(61,*) 'NestednessAin',NestA
      WRITE(61,*) 'NestednessTin',NestT
      WRITE(61,*) 'NestednessPout',NestPout
      WRITE(61,*) 'NestednessAout',NestAout
      WRITE(61,*) 'NestednessTout',NestTout
      WRITE(61,*) 'ConnectanceIn', Connect
      WRITE(61,*) 'ConnectanceOut', ConnectOut
      WRITE(61,*) 'SinglesIn', single
      WRITE(61,*) 'SinglesOut', singleOut
      WRITE(61,*) 'SingleExtinct', SingleExtinct
      WRITE(61,*) 'excludedPin', excludedP
      WRITE(61,*) 'excludedPout', excludedPout
      WRITE(61,*) 'excludedAin', excludedA
      WRITE(61,*) 'excludedAout', excludedAout
      WRITE(61,*) 'NegAlphaP', NegAlphaP
      WRITE(61,*) 'NegAlphaA', NegAlphaA
      WRITE(61,*) 'Dissipation',Dissipation
      WRITE(61,*) 'AvDissipation',AvDissipation
      WRITE(61,*) 'DissRate',DissRate
      WRITE(61,*) 'AvDissRate',AvDissRate
      WRITE(61,*) 'AvAlphaP',AvP
      WRITE(61,*) 'VarAlphaP',VarP
      WRITE(61,*) 'AvAlphaA',AvA
      WRITE(61,*) 'VarAlphaA',VarA
      WRITE(61,*) 'AvDegreeInP',AvDegreeP
      WRITE(61,*) 'VarDegreeInP',VarDegreeP
      WRITE(61,*) 'AvDegreeInA',AvDegreeA
      WRITE(61,*) 'VarDegreeInA',VarDegreeA
      WRITE(61,*) 'BiomassRate0',midNp
      WRITE(61,*) 'MaxhP: ',hP
      WRITE(61,*) 'MaxhA: ',hA
      WRITE(61,*) 'MaxRhoP: ',rhoP
      WRITE(61,*) 'MaxRhoA: ',rhoA
      WRITE(61,*) 'Seed: ',seedout
      WRITE(61,*) 'EntropyP ',ShannonP
      WRITE(61,*) 'EntropyA ',ShannonA

c$$$
c$$$      PRINT *, ''
c$$$      PRINT *, '*******************************'
c$$$      PRINT *, '* RESULTS:'
c$$$      PRINT *, '* BiomassP ', BiomassP
c$$$      PRINT *, '* BiomassA ',BiomassA
c$$$      PRINT *, '* ExtinctP', ExtinctP
c$$$      PRINT *, '* ExtinctA', ExtinctA
c$$$      PRINT *, '* SurvivingP',sizeP-ExtinctP
c$$$      PRINT *, '* SurvivingA',sizeA-ExtinctA
c$$$      PRINT *, '* NestednessPin',NestP
c$$$      PRINT *, '* NestednessAin',NestA
c$$$      PRINT *, '* NestednessTin',NestT
c$$$      PRINT *, '* NestednessPout',NestPout
c$$$      PRINT *, '* NestednessAout',NestAout
c$$$      PRINT *, '* NestednessTout',NestTout
c$$$      PRINT *, '* ConnectanceIn', Connect
c$$$      PRINT *, '* ConnectanceOut', ConnectOut
c$$$      PRINT *, '* SinglesIn', single
c$$$      PRINT *, '* SinglesOut', singleOut
c$$$      PRINT *, '* SingleExtinct', SingleExtinct
c$$$      PRINT *, '* excludedPin', excludedP
c$$$      PRINT *, '* excludedPout', excludedPout
c$$$      PRINT *, '* excludedAin', excludedA
c$$$      PRINT *, '* excludedAout', excludedAout
c$$$      PRINT *, '* NegAlphaP', NegAlphaP
c$$$      PRINT *, '* NegAlphaA', NegAlphaA
c$$$      PRINT *, '* Dissipation',Dissipation
c$$$      PRINT *, '* AvDissipation',AvDissipation
c$$$      PRINT *, '* DissRate',DissRate
c$$$      PRINT *, '* AvDissRate',AvDissRate
c$$$      PRINT *, '*******************************'
c$$$      PRINT *, ''
      RETURN
      END


*     *********************************************
*     Function iTRIM
*     *********************************************                                                  

      INTEGER FUNCTION iTRIM(String)
*     Returns the index of the last non blank character
      IMPLICIT NONE
      CHARACTER*(*) String
      INTEGER Length,Z

      Length=LEN(String)       ! Get the length,                      
      Z=0
      DO WHILE(Z .ne. Length)   ! Look for leading blanks
         Z=Z+1
         IF(String(Z:Z).ne.' ') EXIT
      ENDDO    
      DO WHILE(Z .ne. Length)   ! Look for trailing blanks (the relevant final position)
         Z=Z+1
         IF((String(Z:Z).eq.' ').or.(String(Z:Z).eq."\t")) EXIT
      ENDDO
      IF(Z.ge.LEN(String))THEN
         Z=LEN(String)
      ELSE
         Z=Z-1
      ENDIF

      iTRIM=Z

      RETURN
      END

*     *********************************************
*     Function iADJUSTL
*     *********************************************                                                  

      INTEGER FUNCTION iADJUSTL(String)
*     Returns the index of the first non blank character
      IMPLICIT NONE
      CHARACTER*(*) String
      INTEGER Length,Z

      Length=LEN(String)       ! Get the length,                   
      Z=0
      DO WHILE(Z .ne. Length)   ! Look for leading blanks
         Z=Z+1
         IF(String(Z:Z).ne.' ') EXIT
      ENDDO    
      iADJUSTL=Z

      RETURN
      END

*     *********************************************
*     Subroutine Warning
*     *********************************************                                                  

      SUBROUTINE Warning()
*     Send a warning exit and stop
      PRINT *, ' *********************'
      PRINT *, ' * WARNING !'
      PRINT *, ' *********************'
      PRINT *, ' * Abrupt program termination'
      PRINT *, ' *********************'
      PRINT *, ' *********************'
      STOP
      END

*     *********************************************
*     SUBROUTINE Normal Exit
*     *********************************************                                                  

      SUBROUTINE Normal_Exit()

      PRINT *, ' ' 
      PRINT *, ' >> FINISHING NORMALLY!'
      PRINT *, '  ** Check the results '
      PRINT *, '  ** See You! ' 
      PRINT *, ' ' 
      PRINT *, '<<<<<<<<<<<<..>>>>>>>>>>>>>>'
      PRINT *, '>>>>>>>><<<<..>>><<<<<<<<<<<'
      PRINT *, '>>>>>>>>>><<...>><<<<<<<<<<<'
      PRINT *, '<<<<<<<<<>>>><<<<>>>>>>>>>>>' 
      PRINT *, '<<<<<<<<<<______>>>>>>>>>>>>' 
      PRINT *, '<<<<<<<<<________>>>>>>>>>>>' 
      PRINT *, '____________________________'
      PRINT *, ' ' 

      RETURN
      END
