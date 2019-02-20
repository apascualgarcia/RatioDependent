*********************************************************************
*                       Beyond-MeanField.f                          *
*********************************************************************
*
*
*     *****
*                                               Zurich, February 2019
*                                   Alberto Pascual-García (CBMSO) 
*                                            alberto.pascual.garcia@gmail.com     
*     *****
*     ADAPTIVE VERSION: In this version we include integration routines aimed to 
*     integrate the ODEs with a variable step. We integrate with a Bulirsch-Stoer Method
*     including the Richardson Extrapolation*      
*     *****
*
*     TREE from root (Main):
*
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
*******************************************************************

*    **********************************************
*     MAIN PROGRAM
*    **********************************************                                                  

      IMPLICIT NONE
      INTEGER*4 today(3),now(3)
      INTEGER unit(30),Nfiles,NfilesIn
      INTEGER sizeA,sizeP,sizeT
      CHARACTER*40 path(30)
      REAL*8 seed
      INTEGER fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     ...COMMONS
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 midRhoP,widthRhoP,midRhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 f0P,f0A
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp,
     &f0P,f0A,Delta
      COMMON/Species/ Sa,Sp

      PRINT *, ' ' 
      PRINT *,'************************************************'
      PRINT *,'* RatioDependent Simulations Beyond Mean Field *'
      PRINT *,'************************************************'
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
      REAL*8 seed
      INTEGER fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     This subroutine is just a trick to define dynamically
*     the matrix sizes, as it is not possible to do it in the main program
*     within the static paradigm of Fortran 77
      REAL degreeP(sizeP),degreeA(sizeA)
      INTEGER excludedA,excludedP,NegAlphaP,NegAlphaA
      REAL Connect,single
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP),seedout
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA),NestA,NestP,NestT
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 ggA(sizeA),ggP(sizeP)
      REAL*8 AvA,AvP,VarA,VarP

      call Read_Generate_Inputs(path,Nfiles,NfilesIn,unit,today,now,seed,fixedPoint,
     &sizeA,sizeP,sizeT,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa,BetaA,BetaP,GammaA,GammaP,Na,Np,u,AlphaA,AlphaP,
     &degreeP,degreeA,NestA,NestP,NestT,Connect,single,
     &excludedP,excludedA,NegAlphaP,NegAlphaA,hhA,hhP,ggA,ggP,seedout,AvA,AvP,VarA,VarP)
      call Integration(Nfiles,unit,sizeA,sizeP,sizeT,BetaA,BetaP,
     &GammaA,GammaP,Na,Np,u,AlphaA,AlphaP,hhA,hhP,ggA,ggP)
      call Output_Globals(unit,Nfiles,sizeA,sizeP,sizeT,BetaA,BetaP,GammaA,GammaP,Na,Np,
     &u,AlphaA,AlphaP,degreeP,degreeA,NestA,NestP,NestT,Connect,single,
     &excludedP,excludedA,NegAlphaP,NegAlphaA,seedout,AvA,AvP,VarA,VarP)

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
      CHARACTER*200 pathSummary,junk
      INTEGER*4 today(3), now(3)
      REAL*8 seed
      INTEGER fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     ...Commons (tragedy of)
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 midRhoP,widthRhoP,midRhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 f0P,f0A
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp,
     &f0P,f0A,Delta
      COMMON/Species/ Sa,Sp
       
      Nfiles=28
      NfilesIn=12
      NfilesOut=Nfiles-NfilesIn
      path0='Beyond-MeanField.in'
c      PRINT *, '  * Reading files from: ',Path0
      OPEN(UNIT=69,STATUS='OLD',ERR=769,FILE=path0)
      
      READ(69,111) readNp,junk  ! One if you read Np from file, 0 if will be generated internally
      READ(69,*) midNp ,junk  ! Distribution of densities centered at midNp.
      READ(69,*) widthNp,junk
      READ(69,111) readNa ,junk       ! One if you read Np from file, 0 if will be generated internally
      READ(69,*) midNa,junk       ! Distribution of densities centered at this value
      READ(69,*) widthNa,junk     ! with this width
      READ(69,111) readAlphaP,junk    ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midAlphaP,junk     ! Distribution of densities centered ..(zero if reading from file)
      READ(69,*) widthAlphaP ,junk
      READ(69,111) readAlphaA ,junk     ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midAlphaA,junk     ! Distribution of densities centered ..(zero if reading from file)
      READ(69,*) widthAlphaA ,junk
      READ(69,111) readBetaP,junk     ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midBetaP ,junk    ! Intraspecific competition for plants
      READ(69,*) widthBetaP ,junk     ! Width of the distribution
      READ(69,*) midRhoP,junk       ! Constant for interspecific competition (plants)
      READ(69,*) widthRhoP,junk     ! Width of the distribution
      READ(69,111) readBetaA ,junk      ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midBetaA,junk      ! Intraspecific competition for animals
      READ(69,*) widthBetaA,junk    ! Width of the distribution
      READ(69,*) midRhoA ,junk      ! Constant for interspecific competition (animals)
      READ(69,*) widthRhoA ,junk    ! Width of the distribution
      READ(69,111) readGammaP,junk    ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGammaP,junk     ! Coupling interactions with plants in rows and animals in columns
      READ(69,*) widthGammaP ,junk
      READ(69,111) readGammaA ,junk     ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGammaA ,junk      ! Coupling interactions with plants in rows and animals in columns
      READ(69,*) widthGammaA ,junk
      READ(69,111) readHp, junk         ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midHp,junk       ! handling time for plants with respect to animals abundances
      READ(69,*) widthHp ,junk      ! Width of the distribution
      READ(69,111) readHa,junk        ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midHa ,junk        ! handling time for animals
      READ(69,*) widthHa,junk       ! Width of the distribution
      READ(69,111) readGp,junk        ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGp  ,junk       ! handling time for plants with respect to plants abundances
      READ(69,*) widthGp ,junk ! Width of the distribution
      READ(69,111) readGa ,junk          ! One if you read from file, 0 if will be generated internally            
      READ(69,*) midGa  ,junk          ! handling time for animals
      READ(69,*) widthGa ,junk         ! Width of the distribution
      READ(69,*) f0P  ,junk            ! Scale of the saturating term (plants)
      READ(69,*) f0A ,junk             ! Scale of the saturating term (animals)      
      READ(69,*) Sp ,junk              ! Number of species of plants
      READ(69,*) Sa ,junk              ! of animals
      READ(69,*) Delta  ,junk          ! Amplitude of the growth rates fluctuations
      READ(69,*) fixedPoint  ,junk     ! Parameter controlling if alphas should be estimated at a fixed point, will override any value given to alpha
      READ(69,'(a)') pathSummary
c      PRINT *, '  * The path for the output summary is: ',pathSummary
      READ(69,*) seed           ! read a random number seed from the perl shuttle
      
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
      path(13)='plantsOut.dat'        
      path(14)='animalsOut.dat'
      path(15)='betaOutP.dat'
      path(16)='betaOutA.dat'
      path(17)='gammaOutP.dat'
      path(18)='gammaOutA.dat'
      path(19)='alphaOutP.dat'
      path(20)='alphaOutA.dat'
      path(21)='hhandlingOutP.dat' ! Four new output  files, I leave the output files for paths the last ones
      path(22)='hhandlingOutA.dat'
      path(23)='ghandlingOutP.dat'
      path(24)='ghandlingOutA.dat'
      path(25)='plantsOut-Paths.dat'        
      path(26)='animalsOut-Paths.dat'
      path(27)='plantsOut-Final.dat'        
      path(28)='animalsOut-Final.dat'
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

      unitTmp=81
      unit(Nfiles+1)=unitTmp
      Z=iTRIM(pathSummary)
      ZZ=iADJUSTL(pathSummary)
      OPEN(UNIT=unitTmp,STATUS='NEW',ERR=773,FILE=pathSummary(ZZ:Z))
      call Header(pathSummary,unitTmp,today,now,
     &readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa,fixedPoint) ! Write the header for the output files
      
      RETURN
 111  FORMAT(I1,A200)
 112  FORMAT(F10.0,A200)
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
     &readHp,readHa,readGp,readGa,fixedPoint)
      IMPLICIT NONE
      INTEGER unit
      CHARACTER*200 pathOut
*     Writes the output file header.
      INTEGER*4 today(3), now(3)
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
      INTEGER fixedPoint
*     ...Commons
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 midRhoP,widthRhoP,midRhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 f0P,f0A
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp,
     &f0P,f0A,Delta
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
      IF(fixedPoint .eq.1)THEN
         WRITE(unit,*) 'AlphaP was estimated at a feasible fixed point '
      ELSE
         IF(readAlphaP .eq. 1)THEN
            WRITE(unit,*) 'AlphaP was read from file '
         ELSE
            WRITE(unit,*) 'midAlphaP ', midAlphaP ! Distribution of densities centered ..
            WRITE(unit,*) 'widthAlphaP ',widthAlphaP
         ENDIF
      ENDIF
      IF(fixedPoint .eq.1)THEN
         WRITE(unit,*) 'AlphaA was estimated at a feasible fixed point '
      ELSE
         IF(readAlphaA .eq. 1)THEN
            WRITE(unit,*) 'AlphaA was read from file '
         ELSE
            WRITE(unit,*) 'midAlphaA ', midAlphaA ! Distribution of densities centered ..
            WRITE(unit,*) 'widthAlphaA ',widthAlphaA
         ENDIF
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
      WRITE(unit,*) 'f0P ',f0P  ! of plants
      WRITE(unit,*) 'f0A ',f0A  ! of animals
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
 
      SUBROUTINE Output_Globals(unit,Nfiles,sizeA,sizeP,sizeT,BetaA,BetaP,GammaA,GammaP,Na,Np,
     &u,AlphaA,AlphaP,degreeP,degreeA,NestA,NestP,NestT,Connect,single,
     &excludedP,excludedA,NegAlphaP,NegAlphaA,seedout,AvA,AvP,VarA,VarP)
      IMPLICIT NONE
      REAL Thr
      PARAMETER(Thr=1e-8)
      INTEGER unit(30),i,j,k,tmp,Nfiles
      INTEGER sizeA,sizeP,sizeT
      REAL degreeP(sizeP),degreeA(sizeA)
      REAL single,SingleExtinct
      INTEGER excludedP,excludedA,NegAlphaP,NegAlphaA
      REAL Connect
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 NestA,NestP,NestT
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 AvA,AvP,VarA,VarP
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 seedout
*     --- Writes a summary computing total biomass and extinctions.
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
      REAL*8 midRhoP,widthRhoP,midRhoA,widthRhoA
      REAL*8 midGammaP,widthGammaP,midGammaA,widthGammaA
      REAL*8 midHp,widthHp,midHa,widthHa
      REAL*8 midGa,widthGa,midGp,widthGp
      REAL*8 f0P,f0A
      REAL*8 Delta
      INTEGER Sa,Sp
      COMMON/Parameters/midNp,widthNp,midNa,widthNa,
     &midAlphaP,widthAlphaP,midAlphaA,widthAlphaA,
     &midBetaP,widthBetaP,midBetaA,widthBetaA,
     &midRhoP,widthRhoP,midRhoA,widthRhoA,
     &midGammaP,widthGammaP,midGammaA,widthGammaA,
     &midHa,widthHa,midHp,widthHp,
     &midGa,widthGa,midGp,widthGp,
     &f0P,f0A,Delta
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
      call Topology(sizeA,sizeP,sizeT,GammaP,GammaA,alive,degreePout,
     &     degreeAout,NestPout,NestAout,NestTout,ConnectOut,
     &     singleOut,excludedPout,excludedAout)
      AvDegreeP=0
      VarDegreeP=0
      ShannonP=0
      ShannonA=0
      tmp=Nfiles-1
      WRITE(unit(tmp),*) '#Final Plant Values'
      WRITE(unit(tmp),*) '#SpecieId, BiomassIn, BiomassOut, Alpha, DegreeIn, 
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
         WRITE(unit(tmp),*) i,Np(i),u(i),AlphaP(i),
     &degreeP(i),degreePout(i),alive(i)
      ENDDO
      ShannonP=EXP(-ShannonP)
      AvDegreeA=0
      VarDegreeA=0
      tmp=Nfiles
      WRITE(unit(tmp),*) '#Final Animals Values'
      WRITE(unit(tmp),*) '#SpecieId, BiomassIn, BiomassOut, Alpha, DegreeIn, 
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
         WRITE(unit(tmp),*) i,Na(i),u(sizeP+i),AlphaA(i),
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
      tmp=81
      WRITE(tmp,*) 'BiomassP ', BiomassP
      WRITE(tmp,*) 'BiomassA ',BiomassA
      WRITE(tmp,*) 'ExtinctP', ExtinctP
      WRITE(tmp,*) 'ExtinctA', ExtinctA
      WRITE(tmp,*) 'SurvivingP',sizeP-ExtinctP
      WRITE(tmp,*) 'SurvivingA',sizeA-ExtinctA
      WRITE(tmp,*) 'NestednessPin',NestP
      WRITE(tmp,*) 'NestednessAin',NestA
      WRITE(tmp,*) 'NestednessTin',NestT
      WRITE(tmp,*) 'NestednessPout',NestPout
      WRITE(tmp,*) 'NestednessAout',NestAout
      WRITE(tmp,*) 'NestednessTout',NestTout
      WRITE(tmp,*) 'ConnectanceIn', Connect
      WRITE(tmp,*) 'ConnectanceOut', ConnectOut
      WRITE(tmp,*) 'SinglesIn', single
      WRITE(tmp,*) 'SinglesOut', singleOut
      WRITE(tmp,*) 'SingleExtinct', SingleExtinct
      WRITE(tmp,*) 'excludedPin', excludedP
      WRITE(tmp,*) 'excludedPout', excludedPout
      WRITE(tmp,*) 'excludedAin', excludedA
      WRITE(tmp,*) 'excludedAout', excludedAout
      WRITE(tmp,*) 'NegAlphaP', NegAlphaP
      WRITE(tmp,*) 'NegAlphaA', NegAlphaA
      WRITE(tmp,*) 'AvAlphaP',AvP
      WRITE(tmp,*) 'VarAlphaP',VarP
      WRITE(tmp,*) 'AvAlphaA',AvA
      WRITE(tmp,*) 'VarAlphaA',VarA
      WRITE(tmp,*) 'AvDegreeInP',AvDegreeP
      WRITE(tmp,*) 'VarDegreeInP',VarDegreeP
      WRITE(tmp,*) 'AvDegreeInA',AvDegreeA
      WRITE(tmp,*) 'VarDegreeInA',VarDegreeA
      WRITE(tmp,*) 'Seed: ',seedout
      WRITE(tmp,*) 'TrueDiversityP ',ShannonP
      WRITE(tmp,*) 'TrueDiversityA ',ShannonA

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
