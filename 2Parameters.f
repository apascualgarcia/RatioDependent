
*     *********************************************
*     SUBROUTINE Read_Generate_Inputs
*     *********************************************                                                  
      
      SUBROUTINE Read_Generate_Inputs(path,Nfiles,NfilesIn,unit,today,now,seed,fixedPoint,
     &sizeA,sizeP,sizeT,readNp,readNa,readAlphaP,readAlphaA,
     &readBetaP,readBetaA,readGammaP,readGammaA,
     &readHp,readHa,readGp,readGa,BetaA,BetaP,GammaA,GammaP,Na,Np,u,AlphaA,AlphaP,
     &degreeP,degreeA,NestA,NestP,NestT,AdjA,AdjP,Connect,single,
     &excludedP,excludedA,NegAlphaP,NegAlphaA,hhA,hhP,seedout,AvA,AvP,VarA,VarP)
      IMPLICIT NONE
      INTEGER sizeA,sizeP,sizeT
      INTEGER unit(30),Nfiles,NfilesIn
      INTEGER*4 today(3), now(3)
      CHARACTER*40 path(30)
      REAL*8 seed,fixedPoint
      INTEGER readNp,readNa,readAlphaP,readAlphaA
      INTEGER readBetaP,readBetaA
      INTEGER readGammaP,readGammaA
      INTEGER readHp,readHa,readGp,readGa
*     This subroutine either reads or generates random parameters, 
*     according to the input metaparameters defined by the user. Parameters
*     are read if the metaparameter controling the correspondent distribution
*     is zero (no distribution should be centered around zero). In the case of alpha, if this parameter
*     is negative, these parameters will be generated according to
*     the constrained imposed by all the other parameters (solving eqs. at the fixed point).
*     Since the population abundances may scale other parameters, we must read it first.
*     Since we consider the possibility of generating alphas restricted
*     in accordance to the other parameters, we compute t at the end
*     in a separately subroutine
*     Finally, it computes some topological measures before running the simulation
*     at the routine topology
      LOGICAL alive(sizeT)
      INTEGER i,j,k,tmp,nkill
      REAL*8 scale,rnd,ran2,seedout,c
      INTEGER*8 idum
      REAL AdjA(sizeA,sizeP),AdjP(sizeP,sizeA)
      REAL degreeP(sizeP),degreeA(sizeA)
      INTEGER excludedP,excludedA,NegAlphaP,NegAlphaA
      REAL Connect,single
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA),gammaMin
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 NestA,NestP,NestT
      REAL*8 AvA,AvP,VarA,VarP
*     ...Commons 
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 rhoP,widthRhoP,rhoA,widthRhoA
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
     &midGa,widthGa,midGp,widthGp,f0P,f0A
      COMMON/Species/ Sa,Sp

*     --- Prepare random numbers and initialize the alive vector

c      seed=float(Now(3))**2*float(Now(2))+float(Now(1)) ! -(Seconds**2*Minutes+Hours)          
      seed=seed*1e8
      seedout=seed
      idum=-INT(seed)           ! We generate just one seed for all randomizations, as the total number of realizations is small
!PRINT *,idum,seed,' idum, seed' ! Debug
      
      DO i=1,Sp+Sa         
         alive(i)=.true.
      ENDDO

*     --- Read/generate vectorial tmp.
      tmp=1
      call ReadGenerateVec(unit,path,idum,Sp,tmp,readNp,midNp,widthNp,Np) ! Plant species abundances
      tmp=2
      call ReadGenerateVec(unit,path,idum,Sa,tmp,readNa,midNa,widthNa,Na) ! Animal species abundances
      tmp=9
      call ReadGenerateVec(unit,path,idum,Sp,tmp,readHp,midHp,widthHp,hhP) ! handling times opposite pool abundances
      tmp=10
      call ReadGenerateVec(unit,path,idum,Sa,tmp,readHa,midHa,widthHa,hhA)
      tmp=11
      call ReadGenerateVec(unit,path,idum,Sp,tmp,readGp,midGp,widthGp,ggP) ! handling times same pool abundances
      tmp=12
      call ReadGenerateVec(unit,path,idum,Sa,tmp,readGa,midGa,widthNGa,ggA)
      tmp=3
      call ReadGenerateTrixComp(unit,path,idum,Sp,tmp,readBetaP,midBetaP,midRhoP, ! Competition matrix, plants
     &widthBetaP,widthRhoP,BetaP)
      tmp=4
      call ReadGenerateTrixComp(unit,path,idum,Sa,tmp,readBetaA,midBetaA,midRhoA, ! Competition matrix, animals
     &widthBetaA,widthRhoA,BetaA)
      tmp=5
      call ReadGenerateTrixInt(unit,path,idum,Sp,Sa,tmp,readGammaP,midGammaP, ! Coupling matrix, plants
     &widthGammaP,GammaP)
      tmp=6
      call ReadGenerateTrixInt(unit,path,idum,Sa,Sp,tmp,readGammaA,midGammaA, ! Coupling matrix, animals
     &widthGammaA,GammaA)                     
      IF(fixedPoint .eq. 0)THEN ! Read or generate values for bare productivities
         tmp=7
         call ReadGenerateVec(unit,path,idum,Sp,tmp,readAlphaP,midAlphaP, ! Plants
     &widthAlphaP,AlphaP)
         tmp=8
         call ReadGenerateVec(unit,path,idum,Sa,tmp,readAlphaA,midAlphaA, ! Animals
     &widthAlphaA,AlphaA)
      ELSE
         ! Note that tmp=7,8 should be fixed within the subroutine
         call Consistent_Alpha(NfilesIn,unit,idum,sizeA,sizeP,sizeT, ! Estimate alpha at steady state
     &BetaA,BetaP,GammaA,GammaP,Na,Np,hhA,hhP,ggA,ggP,
     &AlphaA,AlphaP,AvA,AvP,VarA,VarP)
      ENDIF

*     --- Control parameters      
*     .... First, some quantities to compute from the coupling interaction matrix

            call Topology(sizeA,sizeP,sizeT,GammaP,GammaA,alive,degreeP,degreeA,NestP,NestA,
     &NestT,Connect,single,excludedP,excludedA) ! Compute some topological measures   
      
*     .... Next, we perform some additional controls on the alpha parameters, like the sign and variance
*     to control that you are working in the regime you want (e.g. facultative vs. obligatory)
*     mutualism, and to see the dependence with other parameters if you estimated them at a fixed point      
      AvP=0
      VarP=0
      NegAlphaP=0
      DO i=1,Sp
         AvP=AvP+AlphaP(i)
         VarP=VarP+AlphaP(i)**2
         IF(AlphaP(i).le.0) NegAlphaP=NegAlphaP+1
      ENDDO
      AvA=0
      VarA=0
      NegAlphaA=0
      DO i=1,Sa
         AvA=AvA+AlphaA(i)
         VarA=VarA+AlphaA(i)**2
         IF(AlphaA(i).le.0) NegAlphaA=NegAlphaA+1
      ENDDO
      AvP=AvP/float(Sp)
      VarP=VarP/float(Sp)-AvP**2
      AvA=AvA/float(Sa)
      VarA=VarA/float(Sa)-AvA**2
  
*     --- Finally,rewrite populations to simplify the notation
      k=0
      DO i=1,Sp
         k=k+1
         u(k)=Np(i)
      ENDDO
      DO i=1,Sa
         k=k+1
         u(k)=Na(i)
      ENDDO      
  
      END

*     *********************************************
*     SUBROUTINE ReadGenerateVec
*     *********************************************                                                  
 
      SUBROUTINE ReadGenerateVec(unit,path,idum,Sx,tmp,readX,midX,widthX,X)
      IMPLICIT NONE
      INTEGER unit(30),idum
      CHARACTER*40 path(30)
      INTEGER Sx,tmp
      REAL*8 midX,widthX,readX
*     This routine either reads a vector from file
*     identified by unit tmp, or it generates a random
*     one with cells following uniform distribution
*     centered at "midX" and width "widthX", which is a value
*     between 0 and 1 representing a width between 0 and midX. In the case in
*     which the matrix is read from file, it also randomizes
*     the values using the parameter widthX, but  now centered
*     at the value read. Then it  writes the vector in a new file and then returns it.      
      REAL*8 rnd,ran2,scale
      REAL*8 X(Sx)

      scale=widthX/2           ! Width scaling to use the modified ran2 function
      IF(readX.eq.1)THEN        ! Read from file
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sx
            rnd=ran2(idum)
            READ(unit(tmp),*) X(i)
            X(i)=X(i)*(1+rnd*scale)
            WRITE(unit(NfilesIn+tmp),*) X(i)
         ENDDO
      ELSE                      ! Generate parameters
         DO i=1,Sx
            rnd=ran2(idum)
            X(i)=midX*(1+rnd*scale)
            WRITE(unit(NfilesIn+tmp),*) X(i)
         ENDDO
      ENDIF
      CLOSE(unit(NfilesIn+tmp))

      RETURN
 770  CONTINUE
      PRINT *,'  * * Problems opening unit: ',unit(tmp),' >> ',path(tmp)
      call Warning() 
      END

      
*     *********************************************
*     SUBROUTINE ReadGenerateTrixComp
*     *********************************************                                                  
 
      SUBROUTINE ReadGenerateTrixComp(unit,path,idum,Sx,tmp,readM,midDiag,
     &readOff,midOff,widthDiag,widthOff,M)
      IMPLICIT NONE
      INTEGER unit(30),idum
      CHARACTER*40 path(30)
      INTEGER Sx,tmp,readM
      REAL*8 midDiag,widthDiag,midOff,widthOff
*     This routine either reads a competition matrix from file
*     identified by unit tmp, or it generates a random
*     one with cells following uniform distribution
*     centered at "midDiag" and width "widthDiag" for diagonal
*     cells and "minOff" and "widthOff" for off diagonal cells.
*     In the case in which the matrix is read from file, it also randomizes
*     the values using the parameters width, but  now centered
*     at the value read. Then it
*     writes the matrix in a new file and then returns it.
*     With respect to the routine aimed to read the interaction
*     matrix between pools, this one is a square matrix and has four parameters
*     while that ReadGenerateTrixInt considers a rectangular matrix with two param.
      REAL*8 rnd,ran2,scale
      REAL*8 M(Sx,Sx)
          
      IF(readM.eq.1)THEN  
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sx
            DO j=1,Sx
               rnd=ran2(idum)        
               READ(unit(tmp),*) M(i,j)
               IF(i.eq.j)THEN
                  scale=widthDiag/2
                  M(i,j)=M(i,j)*(1.0d0+rnd*scale)
               ELSE
                  scale=widthOff/2
                  M(i,j)=M(i,j)*(1.0d0+rnd*scale)
               ENDIF
               WRITE(unit(NfilesIn+tmp),*) M(i,j)
            ENDDO
         ENDDO
      ELSE
         DO i=1,Sx
            DO j=1,Sx
               rnd=ran2(idum)        
               IF(i.eq.j)THEN
                  scale=widthDiag/2
                  M(i,j)=midDiag*(1.0d0+rnd*scale)
               ELSE
                  scale=widthOff/2
                  M(i,j)=minOff*(1.0d0+rnd*scale)
               ENDIF
               WRITE(unit(NfilesIn+tmp),*) M(i,j)
            ENDDO
         ENDDO 
      ENDIF
      CLOSE(unit(NfilesIn+tmp))
      RETURN
 770  CONTINUE
      PRINT *,'  * * Problems opening unit: ',unit(tmp),' >> ',path(tmp)
      call Warning() 
      END

      
*     *********************************************
*     SUBROUTINE ReadGenerateTrixInt
*     *********************************************                                                  
 
      SUBROUTINE ReadGenerateTrixInt(unit,path,idum,Sx,Sy,tmp,readM,midM,
     &     widthM,M)
      IMPLICIT NONE
      INTEGER unit(30),idum
      CHARACTER*40 path(30)
      INTEGER Sx,Sy,tmp,readM
      REAL*8 midM,widthM
*     This routine either reads a matrix coupling pools of
*     species from a file,  identified by unit tmp, or it generates a random
*     one with cells following a uniform distribution
*     centered at "midM" and width "widthM".
*     In the case in which the matrix is read from file, it also randomizes
*     the values using the parameters width, but  now centered
*     at the value read. Then it
*     writes the matrix in a new file and then returns it.
*     With respect to the routine aimed to read the interaction
*     matrix between pools, this one is a square matrix and has four parameters
*     while that ReadGenerateTrixInt considers a rectangular matrix with two param.
      REAL*8 rnd,ran2,scale
      REAL*8 M(Sx,Sy)
      
      scale=widthM/2         
      IF(readM.eq.1)THEN  
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sx
            DO j=1,Sy
               rnd=ran2(idum)        
               READ(unit(tmp),*) M(i,j)
               M(i,j)=M(i,j)*(1.0d0+rnd*scale)
               WRITE(unit(NfilesIn+tmp),*) M(i,j)
            ENDDO
         ENDDO
      ELSE
         DO i=1,Sx
            DO j=1,Sy
               rnd=ran2(idum)        
               M(i,j)=midM*(1.0d0+rnd*scale)
               WRITE(unit(NfilesIn+tmp),*) M(i,j)
            ENDDO
         ENDDO 
      ENDIF
      CLOSE(unit(NfilesIn+tmp))
      RETURN
 770  CONTINUE
      PRINT *,'  * * Problems opening unit: ',unit(tmp),' >> ',path(tmp)
      call Warning() 
      END

      
*     *********************************************
*     SUBROUTINE Consistent_Alpha
*     *********************************************                                                  
      
      SUBROUTINE Consistent_Alpha(NfilesIn,unit,idum,sizeA,sizeP,sizeT,
     &BetaA,BetaP,GammaA,GammaP,Na,Np,hhA,hhP,ggA,ggP,
     &AlphaA,AlphaP,AvA,AvP,VarA,VarP)
      IMPLICIT NONE
      INTEGER sizeA,sizeP,sizeT
      INTEGER*8 idum
      INTEGER NfilesIn,unit(30)
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 ggA(sizeA),ggP(sizeP)
*     This subroutine determines the parameter alpha compatible
*     with all the other parameters in order to get negative
*     alpha for animals and positive for plants (obligate mutualism).
*     We expect to get stable systems using this criteria.
      INTEGER i,j,k,tmp
      REAL*8 rnd,Comp,Int,ran2
      REAL*8 HollingA(sizeA),HollingP(sizeP)
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 AvA,AvP,VarA,VarP
*     ...Commons
      REAL*8 midNp,widthNp,midNa,widthNa
      REAL*8 midAlphaP,widthAlphaP,midAlphaA,widthAlphaA
      REAL*8 midBetaP,widthBetaP,midBetaA,widthBetaA
      REAL*8 rhoP,widthRhoP,rhoA,widthRhoA
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
     &midGa,widthGa,midGp,widthGp,f0P,f0A
      COMMON/Species/ Sa,Sp
      
      DO i=1,Sp                 ! Compute first a Holling like term to saturate the ODEs
         HollingP(i)=0.0d0
         GollingP(i)=0.0d0 ! New Holling term, a hybrid between Holling and Gollum
         IF(hhP(i).gt.0)THEN
            DO k=1,Sa
               HollingP(i)=HollingP(i)+GammaP(i,k)*Na(k)
            ENDDO
         ENDIF
         IF(ggP(i).gt.0)THEN
            DO k=1,Sa
               GollingP(i)=GollingP(i)+GammaA(k,i)*Na(k)
            ENDDO
         ENDIF
      ENDDO
      DO i=1,Sa
         HollingA(i)=0.0d0
         GollingA(i)=0.0d0
         IF(hhA(i).gt.0)THEN
            DO k=1,Sp
               HollingA(i)=HollingA(i)+GammaA(i,k)*Np(k)
            ENDDO
         ENDIF
         IF(ggA(i).gt.0)THEN
            DO k=1,Sp
               GollingA(i)=GollingA(i)+GammaP(k,i)*Np(k)
            ENDDO
         ENDIF
      ENDDO     
      

      tmp=7
      DO i=1,Sp                 ! Compute first Holling term to saturate the ODEs
         Comp=0.0d0
         Int=0.0d0
         rnd=ran2(idum)
         DO j=1,Sp
            Comp=Comp+BetaP(i,j)*Np(j)
         ENDDO         
         DO k=1,Sa
            Int=Int+GammaP(i,k)*Na(k)/(f0P+hhP(i)*HollingP(i)+ggA(k)*GollingA(k))
         ENDDO         
         AlphaP(i)=(Comp-Int)*(1+2*Delta*rnd) ! The rnd scaling factor is 2 here

         WRITE(unit(NfilesIn+tmp),*) AlphaP(i),Comp,Int
      ENDDO
      CLOSE(unit(NfilesIn+tmp))
      tmp=8
      DO i=1,Sa                 ! Compute first Holling term to saturate the ODEs
         Comp=0.0d0
         Int=0.0d0
         rnd=ran2(idum)
         DO j=1,Sa
            Comp=Comp+BetaA(i,j)*Na(j)
         ENDDO
         DO k=1,Sp
            Int=Int+GammaA(i,k)*Np(k)/(f0A+hhA(i)*HollingA(i)+ggP(k)*GollingP(k))
         ENDDO         
         AlphaA(i)=(Comp-Int)*(1+2*Delta*rnd) ! The rnd scaling factor is 2 here
         WRITE(unit(NfilesIn+tmp),*) AlphaA(i),Comp,Int
      ENDDO
      CLOSE(unit(NfilesIn+tmp))

      RETURN
      END

*     *********************************************
*     SUBROUTINE Topology
*     *********************************************                                                  
 
      SUBROUTINE Topology(sizeA,sizeP,sizeT,AdjP,AdjA,alive,degreeP,
     &degreeA,NestP,NestA,NestT,Connect,single,excludedP,excludedA)
      IMPLICIT NONE
      INTEGER sizeA,sizeP,sizeT
      REAL AdjA(sizeA,sizeP),AdjP(sizeP,sizeA)
      LOGICAL alive(sizeT)
*     This subroutine computes some topological properties. It has as input the
*     adjacency matrices for plants and animals and a logical array (alive) which
*     allow you to reuse the routine when some species go extinct (it must be
*     defined as true if the specie has went extincted). It returns the
*     degrees, the nestedness, the connectance, the number of species
*     connected just once (single), and the number of species not connected (excluded)      
      LOGICAL compute
      INTEGER i,j,k,last
      REAL Ni,Nj,Nij,pairs,minNiNj,minTmp
      REAL Connect,NaliveP,NaliveA
      REAL degreeP(sizeP),degreeA(sizeA)
      INTEGER excludedP,excludedA
      REAL single
      REAL*8 NestA,NestP,NestT

      NestP=0
      minNiNj=0
      pairs=0
      Connect=0
      excludedP=0
      single=0
      NaliveP=0
      last=0
      DO i=1,sizeP-1 ! For each plant
         degreeP(i)=0
         compute=.true.         ! Control that degrees and connectance are computed just once
         IF(alive(i).eqv..true.)THEN ! If it is alive           
            NaliveP=NaliveP+1
            last=i !  This will be the last element alive if i=sizeP has died
            DO j=i+1,sizeP ! Look for another plant
               Ni=0
               Nj=0
               Nij=0
               IF(alive(j).eqv..true.)THEN ! If it is also alive           
                  DO k=1,sizeA ! Look for the animals
                     IF(alive(sizeP+k).eqv..true.)THEN ! If the animal is alive
                        IF(AdjP(i,k).gt.0)THEN ! Look whether the first plant interact
                           Ni=Ni+1 ! And sum up
                           IF(AdjP(j,k).gt.0)THEN ! If also the other plant interact
                              Nj=Nj+1 ! Sum up
                              Nij=Nij+1 ! And count a shared interaction
                           ENDIF
                        ELSEIF(AdjP(j,k).gt.0)THEN ! If the first plant does not interact maybe the second does     
                           Nj=Nj+1                         
                        ENDIF
                     ENDIF
                  ENDDO
                  IF((Ni.gt.0).and.(Nj.gt.0))THEN ! This is true unless the specie has been disconnected
                     pairs=pairs+1
                     minNiNj=minNiNj+MIN(Ni,Nj)
                     NestP=NestP+Nij
                  ENDIF
                  IF(compute.eqv..true.)THEN ! This must be done just once
                     degreeP(i)=Ni 
                     Connect=Connect+Ni
                     IF(degreeP(i).eq.0) excludedP=excludedP+1
                     IF(degreeP(i).eq.1) single=single+1
                     compute=.false.
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IF(alive(sizeP).eqv..true.)THEN ! Check whether the last element is alive, otherwise it will be the one labelled as "last"
         last=sizeP
         NaliveP=NaliveP+1 ! If it was not the last one it was already counted
      ELSE
         degreeP(sizeP)=0
      ENDIF
c$$$      !PRINT *, '> debug',last,NaliveP
c$$$      IF(NaliveP.eq.0)THEN
c$$$         DO i=1,sizeP
c$$$            PRINT *, i,alive(i),' debug alive i'
c$$$         ENDDO
c$$$      ENDIF
      Ni=0
      IF(last.ne.0)THEN
         DO k=1,sizeA
!PRINT *, '> debug',last,k,NaliveP
            IF(adjP(last,k).gt.0) Ni=Ni+1
         ENDDO       
         degreeP(last)=Ni

         Connect=Connect+Ni
         IF(degreeP(last).eq.0) excludedP=excludedP+1
         IF(degreeP(last).eq.1) single=single+1
      ENDIF
      NestT=NestP
      minTmp=minNiNj
      IF(minNiNj.ne.0) NestP=NestP/minNiNj             

      NestA=0
      NaliveA=0
      minNiNj=0
      pairs=0
      excludedA=0
      DO i=1,sizeA-1 ! For each animal
         degreeA(i)=0
         compute=.true.
         IF(alive(sizeP+i).eqv..true.)THEN ! If it is alive
            NaliveA=NaliveA+1
            last=i ! This will be the last element alive if i=sizeA has died
            DO j=i+1,sizeA               ! Compare with the other animals
               IF(alive(sizeP+j).eqv..true.)THEN ! And if the other is also alive
                  Ni=0
                  Nj=0
                  Nij=0
                  DO k=1,sizeP ! Look at the plants
                     IF(alive(k).eqv..true.)THEN   ! If these plants are alive                                     
                        IF(AdjA(i,k).gt.0)THEN ! If the first animal is connected with this plant
                           Ni=Ni+1 ! Count it as an interaction
                           IF(AdjA(j,k).gt.0)THEN ! But if the other animal is also connected
                              Nj=Nj+1 ! Count as an interaction for the other animal
                              Nij=Nij+1 ! And as a shared interactions
                           ENDIF
                        ELSEIF(AdjA(j,k).gt.0)THEN ! If it is not connected with the plant but it is the other animal
                           Nj=Nj+1 ! Take it into account
                        ENDIF
                     ENDIF
                  ENDDO
                  IF((Ni.gt.0).and.(Nj.gt.0))THEN ! If both have partners
                     pairs=pairs+1
                     NestA=NestA+Nij ! Compute the nestedness
                     minNiNj=minNiNj+MIN(Ni,Nj)
                  ENDIF
                  IF(compute.eqv..true.)THEN ! For each i we compute the degrees just once
                     degreeA(i)=Ni 
                     IF(degreeA(i).eq.0) excludedA=excludedA+1
                     IF(degreeA(i).eq.1) single=single+1
                     compute=.false.
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      IF(alive(sizeP+sizeA).eqv..true.)THEN
         last=sizeA             ! The last one is also alive 
         NaliveA=NaliveA+1
      ELSE
         degreeA(sizeA)=0
      ENDIF
      Ni=0
      IF(last.ne.0)THEN
         DO k=1,sizeP
            IF(adjA(last,k).gt.0) Ni=Ni+1
         ENDDO
         degreeA(last)=Ni
         IF(degreeA(last).eq.0) excludedA=excludedA+1
         IF(degreeA(last).eq.1) single=single+1
      ENDIF
      
      NestT=NestT+NestA
      IF(minNiNj.ne.0) NestA=NestA/minNiNj       

      minTmp=minTmp+minNiNj
      IF(minTmp.ne.0) NestT=NestT/minTmp
      IF((NaliveP.gt.0).and.(NaliveA.gt.0))THEN
         Connect=Connect/(NaliveP*NaliveA)
      ELSE
         Connect=0
      ENDIF

      RETURN
      END

*     *********************************************
*     FUNCTION ran2  
*     *********************************************          \xg\N = 0 \xr\f{} = 0.23                                     

      REAL*8 FUNCTION ran2(idum)  
      IMPLICIT NONE
*     Improved Minimal Standard rand function with two additional 
*     shufflings. For further details see Numerical Recepies for
*     Fortran 77, chapter 7.      
      INTEGER*8 idum
      INTEGER IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     &     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     &     IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

*     Long period (> 21018) random number generator of L'Ecuyer with Bays-Durham shuffle
*     and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive
*     of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not
*     alter idum between successive deviates in a sequence. RNMX should approximate the largest
*     floating value that is less than 1. 
*     ---Modified July 2012 (A.P-G.) I have modified the function ran2 to 
*     give a number between [-0.5,0.5]
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      
      IF (idum.le.0) THEN       ! Initialize.
         idum=max(-idum,1)      ! Be sure to prevent idum = 0.
         idum2=idum
         DO j=NTAB+8,1,-1       ! Load the shue table (after 8 warm-ups).
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            IF (idum.lt.0) idum=idum+IM1
            IF (j.le.NTAB) iv(j)=idum
         ENDDO
         iy=iv(1)
      ENDIF

      k=idum/IQ1                ! Start here when not initializing.
      idum=IA1*(idum-k*IQ1)-k*IR1 ! Compute idum=mod(IA1*idum,IM1) without overif
      IF(idum.lt.0) idum=idum+IM1 ! flows by Schrage's method.
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2 ! Compute idum2=mod(IA2*idum2,IM2)  likewise.
      IF(idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV               ! Will be in the range 1:NTAB.
      iy=iv(j)-idum2            ! Here idum is shuffled, idum and idum2 are combined to generate output
      iv(j)=idum                   
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)      ! Because users don't expect endpoint values.
      IF(ran2.ge.0.5d0)THEN ! Map the interval (0-1) to (-0.5,0.5)
         ran2=ran2-0.5d0         
      ELSE
         ran2=-ran2
      ENDIF
      return
      END
            
