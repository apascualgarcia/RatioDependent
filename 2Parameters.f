
*     *********************************************
*     SUBROUTINE Read_Generate_Inputs
*     *********************************************                                                  
      
      SUBROUTINE Read_Generate_Inputs(path,Nfiles,NfilesIn,unit,today,now,seed,
     &sizeA,sizeP,sizeT,BetaA,BetaP,GammaA,GammaP,Na,Np,u,AlphaA,AlphaP
     &,degreeP,degreeA,NestA,NestP,NestT,AdjA,AdjP,Connect,single,
     &excludedP,excludedA,NegAlphaP,NegAlphaA,hhA,hhP,seedout,AvA,AvP,VarA,VarP)
      IMPLICIT NONE
      INTEGER sizeA,sizeP,sizeT
      INTEGER unit(30),Nfiles,NfilesIn
      INTEGER*4 today(3), now(3)
      CHARACTER*40 path(30)
      REAL*8 seed
*     This subroutine either reads or generates random parameters, 
*     according to the input metaparameters defined by the user. Parameters
*     are readed if the metaparameter controling the correspondent distribution
*     is zero (no distribution should be centered around zero). In the case of alpha, if this parameter
*     is negative, these parameters will be generated according to
*     the constrained imposed by all the other parameters (solving eqs. at the fixed point).
*     Since the population abundances may scale other parameters, we must read it first.
*     Since we consider the possibility of generating alphas restricted
*     in accordance to the other parameters, we compute t at the end
*     in a separately subroutine
*     Finally, it computes some topological measures before running the simulation
*     at the routine topology
      LOGICAL alive(sizeT),binary
      INTEGER i,j,k,tmp,nkill
      REAL*8 scale,rnd,ran2,seedout,c
      INTEGER*8 idum
      REAL AdjA(sizeA,sizeP),AdjP(sizeP,sizeA)
      REAL degreeP(sizeP),degreeA(sizeA)
      INTEGER excludedP,excludedA,NegAlphaP,NegAlphaA
      REAL Connect,single
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 bup,blow,rhomaxA,rhomaxP
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA),gammaMin
      REAL*8 Gamma0P,Gamma0A
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT),Nrate,InvNrate
      REAL*8 hhA(sizeA),hhP(sizeP),hmax,hmaxA,hmaxP
      REAL*8 NestA,NestP,NestT
      REAL*8 AvA,AvP,VarA,VarP
*     ...Commons 
      REAL*8 midNa,widthNa,midNp,widthNp
      REAL*8 midAlpha,widthAlpha,Beta0,widthBeta
      REAL*8 rhoA,rhoP
      REAL*8 Gamma0,widthGamma,Delta
      REAL*8 hA,hP
      INTEGER Sa,Sp
      COMMON/Parameters/ midNa,widthNa,midNp,widthNp,
     &midAlpha,widthAlpha,Beta0,widthBeta,rhoA,rhoP,
     &Gamma0,widthGamma,hA,hP,Delta
      COMMON/Species/ Sa,Sp

*     --- Prepare random numbers

c      seed=float(Now(3))**2*float(Now(2))+float(Now(1)) ! -(Seconds**2*Minutes+Hours)          
      seed=seed*1e8
      seedout=seed
      idum=-INT(seed)           ! We generate just one seed for all randomizations, as the total number of realizations is small
      !PRINT *,idum,seed,' idum, seed' ! Debug

*     --- Compute topological measures

      binary=.true. 
      DO i=1,Sp  
         alive(i)=.true.
         READ(unit(Nfiles+1),*) (AdjP(i,k), k=1,Sa) ! Read the adjacency matrix in implicit for
         DO k=1,Sa              ! Compute degrees and transpose the matrix
            IF(i.eq.1)THEN
               alive(Sp+k)=.true.
            ENDIF
            AdjA(k,i)=AdjP(i,k)
            IF(binary.eqv..true.)THEN ! Check if your matrix is binary
               IF((AdjA(k,i).ne.0).or.(AdjA(k,i).ne.1)) binary=.false.
            ENDIF
         ENDDO
      ENDDO
      call Topology(sizeA,sizeP,sizeT,AdjP,AdjA,alive,degreeP,degreeA,NestP,NestA,
     &NestT,Connect,single,excludedP,excludedA) ! Compute some topological measures

*     --- Generate parameters to compute self consistent handling times or rho

      c=0.75d0                  ! Scales how far we set h or rho from hmax and rhomax
      scale=widthBeta/2
      bup=Beta0*(1+scale)
      blow=Beta0*(1-scale)
      gammaMin=0.01             ! BE AWARE of this choice (Obligatory mutualism)
      hmax=100d0
      hmaxA=1.0d0/(bup*((Sa-1)*rhoA+1)) ! Note that rho cannot be negative here
      hmaxP=1.0d0/(bup*((Sp-1)*rhoP+1))   
      hmaxA=c*hmaxA             ! First solution (when h must be below hmax)
      hmaxP=c*hmaxP
      rhomaxP=1.0d0/(Sp-1)*(1.0d0/(hP*bup)-1) ! Note that h cannot be negative here
      rhomaxA=1.0d0/(Sa-1)*(1.0d0/(hA*bup)-1)
      rhomaxP=c*rhomaxP
      rhomaxA=c*rhomaxA
      
      DO i=1,Sp
         IF(hP.gt.0)THEN
            hhP(i)=hP
         ELSE
            hhP(i)=hmaxP        !-InvNrate/(gammaMin*degreeP(i))
            !PRINT *,i,hhP(i),degreeP(i),' hP debug'
         ENDIF
      ENDDO
      DO i=1,Sa
         IF(hA.gt.0)THEN
            hhA(i)=hA
         ELSE
            hhA(i)=hmaxA        !-InvNrate/(gammaMin*degreeA(i))
            IF(hhA(i).lt.hmax) hmax=hhA(i)
            !PRINT *,i,hhA(i),degreeA(i),' hA debug'
         ENDIF
      ENDDO

      !IF(hP.le.0) hP=hmaxP               ! For Nrate and print it as output
      !IF(hA.le.0) hA=hmaxA               ! Print it as output and use it for Nrate
      !IF(rhoA.le.0) rhoA=rhomaxA
      !IF(rhoP.le.0) rhoP=rhomaxP       

*     --- Compute the rate of biomasses needed for an obligatory scenario
            
      IF(midNp.lt.0)THEN        ! We want abundances consistent with h
         hA=hmaxA
         Nrate=(bup*(1+(Sa-1)*rhoA))**2 ! A conservative increase against a realization of biomasses such that Np/Na=Np/Na*(1-0.15)/(1+0.15)
     &        /(gammaMin*blow*(1-hA*(bup*(1+(Sa-1)*rhoA))))**2 ! In the denominator, blow would be actually clow, the mutualistic parameter
         Nrate=Nrate*(1.0d0+0.8d0)
         midNp=Nrate
         midNa=1.0d0
      ENDIF
*     --- Read/generate values for the populations.
      !PRINT*, '>> DEBUG: MidNp',midNp, ' path(1): ',path(1),' NfilesIn ',NfilesIn
      IF(midNp.eq.0)THEN
         tmp=1
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         scale=widthNp/2
         DO i=1,Sp
            rnd=ran2(idum)
            READ(unit(tmp),*) Np(i)
            !Np(i)=100*Np(i) ! Additional abundances parameter scaling
            Np(i)=Np(i)*(1+rnd*scale)
            midNp=midNp+Np(i)
         ENDDO
         midNp=midNp/float(Sp) ! We will need this value to normalize Beta0 and gamma
      ELSE
         scale=widthNp/2
         DO i=1,Sp
            rnd=ran2(idum)
            Np(i)=midNp*(1+rnd*scale)
             !print *,i,Np(i),Nrate,midNp,scale,' i, Np(i), Nrate, debug'
         ENDDO
      ENDIF  
 
      IF(midNa.eq.0)THEN
         tmp=2
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         scale=widthNa/2    ! Scaling to use the modified ran2 function
         DO i=1,Sa
            rnd=ran2(idum)
            READ(unit(tmp),*) Na(i)
            !Na(i)=100*Na(i) ! Additional abundances parameter scaling
            Na(i)=Na(i)*(1+rnd*scale)
            midNa=midNa+Na(i)
         ENDDO
         midNa=midNa/float(Sa)
      ELSE
         scale=widthNa/2    ! Scaling to use the modified ran2 function
         DO i=1,Sa
            rnd=ran2(idum)      ! Modified ran2 function to get a number between [-0.5-0.5]
            Na(i)=midNa*(1+rnd*scale)
            ! print *,i,Na(i),midNa,scale,' i, Na(i), debug'
         ENDDO
      ENDIF

*     --- Cancel the scaling effects, comment otherwise

      midNa=1
      midNp=1

*     --- Read/generate values for competition parameters, scaled with the average of the populations for the latter
        
      IF(Beta0.eq.0)THEN  
         tmp=3
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sp
            DO j=1,Sp
               READ(unit(tmp),*) BetaP(i,j)
               WRITE(unit(NfilesIn+tmp),*) BetaP(i,j)
            ENDDO
         ENDDO
         CLOSE(unit(NfilesIn+tmp))

         tmp=4
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sa
            DO j=1,Sa
               READ(unit(tmp),*) BetaA(i,j)
               WRITE(unit(NfilesIn+tmp),*) BetaP(i,j)
            ENDDO
         ENDDO
         CLOSE(unit(NfilesIn+tmp))

      ELSE
         scale=widthBeta/2
         tmp=NfilesIn+3
         DO i=1,Sp
            DO j=1,Sp
               rnd=ran2(idum)        
               IF(i.eq.j)THEN
                  BetaP(i,j)=Beta0*(1.0d0+rnd*scale)/midNp
               ELSE
                  BetaP(i,j)=rhoP*Beta0*(1.0d0+rnd*scale)/midNp
               ENDIF
               WRITE(unit(tmp),*) BetaP(i,j)
            ENDDO
         ENDDO 
         CLOSE(unit(tmp))
         tmp=NfilesIn+4
         DO i=1,Sa
            DO j=1,Sa
               rnd=ran2(idum)
               IF(i.eq.j)THEN
                  BetaA(i,j)=Beta0*(1.0d0+rnd*scale)/midNa
               ELSE
                  BetaA(i,j)=rhoA*Beta0*(1.0d0+rnd*scale)/midNa
               ENDIF
               WRITE(unit(tmp),*) BetaA(i,j)
            ENDDO
         ENDDO
         CLOSE(unit(tmp))
      ENDIF

*     --- Read/generate values for mutualistic parameters, scaled with the average of the populations for the latter, and the adjacency matrix

      IF(Gamma0.lt.-10000)THEN  ! You need a very low negative value to read the matrices 
         Gamma0=1               ! From now on Gamma0 is just a control parameter, and it is >0 if the system is mutualistic
         tmp=5
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sp
            DO k=1,Sa
               READ(unit(tmp),*) GammaP(i,k)
               WRITE(unit(NfilesIn+tmp),*) GammaP(i,k)
            ENDDO
         ENDDO
         CLOSE(unit(NfilesIn+tmp))
         tmp=6
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sa
            DO k=1,Sp
               READ(unit(tmp),*) GammaA(i,k)
               WRITE(unit(NfilesIn+tmp),*) GammaA(i,k)
            ENDDO            
         ENDDO
         CLOSE(unit(NfilesIn+tmp))
      ELSE
         IF(Gamma0.lt.0)THEN    ! We select here two possible intermediate scenario
!     the first scenario has an increased value of alpha equivalent to the mutualistic strength but without the feedback. To model
!     this case you need to keep Gamma0 with a negative value, what will be a control parameter, and set Gamma0A and Gamma0P
!     to the absolute Gamma0 value. The second scenario sets Gamma0P to zero, and  Gamma0A=ABS(Gamma0), setting also Gamma0 as a positive
!     value in order to proceed as if it were a mutualistic network (again, it works as a control value in the following).
            Gamma0=ABS(Gamma0) ! 
            Gamma0P=0           !ABS(Gamma0) ! the mutualistic strength but without the feedback            
            Gamma0A=Gamma0      !ABS(Gamma0)            
         ELSEIF(Gamma0.gt.0)THEN ! We just need to keep the gamma0 values
            Gamma0P=Gamma0
            Gamma0A=Gamma0
         ENDIF
         scale=widthGamma/2  
         tmp=NfilesIn+5
         DO i=1,Sp
            DO k=1,Sa
               rnd=ran2(idum) ! Modified ran2 function to get also a sign
               IF(binary.eqv..true.)THEN ! We need a randomization of the parameters
                  GammaP(i,k)=AdjP(i,k)*Gamma0P*(1.0d0+rnd*scale)/SQRT(midNa*midNp)
               ELSE             ! The parameters are already randomized or with a structure that we want to conserve
                  GammaP(i,k)=AdjP(i,k)*Gamma0P/SQRT(midNa*midNp)
               ENDIF
               WRITE(unit(tmp),*) GammaP(i,k)
            ENDDO
         ENDDO 
         CLOSE(unit(tmp))
         tmp=NfilesIn+6
         DO i=1,Sa
            DO k=1,Sp
               rnd=ran2(idum)
               IF(binary.eqv..true.)THEN ! We need a randomization of the parameters
                  GammaA(i,k)=AdjA(i,k)*Gamma0A*(1.0d0+rnd*scale)/SQRT(midNa*midNp)
               ELSE             ! The parameters are already randomized    
                  GammaA(i,k)=AdjA(i,k)*Gamma0A/SQRT(midNa*midNp)
               ENDIF
               WRITE(unit(tmp),*) GammaA(i,k)
            ENDDO
         ENDDO
         CLOSE(unit(tmp))
      ENDIF

*     --- Read/generate values for bare productivities

c      PRINT *, '   -- Alphas..  '
      IF(midAlpha.eq.0)THEN   
         tmp=7
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         AvP=0
         VarP=0
         AvA=0
         VarA=0
         DO i=1,Sp
            rnd=ran2(idum)
            READ(unit(tmp),*) AlphaP(i) 
            AlphaP(i)=AlphaP(i)*(1+2*Delta*rnd) ! Careful with this, feasibility is not guaranteed, just adding extra noise
            AvP=AvP+AlphaP(i)
            VarP=VarP+AlphaP(i)**2
         ENDDO
         tmp=8
         OPEN(UNIT=unit(tmp),STATUS='OLD',ERR=770,FILE=path(tmp))
         DO i=1,Sa
            rnd=ran2(idum)
            READ(unit(tmp),*) AlphaA(i)
            AlphaA(i)=AlphaA(i)*(1+2*Delta*rnd)
            AvA=AvA+AlphaA(i)
            VarA=VarA+AlphaA(i)**2
         ENDDO
         AvP=AvP/float(Sp)
         VarP=VarP/float(Sp)-AvP**2
         AvA=AvA/float(Sa)
         VarA=VarA/float(Sa)-AvA**2
      ELSEIF(midAlpha.gt.0)THEN
         scale=widthAlpha/2
         tmp=NfilesIn+7
         AvP=0
         VarP=0
         AvA=0
         VarA=0
         DO i=1,Sp
            rnd=ran2(idum)
            AlphaP(i)=midAlpha*(1+2*Delta*rnd)!(1.0d0+rnd*scale) ! Same happens here
            AvP=AvP+AlphaP(i)
            VarP=VarP+AlphaP(i)**2
            !WRITE(unit(tmp),*) AlphaP(i)
         ENDDO
         tmp=NfilesIn+8
         DO i=1,Sa
            rnd=ran2(idum)
            AlphaA(i)=midAlpha*(1+2*Delta*rnd)!(1.0d0+rnd*scale)
            AvA=AvA+AlphaA(i)
            VarA=VarA+AlphaA(i)**2
            !WRITE(unit(tmp),*) AlphaA(i)
         ENDDO
         AvP=AvP/float(Sp)
         VarP=VarP/float(Sp)-AvP**2
         AvA=AvA/float(Sa)
         VarA=VarA/float(Sa)-AvA**2
      ELSE                      ! Generate alpha according to the constrains imposed by the other parameters (feasibility imposed)
c         PRINT *, '   --- Looking for consistent alpha.. '
         call Consistent_Alpha(NfilesIn,unit,idum,sizeA,sizeP,sizeT,
     &BetaA,BetaP,GammaA,GammaP,Na,Np,hhA,hhP,AlphaA,AlphaP,AvA,AvP,VarA,VarP)
      ENDIF
*     --- Further optional operations (to be compiled
*     -- Select an specie and kill it:
c$$$      rnd=ran2(idum)           ! Select an specie and kill it!
c$$$      nkill=INT((rnd+0.5)*Sp)
c$$$      Np(nkill)=0.0d0
c$$$      AlphaP(nkill)=0

*     -- Add noise to the abundances after you reach the fixed point (global stability test)
c$$$      scale=0.7/2
c$$$      DO i=1,Sp
c$$$         rnd=ran2(idum)
c$$$         Np(i)=Np(i)*(1+rnd*scale)
c$$$      ENDDO
c$$$      scale=0.7/2
c$$$      DO i=1,Sa
c$$$         rnd=ran2(idum)
c$$$         Na(i)=Na(i)*(1+rnd*scale)
c$$$      ENDDO

*     --- Check the regime you are working in
      NegAlphaP=0
      DO i=1,Sp
         IF(AlphaP(i).le.0) NegAlphaP=NegAlphaP+1
      ENDDO
      NegAlphaA=0
      DO i=1,Sa
         IF(AlphaA(i).le.0) NegAlphaA=NegAlphaA+1
      ENDDO   
*     --- Rewrite populations to simplify the notation
      k=0
      DO i=1,Sp
         k=k+1
         u(k)=Np(i)
      ENDDO
      DO i=1,Sa
         k=k+1
         u(k)=Na(i)
      ENDDO      
c      PRINT *, '<< Parameters generated/readed!'
      RETURN
 770  CONTINUE
      PRINT *,'  * * Problems opening unit: ',unit(tmp),' >> ',path(tmp)
      call Warning()      
      END

*     *********************************************
*     SUBROUTINE Consistent_Alpha
*     *********************************************                                                  
      
      SUBROUTINE Consistent_Alpha(NfilesIn,unit,idum,sizeA,sizeP,sizeT,
     &BetaA,BetaP,GammaA,GammaP,Na,Np,hhA,hhP,AlphaA,AlphaP,AvA,AvP,VarA,VarP)
      IMPLICIT NONE
      INTEGER sizeA,sizeP,sizeT
      INTEGER*8 idum
      INTEGER NfilesIn,unit(30)
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 hhA(sizeA),hhP(sizeP)
*     This subroutine determines the parameter alpha compatible
*     with all the other parameters in order to get negative
*     alpha for animals and positive for plants (obligate mutualism).
*     We expect to get stable systems using this criteria.
      INTEGER i,j,k,tmp
      REAL*8 rnd,Comp,Mut,ran2
      REAL*8 HollingA(sizeA),HollingP(sizeP)
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 AvA,AvP,VarA,VarP
*     ...Commons
      REAL*8 midNa,widthNa,midNp,widthNp
      REAL*8 midAlpha,widthAlpha,Beta0,widthBeta
      REAL*8 rhoA,rhoP
      REAL*8 Gamma0,widthGamma,Delta
      REAL*8 hA,hP
      INTEGER Sa,Sp
      COMMON/Parameters/ midNa,widthNa,midNp,widthNp,
     &midAlpha,widthAlpha,Beta0,widthBeta,rhoA,rhoP,
     &Gamma0,widthGamma,hA,hP,Delta
      COMMON/Species/ Sa,Sp

      IF(Gamma0.ne.0)THEN
         DO i=1,Sp              ! Compute first a Holling like term to saturate the ODEs
            HollingP(i)=0.0d0
            DO k=1,Sa
               HollingP(i)=HollingP(i)+GammaP(i,k)*Na(k)
            ENDDO
         ENDDO
         DO i=1,Sa
            HollingA(i)=0.0d0
            DO k=1,Sp
               HollingA(i)=HollingA(i)+GammaA(i,k)*Np(k)
            ENDDO
         ENDDO     
      ENDIF
      AvP=0
      VarP=0
      tmp=NfilesIn+7
      DO i=1,Sp                 ! Compute first Holling term to saturate the ODEs
         Comp=0.0d0
         Mut=0.0d0
         rnd=ran2(idum)
         DO j=1,Sp
            Comp=Comp+BetaP(i,j)*Np(j)
         ENDDO
         IF(Gamma0.ne.0)THEN
            DO k=1,Sa
               Mut=Mut+GammaP(i,k)*Na(k)/(1+hhP(i)*HollingP(i))
            ENDDO
         ENDIF
         AlphaP(i)=(Comp-Mut)*(1+2*Delta*rnd) ! The rnd scaling factor is 2 here
         IF(Gamma0.lt.0)THEN    ! A second possible type of competition with Gamma0=0 and Alphas increased by Mut
            AlphaP(i)=AlphaP(i)+Mut
         ENDIF
         AvP=AvP+AlphaP(i)
         VarP=VarP+AlphaP(i)**2
         WRITE(unit(tmp),*) AlphaP(i),Comp,Mut
      ENDDO
      CLOSE(unit(tmp))
      AvA=0
      VarA=0
      tmp=NfilesIn+8
      DO i=1,Sa                 ! Compute first Holling term to saturate the ODEs
         Comp=0.0d0
         Mut=0.0d0
         rnd=ran2(idum)
         DO j=1,Sa
            Comp=Comp+BetaA(i,j)*Na(j)
         ENDDO
         IF(Gamma0.ne.0)THEN    ! There is no pure competition
            DO k=1,Sp
               Mut=Mut+GammaA(i,k)*Np(k)/(1+hhA(i)*HollingA(i))
            ENDDO
         ENDIF         
         AlphaA(i)=(Comp-Mut)*(1+2*Delta*rnd) ! The rnd scaling factor is 2 here
         IF(Gamma0.lt.0)THEN    ! We model competition with Gammas=0 and Alphas increased by Mut
            AlphaA(i)=AlphaA(i)+Mut
         ENDIF
         AvA=AvA+AlphaA(i)
         VarA=VarA+AlphaA(i)**2
         WRITE(unit(tmp),*) AlphaA(i),Comp,Mut
      ENDDO
      CLOSE(unit(tmp))
      IF(Gamma0.lt.0)THEN       ! After modifying the alphas we declare that we model comptition.
         Gamma0=0
      ENDIF
      AvP=AvP/float(Sp)
      VarP=VarP/float(Sp)-AvP**2
      AvA=AvA/float(Sa)
      VarA=VarA/float(Sa)-AvA**2

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
*     ---Mod¡fied July 2012 (A.P-G.) I have modified the function ran2 to 
*     give a number between [-0.5-0.5]
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
      iy=iv(j)-idum2            ! Here idum is shued, idum and idum2 are combined to generate output
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
            
