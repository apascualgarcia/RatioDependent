*     *********************************************
*     SUBROUTINE Integration
*     *********************************************                                                  
      
      SUBROUTINE Integration(Nfiles,unit,sizeA,sizeP,sizeT,BetaA,BetaP,
     &GammaA,GammaP,Na,Np,u,AlphaA,AlphaP,hhA,hhP,ggA,ggP)
      IMPLICIT NONE
      REAL*8 TOL,dt,TINY,Maxsteps,VERYTINY
      PARAMETER(TOL=1.0e-8,dt=10e-2,TINY=1.D-30,MaxSteps=1.e5,VERYTINY=1.D-90) ! Sometimes it is necessary to reduce dt, when the populations grow a lot
      INTEGER sizeA,sizeP,sizeT,Nfiles,unit(30)
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 Np(sizeP),Na(sizeA),u(sizeT)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 ggA(sizeA),ggP(sizeP)
*     This subroutine includes an integration loop calling to Runge Kutta
*     routine until a tolerance (TOL) has been achieved. We plan to consider
*     a fifth order Runge-Kutta version with adaptive stepside.
      INTEGER i,j,Nsteps,nok,nbad
      REAL*8 t,h,htry,hdid,hnext,hmin
      REAL*8 umax,uscal(sizeT),u1(sizeT),dydx(sizeT)
      REAL*8 error,errorTmp,test,eps,SteadyError
      REAL*8 Biomass

      !PRINT *, '>> Into integration routine..'  

      umax=0
      DO i=1,sizeT
         IF(u(i).gt.umax) umax=u(i)
      ENDDO
      t=0.01
      WRITE(unit(Nfiles-3),222) t,(u(i),i=1,sizeP)
      WRITE(unit(Nfiles-2),222) t,(u(sizeP+i),i=1,sizeA)
      hmin=TINY
      htry=dt
      eps=1.0e-10
      Biomass=0
      error=1.0d0               ! Init error just to enter into the loop      
      Nsteps=0
      SteadyError=0             ! It controls whether it is possible to improve the error
      errorTmp=0
      DO i=1,sizeP
         uscal(i)=Np(i)
      ENDDO
      DO i=1,sizeA
         uscal(sizeP+i)=Na(i)
      ENDDO
      DO WHILE(error.gt.TOL)
         error=-1.0d0            ! Set to zero once into the loop, in order to look for maximum error
         call Derivatives(u,sizeA,sizeP,sizeT,BetaA,BetaP,
     &     GammaA,GammaP,AlphaA,AlphaP,hhA,hhP,ggA,ggP,dydx)
         DO i=1,sizeT          
c$$$            IF(u(i).lt.VERYTINY)THEN ! I fix to zero here extincted species, comment this otherwise and...
c$$$               u(i)=0.0d0
c$$$               dydx(i)=0.0d0
c$$$               uscal(i)=0.0d0
c$$$            ENDIF
            u1(i)=u(i)
c            IF(u(i).gt.0)THEN ! this condition should be also cancelled, compute uscal for all of them (A.P-G)
               uscal(i)=ABS(u(i))+abs(dydx(i))+TINY
c            ENDIF
         ENDDO
         call bsstep(u,dydx,sizeT,t,htry,eps,uscal,
     &    sizeA,sizeP,BetaA,BetaP,GammaA,GammaP,
     &    AlphaA,AlphaP,hhA,hhP,ggA,ggP,hdid,hnext)

         Nsteps=Nsteps+1
         umax=0
         DO i=1,sizeT
c$$$            IF(MOD(Nsteps,100).eq.0) WRITE(*,*) u1(i),u(i),test,error,hnext,i,Nsteps ! debug
c$$$            IF(MOD(Nsteps-1,100).eq.0) WRITE(*,221) u1(i),u(i),test,error,hnext,i,Nsteps ! debug
            IF(i.le.sizeP)THEN
               test=ABS((u1(i)-u(i))/Np(i))
            ELSE
               test=ABS((u1(i)-u(i))/Na(i-sizeP))
            ENDIF
            u1(i)=u(i)
            IF(u(i).gt.umax) umax=u(i)
            IF(u(i).gt.TOL) Biomass=Biomass+u(i)
            IF(test.gt.error) error=test
         ENDDO
         IF(hdid.eq.h)then
            nok=nok+1
         else
            nbad=nbad+1
         ENDIF         
         IF(abs(hnext).lt.hmin)THEN
            WRITE(*,*) '*WARNING* stepsize smaller 
     & than minimum in odeint'
            call wait()
         ENDIF
         IF(MOD(Nsteps,50).eq.0)THEN
            PRINT *, ' -- ', Nsteps,error,errorTmp,ABS(errorTmp-error)
            WRITE(unit(Nfiles-3),222) t,(u(i),i=1,sizeP) ! Check if this print format works
            WRITE(unit(Nfiles-2),222) t,(u(sizeP+i),i=1,sizeA)
         ENDIF
         htry=hnext        
         IF(ABS(errorTmp-error).lt.eps)THEN
            SteadyError=SteadyError+1
         ELSE
            SteadyError=0
         ENDIF
         IF(SteadyError.gt.100)THEN            
            PRINT *, '* WARNING * Error improvement cannot be achieved'
            PRINT *, '* Error improvement: ',ABS(errorTmp-error)
            PRINT *, '* is lower than ',eps
            PRINT *, '* for more than 100 steps'
            PRINT *, '* Leaving the integration loop...'
            EXIT
         ENDIF
         errortmp=error
         IF(Nsteps.gt.MaxSteps)THEN
            PRINT *, '* WARNING * #steps: ',Nsteps,' larger than: ',MaxSteps
            PRINT *, '* Leaving the integration loop...'
            EXIT
         ENDIF
      ENDDO      
c      PRINT *,' -- Convergence achieved after ',Nsteps,' steps'
c      PRINT *,' -- Tolerance/Error: ',TOL,error
c      PRINT *, '<< Leaving Integration routine..'
 222  FORMAT(400(f12.4))
 221  FORMAT(2(f15.6),' ',2(f12.8),' ',f15.8,' ',i3,' ',i5)
      RETURN 
      END


*     *********************************************
*     SUBROUTINE bsstep
*     *********************************************                                                  
  
      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,
     &    sizeA,sizeP,BetaA,BetaP,GammaA,GammaP,
     &    AlphaA,AlphaP,hhA,hhP,ggA,ggP,hdid,hnext)
      IMPLICIT NONE
      INTEGER nv,NMAX,KMAXX,IMAX
      REAL*8 eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv)
      REAL*8 SAFE1,SAFE2,REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25,SAFE2=.7)
      PARAMETER(REDMAX=1.e-5,REDMIN=.7,TINY=1.e-30,SCALMX=.1) 
      INTEGER sizeA,sizeP
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 ggA(sizeA),ggP(sizeP)
C     USES Derivatives,mmid,pzextr
*     From Numerical Recipies chapters 16.4 and 16.3
*     Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust
*     stepsize. Input are the dependent variable vector y(1:nv) and its derivative dydx(1:nv)
*     at the starting value of the independent variable x. Also input are the stepsize to be attempted
*     htry, the required accuracy eps, and the vector yscal(1:nv) against which the
*     error is scaled. On output, y and x are replaced by their new values, hdid is the stepsize
*     that was actually accomplished, and hnext is the estimated next stepsize. derivs is the
*     user-supplied subroutine that computes the right-hand side derivatives. Be sure to set htry
*     on successive steps to the value of hnext returned from the previous step, as is the case
*     if the routine is called by odeint.
*     Parameters: NMAX is the maximum value of nv; KMAXX is the maximum row number used
*     in the extrapolation; IMAX is the next row number; SAFE1 and SAFE2 are safety factors;
*     REDMAX is the maximum factor used when a stepsize is reduced, REDMIN the minimum;
*     TINY prevents division by zero; 1/SCALMX is the maximum factor by which a stepsize can
*     be increased.
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      REAL*8 eps1,epsold,errmax,fact,h,red,scale,work,wrkmin
      REAL*8 values,h1
      REAL*8 xest,xnew,a(IMAX),alf(KMAXX,KMAXX),err(KMAXX)
      REAL*8 yerr(nv),ysav(nv),yseq(nv)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      DATA first/.true./,epsold/-1./
      DATA nseq /2,4,6,8,10,12,14,16,18/

      IF(eps.ne.epsold)THEN     ! A new tolerance, so reinitialize.
         hnext=-1.e29           !"Impossible" xnew.
         values=-1.e29
         eps1=SAFE1*eps
         a(1)=nseq(1)+1         ! Compute work coeficients Ak.
         DO k=1,KMAXX
            a(k+1)=a(k)+nseq(k+1)
         ENDDO
         DO iq=2,KMAXX          ! Compute Alpha(k q).
            DO  k=1,iq-1
               alf(k,iq)=eps1**((a(k+1)-a(iq+1))/
     &              ((a(iq+1)-a(1)+1.)*(2*k+1)))
            ENDDO 
         ENDDO 
         epsold=eps
         DO kopt=2,KMAXX-1      ! Determine optimal row number for convergence
            IF(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt)) GOTO 1 
         ENDDO 
 1       kmax=kopt
      ENDIF
      h=htry
      DO i=1,nv                 ! Save the starting values.
         ysav(i)=y(i)
      ENDDO 
      IF(h.ne.hnext.or.x.ne.xnew)THEN ! A new stepsize or a new integration: re-establish
         first=.true.           ! the order winDOw.
         kopt=kmax
      ENDIF
      reduct=.false.
 2    DO  k=1,kmax              ! Evaluate the sequence of modified midpoint
         xnew=x+h ! integrations.
         IF(xnew.eq.x)THEN
            WRITE(*,*) 'step size underflow in bsstep. You must choose a larger 
     &tolerance or set the minimum step size to a larger value.'
            call wait()
         ENDIF
         call mmid(ysav,dydx,nv,x,h,nseq(k),sizeA,sizeP,BetaA,BetaP,
     &        GammaA,GammaP,AlphaA,AlphaP,hhA,hhP,ggA,ggP,yseq)
         xest=(h/nseq(k))**2    ! Squared, since error series is even.
         call pzextr(k,xest,yseq,y,yerr,nv) ! Perform extrapolation.
         IF(k.ne.1)THEN         ! Compute normalized error estimate epsilon(k).
            errmax=TINY
            DO i=1,nv
               IF(yscal(i).gt.0)THEN ! APG, if I fix extincted species to zero ABOVE this condition must be imposed
                  errmax=max(errmax,abs(yerr(i)/yscal(i))) ! Otherwise you will get here NANs
               ENDIF
            ENDDO 
            errmax=errmax/eps   ! Scale error relative to tolerance.
            km=k-1
            err(km)=(errmax/SAFE1)**(1./(2*km+1))
         ENDIF
c         PRINT *, x,xnew,h,errmax,' x,xnew,h,errmax DEBUG' ! Check this if you find underflow above 
         IF(k.ne.1.and.(k.ge.kopt-1.or.first))THEN ! In order window.
            IF(errmax.lt.1.)GOTO 4 ! Converged.
            IF(k.eq.kmax.or.k.eq.kopt+1)THEN ! Check for possible stepsize reduction.
               red=SAFE2/err(km)
               GOTO 3
            ELSE IF(k.eq.kopt)THEN
               IF(alf(kopt-1,kopt).lt.err(km))THEN
                  red=1./err(km)
                  GOTO 3
               ENDIF
            ELSE IF(kopt.eq.kmax)THEN
               IF(alf(km,kmax-1).lt.err(km))THEN
                  red=alf(km,kmax-1)*SAFE2/err(km)
                  GOTO 3
               ENDIF
            ELSE IF(alf(km,kopt).lt.err(km))THEN
               red=alf(km,kopt-1)/err(km)
               GOTO 3
            ENDIF
         ENDIF
      ENDDO 
 3    red=min(red,REDMIN)       ! Reduce stepsize by at least REDMIN and at
      red=max(red,REDMAX)       ! most REDMAX.
      h=h*red
      reduct=.true.
      GOTO 2                    ! Try again.
 4    x=xnew                    ! Successful step taken.
      hdid=h1 
      first=.false.
      wrkmin=1.e35              ! Compute optimal row for convergence and
      DO kk=1,km                ! corresponding stepsize.
         fact=max(err(kk),SCALMX)
         work=fact*a(kk+1)
         IF(work.lt.wrkmin)THEN
            scale=fact
            wrkmin=work
            kopt=kk+1
         ENDIF
      ENDDO 
      hnext=h/scale
      IF(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)THEN ! Check for possible order increase, but not IF stepsize
         fact=max(scale/alf(kopt-1,kopt),SCALMX) ! was just reduced.
         IF(a(kopt+1)*fact.le.wrkmin)THEN
            hnext=h/fact
            kopt=kopt+1
         ENDIF
      ENDIF
      
      RETURN
      END

*     *********************************************
*     SUBROUTINE pzextr
*     *********************************************          
      
      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      IMPLICIT NONE
      INTEGER iest,nv,IMAX,NMAX
      REAL*8 xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=25,NMAX=500)
*     Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a
*     sequence of estimates with progressively smaller values x = xest, and corresponding function
*     vectors yest(1:nv). This call is number iest in the sequence of calls. Extrapolated
*     function values are output as yz(1:nv), and their estimated error is output as dy(1:nv).
*     Parameters: Maximum expected value of iest is IMAX; of nv is NMAX.
      INTEGER j,k1
      REAL*8 delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x

      x(iest)=xest              ! Save current independent variable.
      DO j=1,nv
         dy(j)=yest(j)
         yz(j)=yest(j)
      ENDDO
      IF(iest.eq.1)THEN        ! Store First estimate in first column.
         DO j=1,nv
            qcol(j,1)=yest(j)
         ENDDO
      ELSE
         DO j=1,nv
            d(j)=yest(j)
         ENDDO 
         DO  k1=1,iest-1
            delta=1./(x(iest-k1)-xest)
            f1=xest*delta
            f2=x(iest-k1)*delta
            DO j=1,nv ! Propagate tableau 1 diagonal more.
               q=qcol(j,k1)
               qcol(j,k1)=dy(j)
               delta=d(j)-q
               dy(j)=f1*delta
               d(j)=f2*delta
               yz(j)=yz(j)+dy(j)
            ENDDO
         ENDDO 
         DO j=1,nv
            qcol(j,iest)=dy(j)
         ENDDO
      ENDIF
      RETURN
      END

*     *********************************************
*     SUBROUTINE mmid
*     *********************************************      
      
      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,sizeA,sizeP,BetaA,BetaP,
     &        GammaA,GammaP,AlphaA,AlphaP,hhA,hhP,ggA,ggP,yout)
      IMPLICIT NONE
      INTEGER nstep,nvar,NMAX
      INTEGER sizeA,sizeP
      REAL*8 htot,xs,dydx(nvar),y(nvar),yout(nvar)
      PARAMETER (NMAX=50)
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 ggA(sizeA),ggP(sizeP)
*     Modified midpoint step. Dependent variable vector y(1:nvar) and its derivative vector
*     dydx(1:nvar) are input at xs. Also input is htot, the total step to be made, and nstep,
*     the number of substeps to be used. The output is returned as yout(1:nvar), which need
*     not be a distinct array from y; if it is distinct, however, then y and dydx are returned
*     undamaged.
      INTEGER i,n
      REAL*8 h,h2,swap,x,ym(nvar),yn(nvar)
      h=htot/nstep              ! Stepsize this trip.
      DO i=1,nvar
         ym(i)=y(i)
         yn(i)=y(i)+h*dydx(i)   ! First step.
      ENDDO 
      x=xs+h
      call Derivatives(yn,sizeA,sizeP,nvar,BetaA,BetaP,
     &     GammaA,GammaP,AlphaA,AlphaP,hhA,hhP,ggA,ggP,yout)
      h2=2.*h
      DO n=2,nstep ! General step.
         DO i=1,nvar
            swap=ym(i)+h2*yout(i)
            ym(i)=yn(i)
            yn(i)=swap
         ENDDO          
         x=x+h
         call Derivatives(yn,sizeA,sizeP,nvar,BetaA,BetaP,
     &     GammaA,GammaP,AlphaA,AlphaP,hhA,hhP,ggA,ggP,yout)
      ENDDO
      
      DO i=1,nvar ! Last step.                  
         yout(i)=0.5*(ym(i)+yn(i)+h*yout(i))
c         PRINT *,i,ym(i),yn(i),yout(i),' !debug'
      ENDDO 
      RETURN

      END

*     *********************************************
*     SUBROUTINE Derivatives
*     *********************************************                                                  
      
      SUBROUTINE Derivatives(u,sizeA,sizeP,sizeT,BetaA,BetaP,
     &     GammaA,GammaP,AlphaA,AlphaP,hhA,hhP,ggA,ggP,f)
      IMPLICIT NONE
      INTEGER sizeA,sizeP,sizeT,idum
      REAL*8 BetaA(sizeA,sizeA),BetaP(sizeP,sizeP)
      REAL*8 GammaA(sizeA,sizeP),GammaP(sizeP,sizeA)
      REAL*8 AlphaA(sizeA),AlphaP(sizeP)
      REAL*8 u(sizeT)
      REAL*8 hhA(sizeA),hhP(sizeP)
      REAL*8 ggA(sizeA),ggP(sizeP)
*     Evaluates the ODEs for given values of the abundances at time t
      INTEGER i,j,k
      REAL*8 HollingA(sizeA),HollingP(sizeP)
      REAL*8 GollingA(sizeA),GollingP(sizeP)
      REAL*8 f(sizeT),Comp,Int
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

      DO i=1,Sp                 ! Compute first a Holling like term to saturate the ODEs
         HollingP(i)=0.0d0
         GollingP(i)=0.0d0 ! New Holling term, a hybrid between Holling and Gollum
         IF(hhP(i).gt.0)THEN
            DO k=1,Sa
               HollingP(i)=HollingP(i)+ABS(GammaP(i,k))*u(Sp+k)
            ENDDO
         ENDIF
         IF(ggP(i).gt.0)THEN
            DO k=1,Sa
               GollingP(i)=GollingP(i)+ABS(GammaA(k,i))*u(Sp+k)
            ENDDO
         ENDIF
      ENDDO
      DO i=1,Sa
         HollingA(i)=0.0d0
         GollingA(i)=0.0d0
         IF(hhA(i).gt.0)THEN
            DO k=1,Sp
               HollingA(i)=HollingA(i)+ABS(GammaA(i,k))*u(k)
            ENDDO
         ENDIF
         IF(ggA(i).gt.0)THEN
            DO k=1,Sp
               GollingA(i)=GollingA(i)+ABS(GammaP(k,i))*u(k)
            ENDDO
         ENDIF
      ENDDO
*     --- Compute the derivatives
      DO i=1,Sp                 ! Compute first Holling term to saturate the ODEs
         Comp=0.0d0
         Int=0.0d0
         DO j=1,Sp
            Comp=Comp+BetaP(i,j)*u(j)
         ENDDO
         DO k=1,Sa
            Int=Int+GammaP(i,k)*u(Sp+k)/(f0P+hhP(i)*HollingP(i)+ggA(k)*GollingA(k))
         ENDDO         
         f(i)=u(i)*(AlphaP(i)-Comp+Int)
      ENDDO
      DO i=1,Sa                 ! Compute first Holling term to saturate the ODEs
         Comp=0.0d0
         Int=0.0d0
         DO j=1,Sa
            Comp=Comp+BetaA(i,j)*u(Sp+j)
         ENDDO         
         DO k=1,Sp
            Int=Int+GammaA(i,k)*u(k)/(f0A+hhA(i)*HollingA(i)+ggP(k)*GollingP(k))
         ENDDO         
         f(Sp+i)=u(Sp+i)*(AlphaA(i)-Comp+Int)
      ENDDO

      RETURN
      END

*     *********************************************
*     SUBROUTINE Wait
*     *********************************************                                                  
*     Subroutine to substitute the obsolete function PAUSE
      
      Subroutine wait()
      implicit none
      character(1) Key

      write(*,*) 'press any key to continue'
      READ(*,*)
      !Key = GETCHARQQ()

      end 
