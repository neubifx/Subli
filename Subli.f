      Program Sublimation
      
      implicit none

      real*8 :: pi, plan, bol, uma, avo, gas, gaskcal, bohr, luz,
     *zpe1
      
      integer :: dim, i, j, k, m, nm
      parameter (dim=200)
      

      real*8 :: ene1, mol1, nsim1, mom1(3), prodI1, frot1,
     *vk1(dim), deom1(dim), freq1(dim), Eexc1(dim), prodQ1(dim),
     *h1(dim), s1(dim), gib1(dim), qrimp1(dim,dim), qrlivre1(dim,dim),
     *qoh1(dim,dim), qoh21(dim,dim), qmo1(dim,dim), b(dim),
     *Ggas(dim,dim), Hgas(dim)
      real*8 :: qv, tvh, tvs, xix, qrimpt, xe, zpem, bvo, expbvo,
     *expxe, emor, termo1, sumtermo1, sumtermo2, qmorse, hmorse, smorse
      real*8 :: qtrans, htrans, strans, qrot, hrot, srot, qe, somae,
     *hele, sele,  qvib, stvh, stvs, hvib, svib,  hrotit, srotat, ired
      real*8 :: enes, pt, conf, Elat, Ecoe
      real*8 :: temp(dim), pres(dim), Gsol(dim,dim), Ssol(dim,dim),
     *Hsol(dim,dim), DeltaG(dim,dim), DeltaH(dim,dim)

      real*8 :: Pref, Tsubref, DeltaHsub, Aanto, Banto
      
      real*8 :: x1(dim), y1(dim), sX, sY, sXX, sXY, slope, clin,
     *Tsub(dim)
      
      real*8 :: BESSI
      
      integer :: lin1, nfreq1, nrotor1, irotor1(200),
     *sroti1(200), groti1(200), nom1, iom1(200), nest1, g1(200)
      
      integer :: ntemp, npres, pos(dim), posref
      
      character*30 fsol,fgas,fout
      character :: bla

      pi=3.14159265
      plan=6.626068D-34
      bol=1.3806503D-23
      uma=1.660538782D-27
      avo=6.0221415D+23
      gas=0.00008314472
      gaskcal=1.987D-3
      bohr=5.292D-11
      luz=2.998D+10
38    format (A65)
40    format (' ',125('*'))
41    format (A50)
42    format (A8,200(ES14.4))
43    format (F8.2,200(F14.4))
44    format (F10.3, ES17.3, ES20.3)
45    format (A44, A37, ES10.4, A7)
46    format (A7, F8.3, A16, F8.3, A7)
47    format (A9, F8.3, A7)
48    format (A43, A9, A9)

      write(*,*) 'Arquivo de entrada - solido ?'
      read(*,'(a)') fsol
      write(*,*) 'Arquivo de entrada - gas ?'
      read(*,'(a)') fgas
      write(*,*) 'Arquivo de sa¡da ?'
      read(*,'(a)') fout

      open(5,file=fgas,status='old')
      open(6,file=fout)
C     citation
      write(6,*)
      write(6,*) '                   Sublimation Program - 2021'
      write(6,*) '            Laboratorio de Cinetica Quimica - UFRRJ'
      write(6,*) '                      Version: May, 2021'
      write(6,*)
      write(6,*) 'Please cite: '
      write(6,38) ' ',
     *'*************************************************************',
     *'  Xavier, N.F., Bauerfeldt, G.F. Determination of the        ',
     *'  Cohesive Properties and Sublimation Temperatures of        ',
     *'  Glycine Polymorphs. Crystal Growth & Design, v. 21,        ',
     *'  p. 6266-6275, 2021. doi.org/10.1021/acs.cgd.1c00724        ',
     *'*************************************************************',
     *' '
      write(6,*)
      write(6,48) 'Calculations are based on results from the ', fsol,
     *' equation'
      write(6,*)
C
      read(5,*)
      read(5,*) ene1, enes, pt, conf
      Elat=enes-ene1
      Ecoe=Elat+pt+conf
      write(6,*) 'Lattice Energy is', Elat, 'kJ/mol'
      write(6,*)
      write(6,*) 'Cohesion Energy is', Ecoe, 'kJ/mol'
      write(6,*)
      read(5,*) mol1
      read(5,*) lin1, nsim1, mom1(1), mom1(2), mom1(3)
      if(lin1.ne.-1) then
        prodI1=mom1(1)*mom1(2)*mom1(3)
        frot1=(((pi*prodI1*(uma*bohr**2)**3)**0.5)/nsim1)
      endif
      read(5,*) nfreq1, nrotor1, (b(i), i=1,4*nrotor1)
      j=0
      do i=1,nrotor1
        j=j+1
        irotor1(i)=b(j)
        j=j+1
        vk1(i)=b(j)
        j=j+1
        sroti1(i)=b(j)
        j=j+1
        groti1(i)=b(j)
      enddo

      read(5,*) nom1, (b(i), i=1,2*nom1)
      j=0
      do i=1,nom1
        j=j+1
        iom1(i)=b(j)
        j=j+1
        deom1(i)=b(j)
      enddo

      if (nfreq1.eq.-1) then
        read(5,*)
      else
        read(5,*) (freq1(i), i=1,nfreq1)
        do i=1,nfreq1
           zpe1=zpe1+freq1(i)
        enddo
        zpe1=zpe1/699.52
      endif
      read(5,*) nest1
      read(5,*) (g1(i), i=1, nest1)
      read(5,*) (Eexc1(i), i=1, nest1)
      read(5,*) ntemp, npres, posref
      
      close(5)
      open(5,file=fsol,status='old')
      do i=1,(13+2*ntemp)
      read(5,*)
      enddo
      read(5,*) bla, (pres(j), j=1,npres)
      do i=1,ntemp
      read(5,*) temp(i), (Gsol(i,j), j=1,npres)
      enddo
      do j=1,npres
        do i=1,ntemp
           if (i.eq.1) then
              Ssol(i,j)=-(Gsol(i+1,j)-Gsol(i,j))/(temp(i+1)-temp(i))
           elseif (i.eq.ntemp) then
              Ssol(i,j)=-(Gsol(i,j)-Gsol(i-1,j))/(temp(i)-temp(i-1))
           else
              Ssol(i,j)=-(Gsol(i+1,j)-Gsol(i-1,j))/(temp(i+1)-temp(i-1))
           endif
           Hsol(i,j)=Gsol(i,j)+temp(i)*Ssol(i,j)
        enddo
      enddo
C      Termodinƒmica do g s

      do i=1,ntemp
C
        qtrans=((2*pi*mol1*uma*bol*temp(i)/plan**2)**1.5)
     **(temp(i)*gas/avo)
        htrans=2.5*gaskcal*temp(i)
        strans=(gaskcal*1000)*(log(qtrans)+2.5)
C
        if (lin1.eq.1) then
           qrot=frot1*(8*pi**2*bol*temp(i)/plan**2)**1.5
           hrot=1.5*gaskcal*temp(i)
           srot=(gaskcal*1000)*(log(qrot)+1.5)
        elseif (lin1.eq.0) then
           qrot=8*pi**2*mom1(3)*uma*bohr**2*bol*temp(i)/plan**2
           hrot=gaskcal*temp(i)
           srot=(gaskcal*1000)*(log(qrot)+1.0)
        else
           qrot=1
           hrot=0
           srot=0
        endif
C
        qe=0
        somae=0
        do j=1,nest1
           qe=qe+g1(j)*exp(-Eexc1(j)/(gaskcal*temp(i)))
           somae=somae+g1(j)*(-Eexc1(j)/(gaskcal*temp(i)))*
     *exp(-Eexc1(j)/(gaskcal*temp(i)))
        enddo
        hele=gaskcal*temp(i)*somae/qe
        sele=gaskcal*1000*((log(qe))+somae/qe)

        prodQ1(i)=qrot*qtrans*qe
        h1(i)=htrans+hrot+hele+(enes-Ecoe)/4.184
        s1(i)=strans+srot+sele

        qvib=1
        hvib=0
        svib=0
        if (nfreq1.ne.-1) then
           stvh=0
           stvs=0
           do j=1,nfreq1
              qvib=qvib/(1-exp(-plan*freq1(j)*luz/
     *(bol*temp(i))))
              stvh=stvh+plan*freq1(j)*luz/(bol*
     *(exp(plan*freq1(j)*luz/bol/temp(i))-1))
              stvs=stvs+(plan*freq1(j)*luz/(bol*temp(i)*
     *(exp(plan*freq1(j)*luz/bol/temp(i))-1)))-log(1-
     *exp(-plan*luz*freq1(j)/bol/temp(i)))
           enddo
           hvib=gaskcal*stvh+zpe1
           svib=stvs*gaskcal*1000
        endif
        if (nrotor1.gt.0) then
           qrimpt=1
           hrotit=0
           srotat=0
           do k=1,nrotor1
              qv=1/(1-exp(-plan*freq1(irotor1(k))*luz/
     *(bol*temp(i))))
              tvh=plan*freq1(irotor1(k))*luz/
     *(bol*(exp(plan*freq1(irotor1(k))*luz/bol/temp(i))-1))
              tvs=(plan*freq1(irotor1(k))*luz/(bol*temp(i)*
     *(exp(plan*freq1(irotor1(k))*luz/bol/temp(i))-1)))-
     *log(1-exp(-plan*luz*freq1(irotor1(k))/bol/temp(i)))

              qvib=qvib/qv
              stvh=stvh-tvh
              stvs=stvs-tvs

              ired=(vk1(k)*4184/avo)*(((sroti1(k))**2)/8/
     *(pi**2))/((freq1(irotor1(k))*luz)**2)
              qrlivre1(i,k)=groti1(k)*((8*(pi**3)*bol*temp(i)*
     *ired/((plan*sroti1(k))**2))**0.5)

              qoh1(i,k)=qv

              xix=vk1(k)/2/(gaskcal*temp(i))
              qrimp1(i,k)=qrlivre1(i,k)*(exp(-xix))*BESSI(0,xix)
              qrimpt=qrimpt*qrimp1(i,k)
              hrotit=hrotit+1.5*gaskcal*temp(i)+
     *(vk1(k)/2)*(1-BESSI(1,xix)/BESSI(0,xix))
              srotat=srotat+gaskcal*1000*(log(qrimp1(i,k)))+
     *gaskcal*1000/2+(vk1(k)*1000/2/temp(i))*(1-BESSI(1,xix)/
     *BESSI(0,xix))
           enddo
           hvib=gaskcal*stvh+zpe1
           svib=stvs*gaskcal*1000
           prodQ1(i)=prodQ1(i)*qvib*qrimpt
           h1(i)=h1(i)+hvib+hrotit
           s1(i)=s1(i)+svib+srotat
        else
           prodQ1(i)=prodQ1(i)*qvib
           h1(i)=h1(i)+hvib
           s1(i)=s1(i)+svib
        endif
C
        if (nom1.gt.0) then

           do k=1,nom1

              qv=1/(1-exp(-plan*freq1(iom1(k))*luz/
     *(bol*temp(i))))
              tvh=freq1(iom1(k))/699.52+gaskcal*plan*
     *freq1(iom1(k))*luz/(bol*(exp(plan*freq1(iom1(k))*
     *luz/bol/temp(i))-1))
              tvs=gaskcal*1000*(plan*freq1(iom1(k))*luz/
     *(bol*temp(i)*(exp(plan*freq1(iom1(k))*luz/bol/temp(i))-1)))-
     *log(1-exp(-plan*luz*freq1(iom1(k))/bol/temp(i)))
C
              if (deom1(k).gt.0) then
                 xe=freq1(iom1(k))/(4*deom1(k))
              else
                 xe=-deom1(k)
              endif
              zpem=freq1(iom1(k))/2-freq1(iom1(k))*xe/4
              bvo=freq1(iom1(k))/(0.695*temp(i))
              nm=int((1-xe)/(2*xe))
              expbvo=exp(-bvo*(1-xe))
              expxe=exp(bvo*xe)
              sumtermo1=1
              sumtermo2=1
              do m=1,nm
      emor=freq1(iom1(k))*(m+0.5)-freq1(iom1(k))*xe*(m+0.5)**2
      termo1=(expbvo**m)*(expxe**(m**2))
      sumtermo1=sumtermo1+termo1
      sumtermo2=sumtermo2+termo1*emor
              enddo
              qmorse=sumtermo1
              qoh21(i,k)=qv
              qmo1(i,k)=qmorse
              hmorse=((1/qmorse)*sumtermo2+zpem)/349.76
              smorse=gaskcal*1000*(log(qmorse))+(1/
     *(qmorse*0.695*temp(i)))*(sumtermo2)
C
              prodQ1(i)=prodQ1(i)*qmorse/qv
              h1(i)=h1(i)-tvh+hmorse
              s1(i)=s1(i)-tvs+smorse
           enddo
        endif
C
        gib1(i)=h1(i)-temp(i)*s1(i)/1000
        do j=1,npres
        Ggas(i,j)=4.184*(gib1(i)+gaskcal*temp(i)*LOG(pres(j)*1000))
        enddo
        Hgas(i)=4.184*h1(i)
C      write(6,*) h1(i), s1(i), gib1(i)
      enddo

      write(6,40)
      write(6,*) 'Gibss Free Energy - Solid'
      write(6,*)
      write(6,41) 'P(kbar)'
      write(6,42) 'T(K)', (pres(i),i=1,npres)
      do i=1,ntemp
      write(6,43) temp(i), (Gsol(i,j), j=1,npres)
      enddo
      write(6,40)
      write(6,*) 'Gibss Free Energy - Gas'
      write(6,*)
      write(6,41) 'P(kbar)'
      write(6,42) 'T(K)', (pres(i),i=1,npres)
      do i=1,ntemp
      write(6,43) temp(i), (Ggas(i,j), j=1,npres)
      enddo
      
      do i=1,ntemp
         do j=1,npres
            DeltaH(i,j)=Hgas(i)-Hsol(i,j)
            DeltaG(i,j)=Ggas(i,j)-Gsol(i,j)
         enddo
      enddo
      write(6,40)
      write(6,*) 'Delta Gibss Free Energy'
      write(6,*)
      write(6,41) 'P(kbar)'
      write(6,42) 'T(K)', (pres(i),i=1,npres)
      do i=1,ntemp
      write(6,43) temp(i), (DeltaG(i,j), j=1,npres)
      enddo
      do j=1,npres
         do i=1,ntemp
            if (i.gt.1) then
               if ((DeltaG(i,j)/DeltaG(i-1,j)).lt.0) then
                   pos(j)=i
                   goto 307
               endif
            endif
         enddo
307      if (pos(j).gt.(ntemp-3)) then
            pos(j)=ntemp-3
         endif
         continue
      enddo
      write(6,40)
      write(6,*) 'Sublimation Properties'
      write(6,*)
      write(6,*) 'Tsub (K)      Pressure (kbar)     Pressure (Pa)'
C     proceeds linear regression of DeltaG
      do j=1,npres
      do i=1,5
         x1(i)=temp(i+pos(j)-3)
         y1(i)=DeltaG(i+pos(j)-3,j)
      enddo
        sX=0
        sY=0
        sXX=0
        sXY=0
        do i=1,5
         sX=sX+x1(i)
         sy=sy+y1(i)
         sXX=sXX+x1(i)*x1(i)
         sXY=sXY+x1(i)*y1(i)
        enddo
        slope=((sX*sY)-(5*sXY))/((sX*sX)-(5*sXX))
        clin=(sY-((((sX*sY)-(5*sXY))/((sX*sX)-(5*sXX)))*sX))/5
        Tsub(j)=-clin/slope
        write(6,44) Tsub(j), pres(j), pres(j)*1E8
      enddo
C
      do i=1,npres
         x1(i)=1/Tsub(i)
         y1(i)=LOG(pres(i))
      enddo
        sX=0
        sY=0
        sXX=0
        sXY=0
        do i=1,npres
         sX=sX+x1(i)
         sy=sy+y1(i)
         sXX=sXX+x1(i)*x1(i)
         sXY=sXY+x1(i)*y1(i)
        enddo
        slope=((sX*sY)-(npres*sXY))/((sX*sX)-(npres*sXX))
        DeltaHsub=-8.314E-3*slope
      write(6,*)
      write(6,47) 'DHsub = ', DeltaHsub, ' kJ/mol'
      write(6,*)
      write(6,40)
      write(6,*)
C      posref=7
      do j=1,npres
         if (j.eq.posref) then
            Tsubref=Tsub(j)
            Pref=pres(j)
            goto 360
         endif
      enddo
360   continue
      do j=1,npres
         do i=1,ntemp
            if (pres(j).eq.Pref) then
               if (temp(i).gt.Tsubref) then
                  DeltaHsub=DeltaH(i,j)-(DeltaH(i,j)-DeltaH(i-1,j))*
     *(temp(i)-Tsubref)/(temp(i)-temp(i-1))
               goto 364
               endif
            endif
         enddo
      enddo
364   continue
      write(6,45) ' Sublimation Properties will be recalculated',
     *' assuming the reference pressure as ', Pref, ' kbar :'
      write(6,*)
      write(6,46) 'Tsub =', Tsubref, ' K and DHsub = ', DeltaHsub,
     *' kJ/mol'
      Aanto=0.434294*DeltaHsub/(4.184*gaskcal)
      Banto=0.434294*(LOG(Pref)+DeltaHsub/(4.184*gaskcal*Tsubref))
      write(6,*)
      write(6,*) 'Parameters of the Antoine Equation',
     *' (logP(kbar) = -A/T + B):'
      write(6,*)
      write(6,*) 'A = ', Aanto
      write(6,*) 'B = ', Banto
      write(6,*)
      write(6,*) 'Final (T,P) values:'
      write(6,*) 'Tsub (K)      Pressure (kbar)     Pressure (Pa)'
      do j=1,npres
        Tsub(j)=-Aanto/(0.434294*LOG(pres(j))-Banto)
        write(6,44) Tsub(j), pres(j), pres(j)*1E8
      enddo
      
      End Program
      ! ----------------------------------------------------------------------
      FUNCTION BESSI(N,X)
!
!     This subroutine calculates the first kind modified Bessel function
!     of integer order N, for any REAL X. We use here the classical
!     recursion formula, when X > N. For X < N, the Miller's algorithm
!     is used to avoid overflows.
!     REFERENCE:
!     C.W.CLENSHAW, CHEBYSHEV SERIES FOR MATHEMATICAL FUNCTIONS,
!     MATHEMATICAL TABLES, VOL.5, 1962.
!
      IMPLICIT NONE
      INTEGER, PARAMETER :: IACC = 40
      REAL*8, PARAMETER :: BIGNO = 1.D10, BIGNI = 1.D-10
      INTEGER N, M, J
      REAL *8 X,BESSI,BESSI0,BESSI1,TOX,BIM,BI,BIP
      IF (N.EQ.0) THEN
      BESSI = BESSI0(X)
      RETURN
      ENDIF
      IF (N.EQ.1) THEN
      BESSI = BESSI1(X)
      RETURN
      ENDIF
      IF(X.EQ.0.D0) THEN
      BESSI=0.D0
      RETURN
      ENDIF
      TOX = 2.D0/X
      BIP = 0.D0
      BI  = 1.D0
      BESSI = 0.D0
      M = 2*((N+INT(SQRT(FLOAT(IACC*N)))))
      DO 12 J = M,1,-1
      BIM = BIP+DFLOAT(J)*TOX*BI
      BIP = BI
      BI  = BIM
      IF (ABS(BI).GT.BIGNO) THEN
      BI  = BI*BIGNI
      BIP = BIP*BIGNI
      BESSI = BESSI*BIGNI
      ENDIF
      IF (J.EQ.N) BESSI = BIP
12    CONTINUE
      BESSI = BESSI*BESSI0(X)/BI
      RETURN
      END
! ----------------------------------------------------------------------
! Auxiliary Bessel functions for N=0, N=1
      FUNCTION BESSI0(X)
      IMPLICIT NONE
      REAL *8 X,BESSI0,Y,P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,
     *Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/1.D0,3.5156229D0,3.0899424D0,
     *1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI0=AX*BX
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
      FUNCTION BESSI1(X)
      IMPLICIT NONE
      REAL *8 X,BESSI1,Y,P1,P2,P3,P4,P5,P6,P7,
     *Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9,AX,BX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *-0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *-0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF(ABS(X).LT.3.75D0) THEN
      Y=(X/3.75D0)**2
      BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
      AX=ABS(X)
      Y=3.75D0/AX
      BX=EXP(AX)/SQRT(AX)
      AX=Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))
      BESSI1=AX*BX
      ENDIF
      RETURN
      END
! ----------------------------------------------------------------------
