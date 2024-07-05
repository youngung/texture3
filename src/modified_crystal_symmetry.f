c
c ***********************************************************************
C     SUBROUTINE CRYSTAL_SYMMETRY   --->   version 09/JAN/2009
c
c *** If IOPTION=1:
c     Reads crystal symmetry 'icrysym' and unit cell parameters.
c     Generates vectors 'cvec(i,n)' of the unit cell.
c     Generates symmetry operators 'h(i,j,nsymop)' for all crystal symmetries.
c *** If IOPTION=2:
c     Reads Miller indices of systems in 3 or 4-index notation 'isn(i)'
c     & 'isb(i)'. Calculates normal & burgers vectors 'sn(i)' & 'sb(i)'
c *** If IOPTION=3:
c     Generates 'nequiv' crystallographically equivalent orientations sneq(i,n)
c     of normal vector sn(i) by applying all the symmetry operations to it.
c     Discards repeated orientations and defines 'nequiv'.
c *** If IOPTION=4:
c     Generates 'nequiv' crystallographically equivalent orientations sneq(i,n)
c     of normal vector sn(i) by applying all the symmetry operations to it.
c     Discards repeated orientations and defines 'nequiv'.
c *** Simmetry parameter ICRYSYM:
c        1: CUBIC
c        2: HEXAGONAL
c        3: TRIGONAL
c        4: TETRAGONAL
c        5: ORTHORHOMBIC
c        6: MONOCLINIC
c        7: TRICLINIC
c ***********************************************************************

      SUBROUTINE CRYSTAL_SYMMETRY(
c          in
     $     ioption,ur1,
c          inout
     $     icrysym,isn,
c          out
     $     sn,sneq,

c          inout
     $     isb,
c          out
     $     sb,nequiv)
      implicit none

      integer, intent(in):: ioption,ur1
      integer, intent(out):: nequiv
      integer, intent(inout):: icrysym,isn(4),isb(4)
      real*8,  intent(out):: sn(3),sneq(3,24),sb(3)

      real*8 h(3,3,24),hx(3,3,6),cdim(3),cang(3),cvec(3,3)
      integer itag(24)

      character crysym*5
      save h,nsymop,cvec

c     locals
      integer i,j,k,m,n,isign,nind,nr,mn,nrot,nsymop
      real*8 sndif,pi,eta,chi,snnor,sbnor,ang
      data pi /3.1415926535898/

c ****************************************************************************

      if(ioption.eq.1) then

        read(ur1,*)
        read(ur1,'(a)') crysym
        icrysym=0
        if(crysym.eq.'cubic' .or. crysym.eq.'CUBIC') icrysym=1
        if(crysym.eq.'hexag' .or. crysym.eq.'HEXAG') icrysym=2
        if(crysym.eq.'trigo' .or. crysym.eq.'TRIGO') icrysym=3
        if(crysym.eq.'tetra' .or. crysym.eq.'TETRA') icrysym=4
        if(crysym.eq.'ortho' .or. crysym.eq.'ORTHO') icrysym=5
        if(crysym.eq.'monoc' .or. crysym.eq.'MONOC') icrysym=6
        if(crysym.eq.'tricl' .or. crysym.eq.'TRICL') icrysym=7
        if(icrysym.eq.0) then
          write(*,*) ' *** CANNOT RECOGNIZE THE CRYSTAL SYMMETRY'
          call xit()
        endif

        READ(UR1,*) (CDIM(i),i=1,3),(CANG(i),i=1,3)

        DO I=1,3
          CANG(I)=CANG(I)*PI/180d0
        ENDDO

c$$$        write(*,*)'pi:',pi
c$$$        write(*,'(3e20.13)')(cdim(i),i=1,3),(cang(i),i=1,3)

c *** assumes 'c' coincident with 'z' and 'a' in the plane 'xz'
        CVEC(1,1)=SIN(CANG(2))
        CVEC(2,1)=0d0
        CVEC(3,1)=COS(CANG(2))
        CVEC(1,2)=(COS(CANG(3))-COS(CANG(1))*COS(CANG(2)))/SIN(CANG(2))
        CVEC(3,2)=COS(CANG(1))
        CVEC(2,2)=SQRT(1d0-CVEC(1,2)**2-CVEC(3,2)**2)
        CVEC(1,3)=0d0
        CVEC(2,3)=0d0
        CVEC(3,3)=1d0

        DO J=1,3
        DO I=1,3
          CVEC(I,J)=CDIM(J)*CVEC(I,J)
        ENDDO
        ENDDO

        DO I=1,3
        DO J=1,3
          DO M=1,6
            HX(I,J,M)=0.d0
          ENDDO
          DO N=1,24
            H(I,J,N)=0.d0
          ENDDO
        ENDDO
        ENDDO

c *** identity operation ---> triclinic & all symmetries
        do i=1,3
           h(i,i,1)=1.d0
        enddo
        nsymop=1

c *** 180 deg rotation around (001) ---> orthorhombic, monoclinic
        if(icrysym.eq.5 .or. icrysym.eq.6) then
           h(1,1,2)= cos(pi)
           h(2,2,2)= cos(pi)
           h(3,3,2)= 1.d0
           h(1,2,2)=-sin(pi)
           h(2,1,2)= sin(pi)
           nsymop=2
        endif

c *** x-mirror & y-mirror ---> orthorhombic
        if(icrysym.eq.5) then
           h(1,1,3)=-1.d0
           h(2,2,3)= 1.d0
           h(3,3,3)= 1.d0

           h(1,1,4)= 1.d0
           h(2,2,4)=-1.d0
           h(3,3,4)= 1.d0
           nsymop=4
        endif

c *** cubic symmetry
        if(icrysym.eq.1) then

c *** rotations of (pi/3) & (2*pi/3) around <111>
           hx(1,3,1)= 1.d0
           hx(2,1,1)= 1.d0
           hx(3,2,1)= 1.d0

           hx(1,2,2)= 1.d0
           hx(2,3,2)= 1.d0
           hx(3,1,2)= 1.d0

           do m=1,2
              do n=1,nsymop
                 mn=m*nsymop+n
                 do i=1,3
                    do j=1,3
                       do k=1,3
                          h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           nsymop=mn

c *** mirror across the plane (110)
           hx(1,2,3)= 1.d0
           hx(2,1,3)= 1.d0
           hx(3,3,3)= 1.d0

           do n=1,nsymop
              mn=nsymop+n
              do i=1,3
                 do j=1,3
                    do k=1,3
                       h(i,j,mn)=h(i,j,mn)+hx(i,k,3)*h(k,j,n)
                    enddo
                 enddo
              enddo
           enddo
           nsymop=mn
c *** rotations of 90, 180, 270 around x3

           do m=1,3
              ang=pi/2.*float(m)
              hx(1,1,m)= cos(ang)
              hx(2,2,m)= cos(ang)
              hx(3,3,m)= 1.0
              hx(1,2,m)=-sin(ang)
              hx(2,1,m)= sin(ang)
              hx(1,3,m)= 0.0
              hx(3,1,m)= 0.0
              hx(2,3,m)= 0.0
              hx(3,2,m)= 0.0
           enddo

           do m=1,3
              do n=1,nsymop
                 mn=m*nsymop+n
                 do i=1,3
                    do j=1,3
                       do k=1,3
                          h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           nsymop=mn
        endif                   !end of condition for icrysym=1

c *** hexagonal, trigonal and tetragonal symmetry

        if(icrysym.ge.2 .and. icrysym.le.4) then
           if(icrysym.eq.2) nrot=6
           if(icrysym.eq.3) nrot=3
           if(icrysym.eq.4) nrot=4

c *** mirror plane at 30 deg or 60 deg or 45 deg with respect to x1
           ang=pi/float(nrot)
           h(1,1,2)= cos(ang)**2-sin(ang)**2
           h(2,2,2)=-h(1,1,2)
           h(3,3,2)= 1.d0
           h(1,2,2)= 2.*cos(ang)*sin(ang)
           h(2,1,2)= h(1,2,2)
           nsymop=2

c *** rotations of 2*pi/6 around axis <001> for hexagonals.
c *** rotations of 2*pi/3 around axis <001> for trigonals.
c *** rotations of 2*pi/4 around axis <001> for tetragonals.
           do nr=1,nrot-1
              ang=nr*2.*pi/nrot
              hx(1,1,nr)= cos(ang)
              hx(2,2,nr)= cos(ang)
              hx(3,3,nr)= 1.d0
              hx(1,2,nr)=-sin(ang)
              hx(2,1,nr)= sin(ang)
           enddo

           do m=1,nrot-1
              do n=1,nsymop
                 mn=m*nsymop+n
                 do i=1,3
                    do j=1,3
                       do k=1,3
                          h(i,j,mn)=h(i,j,mn)+hx(i,k,m)*h(k,j,n)
                       enddo
                    enddo
                 enddo
              enddo
           enddo
           nsymop=mn

        endif                   !end of condition for icrysym= 2,3,4

c     write(10,*)
      write(*,'(''  # of symmetry operations='',i4)') nsymop
      write(*,'(''  symmetry matrices'')')
      write(*,'(i3,9f7.3)') (n,((h(i,j,n),j=1,3),i=1,3),n=1,nsymop)

      endif                     !end of condition for ioption=1

c **************************************************************************
c   Converts Miller-Bravais indices of plane normal and slip direction
c   into normalized vectors sn(i) and sb(i), respectively.
c   Indices for cubic (1), tetragonal (4), orthorhombic (5), monoclinic (6)
c   & triclinic (7) systems are in 3-index notation.
c   For hexagonal (2) & trigonal (3) systems uses 4-index notation.
c **************************************************************************

      if (ioption.eq.2) then

        if(icrysym.eq.2 .or. icrysym.eq.3) then
          isn(3)=isn(4)
          isb(1)=isb(1)-isb(3)
          isb(2)=isb(2)-isb(3)
          isb(3)=isb(4)
        endif

c *** assumes 'c' coincident with 'z' and 'a' in the plane 'xz'
        sn(3)= isn(3)/cvec(3,3)
        sn(1)=(isn(1)-cvec(3,1)*sn(3))/cvec(1,1)
        sn(2)=(isn(2)-cvec(1,2)*sn(1)-cvec(3,2)*sn(3))/cvec(2,2)

c$$$        write(*,*)'cvec'
c$$$        do i=1,3
c$$$           write(*,'(3e20.13)')(cvec(i,j),j=1,3)
c$$$        enddo
c$$$
c$$$        write(*,*)'sn'
c$$$        write(*,'(3e20.13)')(sn(j),j=1,3)

        snnor=sqrt(sn(1)**2+sn(2)**2+sn(3)**2)
        do j=1,3
          sn(j)=sn(j)/snnor
          if(abs(sn(j)).lt.1.e-03) sn(j)=0d0
        enddo

c$$$        write(*,*)'sn'
c$$$        write(*,'(3e20.13)')(sn(j),j=1,3)

c *** this block specific for EPSC & VPSC

        do i=1,3
          sb(i)=isb(1)*cvec(i,1)+isb(2)*cvec(i,2)+isb(3)*cvec(i,3)
        enddo
        sbnor=sqrt(sb(1)**2+sb(2)**2+sb(3)**2)
        do j=1,3
          sb(j)=sb(j)/sbnor
          if(abs(sb(j)).lt.1.e-03) sb(j)=0d0
        enddo

      endif      ! end of if(ioption.eq.2)

c **************************************************************************
c *** generates all symmetry related vectors sneq(i,n) with z>0.
c *** eliminates redundant poles: coincidents and opposites
c **************************************************************************

      IF(IOPTION.EQ.3) THEN
c         write(*,*)'crystal symmetry for ioption.eq.3'
c         write(*,*) 'icrysym:',icrysym
c         write(*,*) 'nsymop:',nsymop
c         call w_chrc(0,'crystal symmetry for ioption.eq.3')
c         call w_ival(0,'icrysym:',icrysym)
         NIND=3
         IF (ICRYSYM.EQ.2 .OR. ICRYSYM.EQ.3) NIND=4
         READ(UR1,*) (ISN(I),I=1,NIND),CHI,ETA
!          IF(NIND.EQ.3) WRITE(10,'(3I4,2F10.1)') (ISN(I),I=1,3),CHI,ETA
!          IF(NIND.EQ.4) WRITE(10,'(4I4,2F10.1)') (ISN(I),I=1,4),CHI,ETA
         ETA=ETA*PI/180.0
         CHI=CHI*PI/180.0
         SB(1)=COS(ETA)*SIN(CHI)
         SB(2)=SIN(ETA)*SIN(CHI)
         SB(3)=          COS(CHI)

        if(nind.eq.4) isn(3)=isn(4)
        SN(1)= ISN(1)/CVEC(1,1)
        SN(2)=(ISN(2)-CVEC(1,2)*SN(1))/CVEC(2,2)
        SN(3)=(ISN(3)-CVEC(1,3)*SN(1)-CVEC(2,3)*SN(2))/CVEC(3,3)
        SNNOR=SQRT(SN(1)**2+SN(2)**2+SN(3)**2)
        DO J=1,3
          SN(J)=SN(J)/SNNOR
          IF(ABS(SN(J)).LT.1.D-03) SN(J)=0.D0
        ENDDO

        DO N=1,NSYMOP
          ITAG(N)=0
          DO I=1,3
          SNEQ(I,N)=0.D0
            DO J=1,3
              SNEQ(I,N)=SNEQ(I,N)+H(I,J,N)*SN(J)
            ENDDO
          ENDDO
        ENDDO

        IF(ICRYSYM.NE.7) THEN      ! NSYMOP=1 FOR TRIGONAL
          DO M=1,NSYMOP-1
            IF(ITAG(M).EQ.0) THEN
              DO N=M+1,NSYMOP
                SNDIF=ABS(SNEQ(1,M)-SNEQ(1,N))+ABS(SNEQ(2,M)-SNEQ(2,N))
     #               +ABS(SNEQ(3,M)-SNEQ(3,N))
                IF(SNDIF .LE. 0.0001) ITAG(N)=1
                SNDIF=ABS(SNEQ(1,M)+SNEQ(1,N))+ABS(SNEQ(2,M)+SNEQ(2,N))
     #               +ABS(SNEQ(3,M)+SNEQ(3,N))
                IF(SNDIF .LE. 0.0001) ITAG(N)=1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

        NEQUIV=0
        DO N=1,NSYMOP
          IF(ITAG(N).EQ.0) THEN
            NEQUIV=NEQUIV+1
            ISIGN=1
            IF(SNEQ(3,N).LT.0.) ISIGN=-1
            SNEQ(1,NEQUIV)=ISIGN*SNEQ(1,N)
            SNEQ(2,NEQUIV)=ISIGN*SNEQ(2,N)
            SNEQ(3,NEQUIV)=ISIGN*SNEQ(3,N)
          ENDIF
        ENDDO

c        call xit()
      ENDIF            !END OF IF(IOPTION=3)


c **************************************************************************
c *** generates all symmetry related vectors sneq(i,n) with z>0.
c *** eliminates redundant poles: coincidents and opposites
c **************************************************************************

      IF(IOPTION.EQ.4) THEN
c         call w_chrc(0,'crystal symmetry for ioption.eq.3')
c         call w_ival(0,'icrysym:',icrysym)
         NIND=3
         IF (ICRYSYM.EQ.2 .OR. ICRYSYM.EQ.3) NIND=4
c         READ(UR1,*) (ISN(I),I=1,NIND),CHI,ETA
!          IF(NIND.EQ.3) WRITE(10,'(3I4,2F10.1)') (ISN(I),I=1,3),CHI,ETA
!          IF(NIND.EQ.4) WRITE(10,'(4I4,2F10.1)') (ISN(I),I=1,4),CHI,ETA
c         ETA=ETA*PI/180.0
c         CHI=CHI*PI/180.0
c         SB(1)=COS(ETA)*SIN(CHI)
c         SB(2)=SIN(ETA)*SIN(CHI)
c         SB(3)=          COS(CHI)

c        if(nind.eq.4) isn(3)=isn(4)
c        SN(1)= ISN(1)/CVEC(1,1)
c        SN(2)=(ISN(2)-CVEC(1,2)*SN(1))/CVEC(2,2)
c        SN(3)=(ISN(3)-CVEC(1,3)*SN(1)-CVEC(2,3)*SN(2))/CVEC(3,3)
        SNNOR=SQRT(SN(1)**2+SN(2)**2+SN(3)**2)
        DO J=1,3
          SN(J)=SN(J)/SNNOR
          IF(ABS(SN(J)).LT.1.D-03) SN(J)=0.D0
        ENDDO

        DO N=1,NSYMOP
          ITAG(N)=0
          DO I=1,3
          SNEQ(I,N)=0.D0
            DO J=1,3
              SNEQ(I,N)=SNEQ(I,N)+H(I,J,N)*SN(J)
            ENDDO
          ENDDO
        ENDDO

        IF(ICRYSYM.NE.7) THEN      ! NSYMOP=1 FOR TRIGONAL
          DO M=1,NSYMOP-1
            IF(ITAG(M).EQ.0) THEN
              DO N=M+1,NSYMOP
                SNDIF=ABS(SNEQ(1,M)-SNEQ(1,N))+ABS(SNEQ(2,M)-SNEQ(2,N))
     #               +ABS(SNEQ(3,M)-SNEQ(3,N))
                IF(SNDIF .LE. 0.0001) ITAG(N)=1
                SNDIF=ABS(SNEQ(1,M)+SNEQ(1,N))+ABS(SNEQ(2,M)+SNEQ(2,N))
     #               +ABS(SNEQ(3,M)+SNEQ(3,N))
                IF(SNDIF .LE. 0.0001) ITAG(N)=1
              ENDDO
            ENDIF
          ENDDO
        ENDIF

        NEQUIV=0
        DO N=1,NSYMOP
          IF(ITAG(N).EQ.0) THEN
            NEQUIV=NEQUIV+1
            ISIGN=1
            IF(SNEQ(3,N).LT.0.) ISIGN=-1
            SNEQ(1,NEQUIV)=ISIGN*SNEQ(1,N)
            SNEQ(2,NEQUIV)=ISIGN*SNEQ(2,N)
            SNEQ(3,NEQUIV)=ISIGN*SNEQ(3,N)
          ENDIF
        ENDDO

c        call xit()
      ENDIF            !END OF IF(IOPTION=4)

      return
      end subroutine crystal_symmetry
