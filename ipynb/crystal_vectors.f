c----------------------------------------------------------------------
c     This program calculates cartesian vectors of miller-indexed crystal
c     plane normals and crystal directions.
c     For compile, it requires below objects, and more ...
c         - voigt.o
c         - chg_basis.o
c         - crystal_symmetry.o
c         - fake_xit.o
      program main
c     arguments
c     ---------
c     single crystal file name (that contains the crystal symmetry) - the one used for VPSC or dEVPSC codes
c
c     usage
c     -----
c     compile procedure:
c         $ gfortran -c crystal_vectors.f -o crystal_vectors.o
c         $ gfortran crystal_vectors.o /Users/youngung/repo/evpsc/objdir/crystal_symmetry.o /Users/youngung/repo/evpsc/objdir/chg_basis.o /Users/youngung/repo/evpsc/objdir/voigt.o /Users/youngung/repo/evpsc/objdir/fake_xit.o /Users/youngung/repo/evpsc/objdir/data_crystal_fe.o /Users/youngung/repo/evpsc/objdir/funcs.o /Users/youngung/repo/evpsc/objdir/lu_libs.o /Users/youngung/repo/evpsc/objdir/thermo_emod_cte.o /Users/youngung/repo/evpsc/objdir/lib_write.o lib7.o -o schmid_factor
c     execution:
c         $ ./crystal_vectors ../../matData/vpscData/CP-Ti-ICN/dd_298.sx 1. 0. 0. 0. 0. 0.
      implicit none
c      real*8, intent(out)::R(3,3),P(3,3),Q(3,3)
c      integer,intent(in)::isn(4),isb(4),isnp(4),isbp(4)
      integer isn(4),isb(4),isnp(4),isbp(4)
      integer i,j,k,l,nphmx,ursx,ur8,iph,npoles,mxndeg,iphel,ntwmmx,
     $     nmodmx,nsysmx,ntwsmx
      parameter(nphmx=1,ursx=1,ur8=8,mxndeg=3,iphel=1,ntwmmx=1,
     $     nmodmx=4,nsysmx=48,ntwsmx=12)
      character*512 filecrys,prosa
      integer icrysym,iaux(4),nmodes(nphmx),nslmod(nphmx),ntwmod(nphmx),
     $     nsyst(nphmx),
     $     nslsys(nphmx),ntwsys(nphmx),nsm(nmodmx,nphmx),
     $     isense(nsysmx,nphmx),itwtype(nsysmx,nphmx),
     $     isectw(nsysmx,nphmx),itwinlaw,icrysymph(nphmx),
     $     ndegs_fit(nphmx,2),isx_ver

      real*8 sn(3),sb(3),st(3),snp(3),sbp(3),stp(3),sneq(3,24),aux3(3),
     $     celccv(6,6),cel(3,3,3,3),sel(3,3,3,3),aux6(6),aux33(3,3),
     $     aux3333(3,3,3,3),aux66(6,6),
     $     emod,celcc(6,6,nphmx),athcc(6,nphmx),selcc(6,6,nphmx),
     $     cijv_fit(nphmx,0:mxndeg,6,6),cte_fit(nphmx,0:mxndeg,6),
     $     dnca(3,nsysmx,nphmx),dbca(3,nsysmx,nphmx),fstar,
     $     schca(5,nsysmx,nphmx),nschca(6,nsysmx,nphmx),
     $     twsh(ntwmmx,nphmx),xmu_mode(nmodmx,nphmx),
     $     scauch_voigt(6),scauch(3,3),dum,energy
      logical verbose
c      parameter(verbose=.true.)
      parameter(verbose=.false.)
c     command line prompt
      integer nargs
      parameter(nargs=9)
c     argument types
c     1) single crystal file name
c     2)
      integer istatus(nargs)
      character*512 com_buf(nargs)
      integer isys,imod,ism,inverr,nind
      real*8 schf
c----------------------------------------------------------------------
c     command line arguments (sx file name, hkl indices)
      istatus(:)=0
      do i=1,iargc()
         istatus(i)=1
         call get_command_argument(i,com_buf(i))
      enddo

      NIND=3

      do i=1,iargc()
         if (i.eq.1) then
            filecrys=trim(com_buf(i))

            open (action='read',unit=ursx,file=filecrys,status='old')
            call sx_version_tester(ursx,isx_ver)
c      write(*,*)'isx_version:',isx_ver

            iph=1
            CALL CRYSTAL_SYMMETRY(1,URsx,icrysymph(iph),
     $           ISN,SN,SNEQ,ISB,SB,NPOLES)


            IF (ICRYSYMph(iph).EQ.2 .OR. ICRYSYMph(iph).EQ.3) NIND=4

c            write(*,*)'icystal:', icrysymph(iph)

         elseif (i.ge.2. .and. i.lt.2+nind) then ! 2, 3, 4
            read(com_buf(i),*)  isn(i-1)
         elseif (i.ge.2+nind .and. i.lt.2+nind*2) then ! 5, 6, 7
            read(com_buf(i),*) isb(i- (1+nind))
         else
            write(*,*)'** Only filecrystal is read'
            call xit()
         endif
      enddo

c      write(*,'(a4,x,4(i6.2))')'isn:',(isn(i),i=1,nind)
c      write(*,'(a4,x,4(i6.2))')'isb:',(isb(i),i=1,nind)
c      call xit()
      CALL CRYSTAL_SYMMETRY(2,URsx,icrysymph(iph),
     $     ISN,SN,SNEQ,ISB,SB,NPOLES)
      write(*,*) (sn(i),i=1,3)
      write(*,*) (sb(i),i=1,3)


      call crystal_symmetry(4,ursx,icrysymph(iph),
     $     ISN,SN,SNEQ,ISB,aux3,NPOLES)
c      write(*,*)'equivalent poles'
      write(*,*)'npoles:',npoles
      do i=1,npoles
         write(*,'(3f7.3)') (sneq(j,i),j=1,3)
      enddo


      call crystal_symmetry(4,ursx,icrysymph(iph),
     $     ISN,Sb,SNEQ,ISB,SB,NPOLES)
c      write(*,*)'equivalent poles'
      write(*,*)'npoles:',npoles
      do i=1,npoles
         write(*,'(3f7.3)') (sneq(j,i),j=1,3)
      enddo

      end program
