c----------------------------------------------------------------------
c     This program calculates schmid factor of given slip system
c     For compile, it requires below objects, and more ...
c         - voigt.o
c         - chg_basis.o
c         - crystal_symmetry.o
c         - fake_xit.o
      program main
c     arguments
c     ---------
c     single crystal file name (that contains the crystal symmetry, elastic modulus)
c
c     usage
c     -----
c     compile procedure:
c         $ gfortran -c schmid_factor.f  -o schmid_factor.o
c         $ gfortran schmid_factor.o /Users/youngung/repo/evpsc/objdir/crystal_symmetry.o /Users/youngung/repo/evpsc/objdir/chg_basis.o /Users/youngung/repo/evpsc/objdir/voigt.o /Users/youngung/repo/evpsc/objdir/fake_xit.o /Users/youngung/repo/evpsc/objdir/data_crystal_fe.o /Users/youngung/repo/evpsc/objdir/funcs.o /Users/youngung/repo/evpsc/objdir/lu_libs.o /Users/youngung/repo/evpsc/objdir/thermo_emod_cte.o /Users/youngung/repo/evpsc/objdir/lib_write.o lib7.o -o schmid_factor
c     execution:
c         $ ./schmid_factor ../../matData/vpscData/CP-Ti-ICN/dd_298.sx 1. 0. 0. 0. 0. 0.
      implicit none
c      real*8, intent(out)::R(3,3),P(3,3),Q(3,3)
c      integer,intent(in)::isn(4),isb(4),isnp(4),isbp(4)
      integer isn(4),isb(4),isnp(4),isbp(4)
      integer i,j,k,l,nphmx,ur1,ur8,iph,npoles,mxndeg,iphel,ntwmmx,
     $     nmodmx,nsysmx,ntwsmx
      parameter(nphmx=1,ur1=1,ur8=8,mxndeg=3,iphel=1,ntwmmx=1,
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
      parameter(nargs=5)
      integer istatus(nargs)
      character*512 com_buf(nargs)
      integer isys,imod,ism,inverr
      real*8 schf
c----------------------------------------------------------------------
c     command line arguments (sx file name, hkl indices)
      istatus(:)=0
      do i=1,iargc()
         istatus(i)=1
         call get_command_argument(i,com_buf(i))
      enddo
      do i=1,iargc()
         if (i.eq.1) then
            filecrys=trim(com_buf(i))
         elseif (i.gt.1 .and. i.lt.8) then
            read(com_buf(i),*) scauch_voigt(i-1)
         else
            write(*,*)'** Only filecrystal is read'
            call xit()
         endif
      enddo
c----------------------------------------------------------------------
c     Convert stress in Voigt notation to that in 3x3 components
      call voigt(scauch_voigt,scauch,aux66,aux3333,1)
      call w_chrc(0,'Stress tensor in crystal axes')
      call w_mdim(0,scauch,3,1d0)
c----------------------------------------------------------------------
c     elastic modulus
      OPEN(action='read',UNIT=UR1,FILE=FILECRYS,STATUS='OLD')
      do i=1,4
         READ(UR1,'(A)') PROSA
      enddo
      READ(UR1,*) ((CELCCV(I,J),J=1,6),i=1,6) ! C: stiffness, S: compliance
      rewind ur1
      call voigt(aux6,aux33,celccv,cel,3)

c----------------------------------------------------------------------
c     Read slip plane normal and slip directions
c     in crystal axes (Cartesian)
      isx_ver=0
c      OPEN(action='read',UNIT=UR1,FILE=FILECRYS,STATUS='OLD')
      call data_crystal_fe(     ! phase specific
c             in
     $     ur1,ur8,isx_ver,mxndeg,iphel,ntwmmx,nmodmx,
     $     nsysmx,ntwsmx,nphmx,
c             out
     $     celcc,athcc,selcc,nmodes,nslmod,ntwmod,nsyst,nslsys,
     $     ntwsys,twsh,nsm,dnca,dbca,schca,nschca,isense,itwtype,
     $     xmu_mode,icrysymph,fstar,ndegs_fit,cijv_fit,cte_fit)
      close(ur1)
c----------------------------------------------------------------------
c     elastic compliance
      aux66(:,:)=selcc(:,:,1)
      call chg_basis(aux6,aux33,aux66,sel,3,6)
      call w_chrc(0,'sel 66')
      call w_mdim(0,aux66,6,1d0)
c----------------------------------------------------------------------
c     elastic energy
      call calc_el_energy(sel,scauch,energy)
      call w_val(0,'Elastic energy:',energy)
c----------------------------------------------------------------------
c     schmid factor for each slip system
      ism=0
      iph=1
      do imod=1, nmodes(iph)
         write(*,*)'** imod:',imod
         do isys=1,nsm(imod,iph)
            ism=ism+1
c            write(*,'(a6,x,i6,x,a6,x,i6)')'isys:',isys,'ism:',ism
            call calc_schfactor(dnca(:,ism,iph),dbca(:,ism,iph),scauch,
     $           schf)
            call proj_young_emod(dnca(:,ism,iph),cel,emod)
            write(*,
     $  '(a6,x,3f6.2,x,
     $    a6,x,3f6.2,x,
     $    a6,x,1f8.4,x,
     $    a6,x,1e10.2)')
     $           'dnca:',dnca(:,ism,iph),
     $           'dbca:',dbca(:,ism,iph),
     $           'schf:',schf,
     $           'emod:',emod
         enddo
         write(*,*)
      enddo
      end program




      subroutine calc_schfactor(sn,sb,stress,factor)
      real*8, intent(in):: sn(3),sb(3),stress(3,3)
      real*8, intent(out):: factor

c     locals
      real*8 s6(6),aux66(6,6),aux3333(3,3,3,3),vnorm,xmag
      integer i,j

      call chg_basis(s6,stress,aux66,aux3333,2,6)
      xmag=vnorm(s6,6)

      factor=0.
      do i=1,3
         do j=1,3
            factor=factor+sn(i)*sb(j)*stress(i,j)
         enddo
      enddo
      factor=factor/xmag
      end subroutine


      subroutine proj_young_emod(sn,cel,emod)
      real*8, intent(in):: sn(3), cel(3,3,3,3)
      real*8, intent(out):: emod
c     locals
      integer i,j,k,l

c     Calculate the projected elastic modulus
      emod=0d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         emod=emod+cel(i,j,k,l)*sn(i)*sn(j)*sn(k)*sn(l)
      enddo
      enddo
      enddo
      enddo
      end subroutine proj_young_emod

c     ----
      subroutine calc_el_energy(sel,stress,energy)
      real*8, intent(in):: sel(3,3,3,3),stress(3,3)
      real*8, intent(out):: energy
c     locals
      integer i,j,k,l
      energy=0.
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         energy=energy + sel(i,j,k,l)*stress(i,j)*stress(k,l)
      enddo
      enddo
      enddo
      enddo

      end subroutine calc_el_energy
