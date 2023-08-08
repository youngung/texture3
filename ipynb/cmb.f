      program cmb
      implicit none
      character(len=256) fin, fout
      integer ngr,i,j
      parameter(ngr=20000)
      real*8 px(ngr,4)


!     fin='/Users/youngung/repo/mymtex_analysis/mooyeong/eddq.odf'
!     fin='/Users/youngung/repo/mtexscripts/steglich/Mg10Gd.odf'

      if (iargc().ne.2) then
         write(*,*)'requires 2 inputs (fin and fout)'
         stop -1
      endif

      call get_command_argument(1,fin)
      fin=trim(fin)
      call get_command_argument(2,fout)
      fout=trim(fout)

      call weighted_grains(ngr,fin,px)

      open(1,file=fout,status='unknown')
      write(1,*)
      write(1,*)
      write(1,*)
      write(1,*)'B',ngr
      do i=1,ngr
         write(1,'(4e13.5)')(px(i,j),j=1,4)
      enddo
      close(1)
      end program cmb


c     ---------------------------------------------------------------
      subroutine weighted_grains(ngrainran,fin,px)
      implicit none

      integer, intent(in)::ngrainran
      character(len=256),intent(in)::fin
      real*8, intent(out)::px(ngrainran,4)

      integer i,j,k,l,nphi,nphi1,nphi2,ix,iy,iz
      integer nphimx,nphi1mx,nphi2mx,ngrmx,jran,nr
      parameter(nphi1mx=72+1,nphimx=18+1,nphi2mx=18+1,ngrmx=1e6)
      real*8 phi,phi1,phi2,phi1_pre,phi1_next,res,val,
     $     cod(nphi1mx,nphimx,nphi2mx),cod_labo(nphi1mx,nphi2mx,nphimx),
     $     pi,omer(ngrmx),phir(ngrmx),ther(ngrmx),wgtr(ngrmx)
      real*8 r(3,2),phi1mx,phimx,phi2mx,c(0:1,0:1,0:1),aux3(3)
      integer ir(3,2)
      real random1,random2,random3
      logical verbose
      parameter(verbose=.false.)

      PI=4d0*dATAN(1d0)

      open(1,file=fin,status='unknown')


ccc   find the resolution
c     heading lines
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
c     res=phi1_next - phi1_pre
      do i=1,2
         read(1,*)phi1_pre,phi,phi2,val
         read(1,*)phi1_next,phi,phi2,val
         res=phi1_next-phi1_pre
      enddo

c     Find the maximum phi1, phi, phi2 value
      k=1
      do while(.true.)
         read(1,*,end=20)phi1mx,phimx,phi2mx,val
         k=k+1
         if (k.gt.100000)then
            write(*,*)'** line exceeded 100000'
            write(*,*)'** Too many lines... something went wrong?'
            stop
         endif
      enddo
 20   rewind 1

      if (verbose) then
         write(*,'(a9,f6.2)')'reso:',res
         write(*,'(a9,f6.2)')'phi1mx:',phi1mx
         write(*,'(a9,f6.2)')'phimx :',phimx
         write(*,'(a9,f6.2)')'phi2mx:',phi2mx
      endif

      nphi1 = int((phi1mx+res) / res) + 1
      nphi  = int(phimx / res) + 1
      nphi2 = int((phi2mx+res) / res) + 1

      if (verbose) then
         write(*,'(a9,i6)')'nphi1:',nphi1
         write(*,'(a9,i6)')'nphi :',nphi
         write(*,'(a9,i6)')'nphi2:',nphi2
      endif
      phi1mx=phi1mx+res
      phi2mx=phi2mx+res

      !! heading lines
      read(1,*)
      read(1,*)
      read(1,*)
      read(1,*)
      do i=1,nphi2-1
         do j=1,nphi
            do k=1,nphi1-1
               read(1,*,end=20)phi1,phi,phi2,val
               cod(k,j,i)=val
            enddo
         enddo
      enddo
      do j=1,nphi
         do k=1,nphi2
            cod(nphi1,j,k)=cod(1,j,k)
         enddo
      enddo
      do i=1,nphi1
         do j=1,nphi
            cod(i,j,nphi2)=cod(i,j,1)
         enddo
      enddo
c$$$      write(*,*)'cod(1,1,1):',cod(1,1,1)
c$$$      write(*,*)'cod(2,1,1):',cod(2,1,1)
c$$$      write(*,*)'cod(1,2,1):',cod(1,2,1)
c$$$      write(*,*)'cod(1,1,2):',cod(1,1,2)
c$$$      write(*,*)'cod(72,19,18):',cod(72,19,18)
c$$$      write(*,*)'cod(73,19,19):',cod(73,19,19)
      JRAN=-1                   ! defines seed for subsequent calls to random generator

      DO NR=1,ngrainran
        call random_number(random1)
        call random_number(random2)
        call random_number(random3)
        PHIR(NR)=360.*RANDOM1
        THER(NR)=ACOS(RANDOM2)*180./PI
        OMER(NR)= 90.*RANDOM3 ! cubic
c$$$        IF(ICRYSYM.EQ.1) OMER(NR)= 90.*RANDOM3 ! cubic
c$$$        IF(ICRYSYM.EQ.2) OMER(NR)= 60.*RANDOM3 ! hexagonal
c$$$        IF(ICRYSYM.EQ.5) OMER(NR)=180.*RANDOM3 ! orthogonal

        if ((phir(nr).eq.phi1mx).or.(phir(nr).eq.0.)) then
           write(*,*)'phir(nr) is on the edge'
           stop -1
        endif
        if ((ther(nr).eq.phimx).or.(ther(nr).eq.0.)) then
           write(*,*)'ther(nr) is on the edge'
           stop -1
        endif
        if ((omer(nr).eq.phi2mx).or.(omer(nr).eq.0.)) then
           write(*,*)'omer(nr) is on the edge'
           stop -1
        endif
      ENDDO

      do nr=1,ngrainran
         if (verbose) then
            write(*,'(a9,f6.2)')'phir(nr):',phir(nr)
            write(*,'(a9,f6.2)')'ther(nr):',ther(nr)
            write(*,'(a9,f6.2)')'omer(nr):',omer(nr)
         endif

         r(1,1)=phir(nr)-mod(phir(nr),res)
         r(1,2)=phir(nr)-mod(phir(nr),res)+res
         r(2,1)=ther(nr)-mod(ther(nr),res)
         r(2,2)=ther(nr)-mod(ther(nr),res)+res
         r(3,1)=omer(nr)-mod(omer(nr),res)
         r(3,2)=omer(nr)-mod(omer(nr),res)+res

         do i=1,3
            if (verbose) write(*,'(a9,2f6.1)')'r(i,1:2):',(r(i,j),j=1,2)
            do j=1,2
               ir(i,j)=int((r(i,j)+res)/res)
            enddo
         enddo

         do i=1,2
            do j=1,2
               do k=1,2
                  ix=ir(1,i)
                  iy=ir(2,j)
                  iz=ir(3,k)
c                  write(*,*)'(i,j,k)',ix,iy,iz
                  phi1=(ix-1)*res
                  phi =(iy-1)*res
                  phi2=(iz-1)*res
                  if (verbose) then
                     write(*,*)'(p1,p,p2,cod)',
     $                    phi1,phi,phi2,cod(ix,iy,iz)
                  endif
                  c(i-1,j-1,k-1)=cod(ix,iy,iz)
               enddo
            enddo
         enddo

         aux3(1)=phir(nr)
         aux3(2)=ther(nr)
         aux3(3)=omer(nr)
         call trilinear(c,aux3,r(1,1),r(1,2),r(2,1),r(2,2),r(3,1),r(3,2)
     $        ,val)
         if (verbose) write(*,'(a9,f9.5)')'val:',val
         px(nr,1)=phir(nr)
         px(nr,2)=ther(nr)
         px(nr,3)=omer(nr)
         px(nr,4)=val
      enddo

      close(1)
      return
      end subroutine

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine trilinear(cijk,xyz,x0,x1,y0,y1,z0,z1,c)
      implicit none
      real*8,intent(in)::cijk(0:1,0:1,0:1),x0,x1,y0,y1,z0,z1,
     $     xyz(3)
      real*8,intent(out)::c
c     locals
      real*8 x,y,z,xd,yd,zd,c00,c01,c10,c11,c0,c1
      x=xyz(1)
      y=xyz(2)
      z=xyz(3)

      xd=(x-x0)/(x1-x0)
      yd=(y-y0)/(y1-y0)
      zd=(z-z0)/(z1-z0)

      c00=cijk(0,0,0)*(1-xd)+cijk(1,0,0)*xd
      c01=cijk(0,0,1)*(1-xd)+cijk(1,0,1)*xd
      c10=cijk(0,1,0)*(1-xd)+cijk(1,1,0)*xd
      c11=cijk(0,1,1)*(1-xd)+cijk(1,1,1)*xd

      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      c=c0*(1-zd)+c1*zd
      return
      end subroutine trilinear
