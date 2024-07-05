c *****************************************************************************
c
      subroutine gauleg(x1,x2,x,w,n)
      implicit none
c     implicit real*8(a-h,p-z)
      integer n,xj,j,iter,xi,i,m,xn
      real*8 x(n),w(n),x1,x2,pi,eps,p1,p2,p3,pp,
     $     z1,z,xl,xm

c      parameter(eps=3.d-14)
c
c     changed by R.L. -  8/2/97
c
      parameter(eps=1.e-07)
      pi=4.d0*datan(1.d0)
      m=(n+1)/2
      xm=0.5d0*(x1+x2)
      xl=0.5d0*(x2-x1)
      xn=n
      do 12 i=1,m
         xi=i
         z=dcos(pi*(xi-.25d0)/(xn+0.5d0))
c
         iter=0
 1       continue
         iter=iter+1
c
c     R.L. 8/2/97
c
         if(iter.gt.10000) then
            write(*,*)
     $           'GAULEG WARNING: TOL 1.e-07 NEVER REACHED - ERR = ',
     $           abs(z-z1)
            return
         endif
c
         p1=1.d0
         p2=0.d0
         do 11 j=1,n
            xj=j
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(xj-1.d0)*p3)/xj
 11      continue
         pp=n*(z*p1-p2)/(z*z-1.d0)
         z1=z
         z=z1-p1/pp

         if(abs(z-z1).gt.eps) go to 1
         x(i)=xm-xl*z
         x(n+1-i)=xm+xl*z
         w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
         w(n+1-i)=w(i)
 12   continue
      return
      end subroutine gauleg
C
C *****************************************************************************
C
      FUNCTION ran2(idum)
      implicit none
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL*8 ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     $     IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     $     NTAB=32,NDIV=1+int(real(IMM1)/real(NTAB)),
     $     EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software #>,13.
c
C *****************************************************************************
c
      SUBROUTINE ludcmp(a,n,np,indx,d,isingular)
      implicit none
      INTEGER n,np,indx(n),NMAX
      REAL*8 d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k,isingular
      REAL*8 aamax,dum,sum,vv(NMAX)
      d=1d0
      do 12 i=1,n
        aamax=0d0
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
c
c        if (aamax.eq.0.) pause 'singular matrix in ludcmp'
c
        if(aamax.eq.0d0) then
           isingular=1
        return
        endif
c
        vv(i)=1d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0d0

        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax

        if(a(j,j).eq.0d0) a(j,j)=TINY

c       if(a(j,j).eq.0.) then
c       isingular=1
c       return
c       endif

        if(j.ne.n)then
          dum=1d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
c
      isingular=0
c
      return
      END
c
c *****************************************************************************
c
      SUBROUTINE lubksb(a,n,np,indx,b)
      implicit none
      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL*8 sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
c
c *****************************************************************************
c
      SUBROUTINE jacobi(a,n,np,d,v,nrot,ier)
      implicit none
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j,ier
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0d0
11      continue
        v(ip,ip)=1d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+dabs(a(ip,iq))

14        continue
15      continue
c
        if(sm.eq.0.)then
        ier=0
        return
        endif
c
        if(i.lt.4)then
           tresh=0.2d0*sm/n**2
        else
          tresh=0d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
             g=100d0*dabs(a(ip,iq))
            if((i.gt.4).and.(dabs(d(ip))+
     *g.eq.dabs(d(ip))).and.(dabs(d(iq))+g.eq.dabs(d(iq))))then
              a(ip,iq)=0d0
            else if(dabs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(dabs(h)+g.eq.dabs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5d0*h/a(ip,iq)
                t=1./(dabs(theta)+dsqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1d0/dsqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                if (dabs(g)<1e-20) g=0d0
                if (dabs(h)<1e-20) h=0d0
                if (dabs(s)<1e-20) s=0d0
                if (dabs(tau)<1e-20) tau=0d0

c$$$                write(*,*)'g:',g
c$$$                write(*,*)'s:',s
c$$$                write(*,*)'h:',h
c$$$                write(*,*)'tau:',tau
c$$$
c$$$                write(*,*)'g*tau',g*tau
c$$$                write(*,*)'h+g*tau',h+g*tau
c$$$                write(*,*)'s*(h+g*tau)',s*(h+g*tau)
c$$$                write(*,*)'g-s*(h+g*tau)',g-s*(h+g*tau)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)

                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0d0
23      continue
24    continue
c      pause 'too many iterations in jacobi'
c
      ier=1
c
      return
      END
c
c *****************************************************************************
c
      SUBROUTINE eigsrt(d,v,n,np)
      implicit none
      INTEGER n,np
      REAL*8 d(np),v(np,np)
      INTEGER i,j,k
      REAL*8 p
      do i=1,n-1
         k=i
         p=d(i)
         do j=i+1,n
            if(d(j).ge.p)then
               k=j
               p=d(j)
            endif
         enddo
         if(k.ne.i)then
            d(k)=d(i)
            d(i)=p
            do j=1,n
               p=v(j,i)
               v(j,i)=v(j,k)
               v(j,k)=p
            enddo
         endif
      enddo
      return
      END

c-----------------------------------------------------------------------
c     Rotate the given 2nd order tensor <a33> using the given
c     transformation matrix (rot) and save to <b33>.

c     The transformation matrix rotates the (new<-old) wise
c     In the explicit index-notation, the opertaion is expressed:
c         **  b_ij = r_ik a_kl r_jl
      subroutine rot_2nd_tensor(a33,rot,b33)
      implicit none
      dimension a33(3,3),rot(3,3),b33(3,3)
      real*8 a33,rot,b33
      integer i,j,k,l
      b33(:,:) = 0d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
         b33(i,j) = b33(i,j) + rot(i,k) * a33(k,l) * rot(j,l)
      enddo
      enddo
      enddo
      enddo
      return
      end subroutine rot_2nd_tensor
c-----------------------------------------------------------------------c
c     Rotate the given 4th order tensor <a3333> using the given
c     transformation matrix (rot) and save to <b3333>.

c     The transformation matrix rotates the (new<-old) wise
c     In the explicit index-notation, the opertaion is expressed:
c         **  b_ijkl =  rot_im rot_jn rot_ko rot_lp a_mnop
      subroutine rot_4th_tensor(a3333,rot,b3333)
      implicit none
      dimension a3333(3,3,3,3),rot(3,3),b3333(3,3,3,3)
      real*8 a3333,rot,b3333
      integer i,j,k,l,m,n,o,p
      b3333(:,:,:,:) = 0d0
      do i=1,3
      do j=1,3
      do k=1,3
      do l=1,3
c     --
      do m=1,3
      do n=1,3
      do o=1,3
      do p=1,3
         b3333(i,j,k,l) = b3333(i,j,k,l) +
     $        rot(i,m)  * rot(j,n) * rot(k,o) * rot(l,p)
     $        * a3333(m,n,o,p)
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end subroutine rot_4th_tensor
