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

c     distance
      xd=(x-x0)/(x1-x0)
      yd=(y-y0)/(y1-y0)
      zd=(z-z0)/(z1-z0)

c     Apply linear interpolation on the top and bottom faces
c     to obtain two separate interpolants on the respective faces.
      c00=cijk(0,0,0)*(1-xd)+cijk(1,0,0)*xd
      c01=cijk(0,0,1)*(1-xd)+cijk(1,0,1)*xd
      c10=cijk(0,1,0)*(1-xd)+cijk(1,1,0)*xd
      c11=cijk(0,1,1)*(1-xd)+cijk(1,1,1)*xd

c     linear interpolation using the two interpolants
c     on the top and bottom faces
      c0=c00*(1-yd)+c10*yd
      c1=c01*(1-yd)+c11*yd

      c=c0*(1-zd)+c1*zd
      return
      end subroutine trilinear
