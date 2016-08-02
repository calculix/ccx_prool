!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine beamextscheme(yil,ndim,nfield,lakonl,npropstart,prop,
     &  field,mi)
!
!     provides the extrapolation scheme for beams with a cross section
!     which is not rectangular nor elliptical
!
      implicit none
!
      character*8 lakonl
!
      integer npropstart,mi(*),ndim,nfield,j,k,l
!
      real*8 prop(*),ratio,ratio2,yil(ndim,mi(1)),yig(nfield,mi(1)),
     &  field(999,20*mi(3)),r,scal,a8(8,8),pmean1(nfield),pmean2(nfield)
!
!     extrapolation from a 2x2x2=8 integration point scheme in a hex to
!     the vertex nodes
!    
      data a8 /2.549,-.683,.183,-.683,-.683,.183,
     &        -.04904,.183,-.683,2.549,-.683,.183,
     &        .183,-.683,.183,-.04904,-.683,.183,
     &        -.683,2.549,.183,-.04904,.183,-.683,
     &        .183,-.683,2.549,-.683,-.04904,.183,
     &        -.683,.183,-.683,.183,-.04904,.183,
     &        2.549,-.683,.183,-.683,.183,-.683,
     &        .183,-.04904,-.683,2.549,-.683,.183,
     &        .183,-.04904,.183,-.683,-.683,.183,
     &        -.683,2.549,-.04904,.183,-.683,.183,
     &        .183,-.683,2.549,-.683/  
!
      if(lakonl(8:8).eq.'P') then
!
!        pipe cross section
!
!        the axis of the pipe is along the local xi-direction
!        the integration points are at positions +-0.57 along
!        the xi-axis. At each of these positions there are 8
!        integration points in the eta-zeta plane at one radial
!        position and
!        equally spaced along the circumferential direction
!
!        ratio of inner radius to outer radius
!
         ratio=(prop(npropstart+1)-prop(npropstart+2))/
     &              prop(npropstart+1)
         ratio2=ratio*ratio
!
!        radial location of integration points
!
         r=dsqrt((ratio2+1.d0)/2.d0)
!
!        scaling factor between the radial location of the integration
!        points and the regular location of C3D20R integration points
!
         scal=dsqrt(2.d0/3.d0)/r
!
!        calculating the mean values at each of the two xi-positions
!
         do k=1,nfield
            pmean1(k)=(yil(k,2)+yil(k,4)+yil(k,6)+yil(k,8))/4.d0
            pmean2(k)=(yil(k,10)+yil(k,12)+yil(k,14)+yil(k,16))/4.d0
         enddo
!
!        translating the results from the integration points of the
!        pipe section to the integration points of the C3D20R element
!
         do k=1,nfield
            yig(k,1)=pmean1(k)+(yil(k,6)-pmean1(k))*scal
            yig(k,2)=pmean2(k)+(yil(k,14)-pmean2(k))*scal
            yig(k,3)=pmean1(k)+(yil(k,8)-pmean1(k))*scal
            yig(k,4)=pmean2(k)+(yil(k,16)-pmean2(k))*scal
            yig(k,5)=pmean1(k)+(yil(k,4)-pmean1(k))*scal
            yig(k,6)=pmean2(k)+(yil(k,12)-pmean2(k))*scal
            yig(k,7)=pmean1(k)+(yil(k,2)-pmean1(k))*scal
            yig(k,8)=pmean2(k)+(yil(k,10)-pmean2(k))*scal
         enddo
!
!        standard extrapolation for the C3D20R element
!
         do j=1,8
            do k=1,nfield
               field(k,j)=0.d0
               do l=1,8
                  field(k,j)=field(k,j)+a8(j,l)*yig(k,l)
               enddo
            enddo
         enddo
!
      endif
!     
      return
      end
