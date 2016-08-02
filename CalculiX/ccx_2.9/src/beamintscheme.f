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
      subroutine beamintscheme(lakonl,mint3d,npropstart,prop,
     &  kk,xi,et,ze,weight)
!
!     provides the integration scheme for beams with a cross section
!     which is not rectangular nor elliptical
!
!     mint3d: number of integration points (returned if kk=0)
!     xi,et,ze: local coordinates of integration point kk
!     weight: weight for integration point kk
!
      implicit none
!
      character*8 lakonl
!
      integer mint3d,npropstart,jj,kk
!
      real*8 prop(*),xi,et,ze,weight,ratio,ratio2,dtheta,theta,r
!
      intent(in) lakonl,npropstart,prop,
     &  kk
!
      intent(inout) weight,xi,et,ze,mint3d
!
      if(lakonl(8:8).eq.'P') then
!
!        pipe cross section
!
         if(kk.eq.0) then
            mint3d=16
            return
         endif
!
!        ratio of inner radius to outer radius
!
         ratio=(prop(npropstart+1)-prop(npropstart+2))/
     &              prop(npropstart+1)
         ratio2=ratio*ratio
!
         if(kk.gt.8) then
            jj=kk-8
            xi=1.d0/dsqrt(3.d0)
         else
            jj=kk
            xi=-1.d0/dsqrt(3.d0)
         endif
!
!        pi/4
!
         dtheta=datan(1.d0)
!
         theta=(jj-1)*dtheta
         r=dsqrt((ratio2+1.d0)/2.d0)
!
         et=r*dcos(theta)
         ze=r*dsin(theta)
         weight=dtheta*(1.d0-ratio2)/2.d0
      endif
!     
      return
      end
