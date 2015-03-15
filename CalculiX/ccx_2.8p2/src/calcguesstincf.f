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
      subroutine calcguesstincf(nface,dmin,vfa,umfa,tincfguess)
!
!     calculates a guess for tincf based on the minimum of:
!       dmin/(local velocity)
!       density*dmin**2/(2*dynamic_viscosity)
!
!       dmin is the smallest edge length of the mesh
!
      implicit none
!
      integer nface,i
!
      real*8 vfa(0:5,*),umax,dmin,umfa(*),tincfguess
!
      tincfguess=1.d30
      do i=1,nface
         umax=dsqrt(vfa(1,i)*vfa(1,i)+
     &              vfa(2,i)*vfa(2,i)+
     &              vfa(3,i)*vfa(3,i))
         if(umax.gt.1.d-30) tincfguess=min(tincfguess,dmin/umax)
         tincfguess=min(tincfguess,vfa(5,i)*dmin*dmin/
     &                  (2.d0*umfa(i)))
      enddo
!     
      return
      end
