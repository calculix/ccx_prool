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
      subroutine extrapolate_pel(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xboun,nef)
!
!     extrapolation of p element values to the faces
!
      implicit none
!
      integer nface,ielfa(4,*),ifabou(*),i,iel1,iel2,nef
!
      real*8 xrlfa(3,*),vel(nef,0:7),vfa(0:7,*),xboun(*),xl1
!     
c      write(*,*) 'extrapolate_pel ',vel(1,4),vel(2,4)
c$omp parallel default(none)
c$omp& shared(nface,ielfa,xrlfa,vfa,vel,ifabou,xboun)
c$omp& private(i,iel1,xl1,iel2)
c$omp do
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            vfa(4,i)=xl1*vel(iel1,4)+xrlfa(2,i)*vel(iel2,4)
         elseif(iel2.lt.0) then
            if(ifabou(-iel2+4).gt.0) then
!
!              boundary face; p given
!               
               vfa(4,i)=xboun(ifabou(-iel2+4))
            elseif(ielfa(3,i).gt.0) then
!
!              boundary face; linear extrapolation
!
               vfa(4,i)=xl1*vel(iel1,4)+xrlfa(3,i)*vel(ielfa(3,i),4)
            else
!
!              boundary face; constant extrapolation
!
               vfa(4,i)=vel(iel1,4)
            endif
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; linear extrapolation
!
            vfa(4,i)=xl1*vel(iel1,4)+xrlfa(3,i)*vel(ielfa(3,i),4)
         else
!
!           boundary face; constant extrapolation
!
            vfa(4,i)=vel(iel1,4)
         endif
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
