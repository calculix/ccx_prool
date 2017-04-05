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
      subroutine extrapolate_tel(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xboun,ipnei,nef)
!
!     inter/extrapolation of the temperature at the center of the elements
!     to the center of the faces
!
      implicit none
!
      integer nface,ielfa(4,*),ifabou(*),iel1,iel2,iel3,i,ipointer,
     &  indexf,ipnei(*),nef
!
      real*8 xrlfa(3,*),vel(nef,0:7),vfa(0:7,*),xboun(*),xl1,xl2
!
c$omp parallel default(none)
c$omp& shared(nface,ielfa,xrlfa,vel,vfa,ipnei,ifabou,xboun)
c$omp& private(i,iel1,xl1,iel2,xl2,indexf,iel3,ipointer)
c$omp do
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            xl2=xrlfa(2,i)
            vfa(0,i)=xl1*vel(iel1,0)+xl2*vel(iel2,0)
         else
            indexf=ipnei(iel1)+ielfa(4,i)
            iel3=ielfa(3,i)
            if(iel2.lt.0) then
               ipointer=-iel2
!
               if(ifabou(ipointer).gt.0) then
!     
!                 temperature given
!     
                  vfa(0,i)=xboun(ifabou(ipointer))
               elseif((ifabou(ipointer+6).ne.0).and.(iel3.ne.0)) then
!     
!                 flux given; linear interpolation
!     
                  vfa(0,i)=xl1*vel(iel1,0)+xrlfa(3,i)*vel(iel3,0)
               else
!     
!                 constant extrapolation
!     
                  vfa(0,i)=vel(iel1,0)
               endif
!     
            else
!     
!              constant extrapolation
!     
               vfa(0,i)=vel(iel1,0)
            endif
         endif
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
