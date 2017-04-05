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
      subroutine extrapolate_oel(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xboun,ipnei,nef,physcon,umfa,dy)
!
!     inter/extrapolation of the temperature at the center of the elements
!     to the center of the faces
!
      implicit none
!
      integer nface,ielfa(4,*),ifabou(*),iel1,iel2,iel3,i,ipointer,
     &  indexf,ipnei(*),nef
!
      real*8 xrlfa(3,*),vel(nef,0:7),vfa(0:7,*),xboun(*),xl1,xl2,
     &  cinf,physcon(*),umfa(*),dy(*)
!
!     5.5 x 10**(-3.5) = 1.7393e-3
!
      cinf=5.5d0*physcon(5)/physcon(8)
!
c$omp parallel default(none)
c$omp& shared(nface,ielfa,xrlfa,vel,vfa,ipnei,ifabou,xboun,
c$omp&        umfa,cinf,dy)
c$omp& private(i,iel1,xl1,iel2,xl2,indexf,iel3,ipointer)
c$omp do
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            xl2=xrlfa(2,i)
            vfa(7,i)=xl1*vel(iel1,7)+xl2*vel(iel2,7)
         else
            indexf=ipnei(iel1)+ielfa(4,i)
            iel3=ielfa(3,i)
            if(iel2.lt.0) then
               ipointer=-iel2
!
               if(ifabou(ipointer+5).gt.0) then
!     
!                 wall
!     
                  vfa(7,i)=dy(ifabou(ipointer+5))*umfa(i)/vfa(5,i)
               elseif((ifabou(ipointer+1).gt.0).and.
     &                (ifabou(ipointer+2).gt.0).and.
     &                (ifabou(ipointer+3).gt.0)) then
!     
!                 inlet
!     
                  vfa(7,i)=cinf*umfa(i)
               else
!     
!                 constant extrapolation
!     
                  vfa(7,i)=vel(iel1,7)
               endif
!     
            else
!     
!              constant extrapolation
!     
               vfa(7,i)=vel(iel1,7)
            endif
         endif
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
