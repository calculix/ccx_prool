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
      subroutine extrapolate_gradvel(nface,ielfa,xrlfa,gradvel,gradvfa)
!
!     inter/extrapolation of p element values to the faces
!
      implicit none
!
      integer nface,ielfa(4,*),i,iel1,iel2,k,l
!
      real*8 xrlfa(3,*),gradvel(3,3,*),gradvfa(3,3,*),xl1
!     
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face in between two elements
!
            do k=1,3
               do l=1,3
                  gradvfa(k,l,i)=xl1*gradvel(k,l,iel1)+
     &                           xrlfa(2,i)*gradvel(k,l,iel2)
c               write(*,*) 'extrapolate_gradvel1',i,xrlfa(1,i),xrlfa(2,i)
               enddo
            enddo
         elseif(ielfa(3,i).gt.0) then
!
!           boundary face; linear extrapolation
!
            do k=1,3
               do l=1,3
                  gradvfa(k,l,i)=xl1*gradvel(k,l,iel1)+
     &                           xrlfa(3,i)*gradvel(k,l,ielfa(3,i))
c               write(*,*) 'extrapolate_gradvel2',i,iel1,iel2,
c     &              xrlfa(1,i),xrlfa(3,i)
               enddo
            enddo
         else
!
!           boundary face; constant extrapolation (one element layer)
!
            do k=1,3
               do l=1,3
                  gradvfa(k,l,i)=gradvel(k,l,iel1)
               enddo
            enddo
         endif
      enddo
!            
      return
      end
