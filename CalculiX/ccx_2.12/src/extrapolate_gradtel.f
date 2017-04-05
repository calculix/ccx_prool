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
      subroutine extrapolate_gradtel(nface,ielfa,xrlfa,gradtel,gradtfa,
     &           icyclic,c,ifatie)
!
!     inter/extrapolation of p element values to the faces
!
      implicit none
!
      integer nface,ielfa(4,*),i,iel1,iel2,l,icyclic,ifatie(*)
!
      real*8 xrlfa(3,*),gradtel(3,*),gradtfa(3,*),xl1,xl2,c(3,3)
!     
c$omp parallel default(none)
c$omp& shared(nface,ielfa,xrlfa,gradtfa,gradtel,icyclic,c,ifatie)
c$omp& private(i,iel1,xl1,iel2,l,xl2)
c$omp do
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
!
!           face in between two elements
!
            xl2=xrlfa(2,i)
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do l=1,3
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+
     &                 xl2*gradtel(l,iel2)
               enddo
            elseif(ifatie(i).gt.0) then
               do l=1,3
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+xl2*
     &                  (gradtel(1,iel2)*c(l,1)+
     &                   gradtel(2,iel2)*c(l,2)+
     &                   gradtel(3,iel2)*c(l,3))
               enddo
            else
               do l=1,3
                  gradtfa(l,i)=xl1*gradtel(l,iel1)+xl2*
     &                  (gradtel(1,iel2)*c(1,l)+
     &                   gradtel(2,iel2)*c(2,l)+
     &                   gradtel(3,iel2)*c(3,l))
               enddo
            endif
         elseif(ielfa(3,i).gt.0) then
!     
!     boundary face; linear extrapolation
!     
            do l=1,3
               gradtfa(l,i)=xl1*gradtel(l,iel1)+
     &              xrlfa(3,i)*gradtel(l,ielfa(3,i))
            enddo
         else
!     
!     boundary face; constant extrapolation (one element layer)
!     
            do l=1,3
               gradtfa(l,i)=gradtel(l,iel1)
            enddo
         endif
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
