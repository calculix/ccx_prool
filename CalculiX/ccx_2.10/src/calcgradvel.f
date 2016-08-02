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
      subroutine calcgradvel(nef,lakonf,ipnei,vfa,area,xxn,gradvel,
     &  neifa,volume)
!
!     calculation of the gradient of the velocity at the center
!     of the elements from the velocity values at the neighboring
!     faces
!
      implicit none
!
      character*8 lakonf(*)
!
      integer nef,ipnei(*),neifa(*),i,j,k,l,indexf,ifa,numfaces
!
      real*8 vfa(0:5,*),area(*),xxn(3,*),gradvel(3,3,*),volume(*)
!
c$omp parallel default(none)
c$omp& shared(nef,ipnei,lakonf,neifa,gradvel,vfa,area,xxn,volume)
c$omp& private(i,indexf,numfaces,j,ifa,k,l)
c$omp do
      do i=1,nef
         indexf=ipnei(i)
         if(lakonf(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakonf(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
            indexf=indexf+1
            ifa=neifa(indexf)
            do k=1,3
               do l=1,3
                  gradvel(k,l,i)=gradvel(k,l,i)+
     &                 vfa(k,ifa)*area(ifa)*xxn(l,indexf)
               enddo
            enddo
         enddo
!     
!     dividing by the volume of the element
!     
         do k=1,3
            do l=1,3
               gradvel(k,l,i)=gradvel(k,l,i)/volume(i)
            enddo
         enddo
      enddo
c$omp end do
c$omp end parallel
!            
      return
      end
