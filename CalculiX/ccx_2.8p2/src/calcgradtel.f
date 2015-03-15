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
      subroutine calcgradtel(nef,lakon,ipnei,vfa,area,xxn,gradtel,neifa,
     &  volume)
!
!     calculation of the gradient of the velocity at the center
!     of the elements from the velocity values at the neighboring
!     faces
!
      implicit none
!
      character*8 lakon(*)
!
      integer nef,ipnei(*),neifa(*),i,j,l,indexf,ifa,numfaces
!
      real*8 vfa(0:5,*),area(*),xxn(3,*),gradtel(3,*),volume(*)
!
      do i=1,nef
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakon(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
            indexf=indexf+1
            ifa=neifa(indexf)
            do l=1,3
               gradtel(l,i)=gradtel(l,i)+
     &              vfa(0,ifa)*area(ifa)*xxn(l,indexf)
            enddo
         enddo
!     
!     dividing by the volume of the element
!     
         do l=1,3
            gradtel(l,i)=gradtel(l,i)/volume(i)
         enddo
      enddo
!            
      return
      end
