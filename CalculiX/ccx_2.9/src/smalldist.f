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
      subroutine smalldist(co,distmin,lakon,ipkon,kon,ne)
!
      implicit none
!
      character*8 lakon(*)
!
      integer ipkon(*),kon(*),i,j,ne,
     &  inoelfree,nope,indexe,node
!
!     New variables for smalldist calculation
!
      integer k,neighbor,noneigh(3,8)
!
      real*8 dist,distmin,co(3,*)
!
!     determining the distance between nodes 
!     
      distmin=1.d6
      inoelfree=1
      do i=1,ne
         if(ipkon(i).lt.0) cycle
         if(lakon(i)(1:1).eq.'F') cycle
         if(lakon(i)(4:4).eq.'2') then
C           nope=20
            nope=8
            noneigh=
     &          reshape((/9,12,17,9,10,18,10,11,19,11,12,20,13,16,17,
     &           13,14,18,14,15,19,15,16,20/),(/3,8/))
         elseif(lakon(i)(1:8).eq.'ESPRNGA1') then
            nope=2
            noneigh=reshape((/2,2,2,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     &           0,0,0/),(/3,8/))
         elseif(lakon(i)(4:4).eq.'8') then
C     nope=8
            nope=8
            noneigh=
     &         reshape((/2,4,5,1,3,6,2,4,7,1,3,8,1,6,8,2,5,7,3,6,8,4,
     &           5,7/),(/3,8/))
         elseif(lakon(i)(4:4).eq.'4') then
C     nope=4
            nope=4
            noneigh=
     &        reshape((/2,3,4,1,3,4,1,2,4,1,2,3,0,0,0,0,0,0,0,0,0,0,
     &           0,0/),(/3,8/))
         elseif(lakon(i)(4:5).eq.'10') then
C     nope=10
            nope=4
            noneigh=
     &         reshape((/5,7,8,5,6,9,6,7,10,8,9,10,0,0,0,0,0,0,0,0,0,
     &           0,0,0/),(/3,8/))
         elseif(lakon(i)(4:4).eq.'6') then
C     nope=6
            nope=6
            noneigh=
     &         reshape((/2,3,4,1,3,5,1,2,6,1,5,6,2,4,6,3,4,5,0,0,0,0,
     &           0,0/),(/3,8/))
         elseif(lakon(i)(4:5).eq.'15') then
C     nope=15
            nope=6
            noneigh=reshape((/7,9,13,7,8,14,8,9,15,10,12,13,10,11,14,11,
     &           12,14,0,0,0,0,0,0/),(/3,8/))
         elseif((lakon(i)(1:2).eq.'ES').or.
     &           (lakon(i)(1:2).eq.'ED')) then
            read(lakon(i)(8:8),'(i1)') nope
            nope=nope+1
         else
            cycle
         endif
         indexe=ipkon(i)
         do j=1,nope
            node=kon(indexe+j)
            do k=1,3
               neighbor=kon(indexe+noneigh(k,j))
               dist=(co(1,node)-co(1,neighbor))**2+(co(2,node)-
     &              co(2,neighbor))**2+(co(3,node)-co(3,neighbor))**2
               distmin=min(dist,distmin)
!     following syntax just for saving node number of smallest distance 
!     if(distmin.eq.dist) then
!     mindistnod(1,1)=node
!     mindistnod(2,1)=neighbor
!     endif
            enddo
         enddo
      enddo
!     
!     1% of the smallest distance used for the variation
!     
      distmin=(distmin**0.5)*0.01
!     
!     Write results in file
!     
!     open(16,file='results_distance.out',status='unknown')
!     write(16,'(1x,i6,1x,i6,1x,f5.4)') mindistnod(1,1),
!     &	mindistnod(2,1),distmin
!     close(16)  
      
      return
      end
