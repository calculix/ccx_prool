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
      subroutine applympc(nface,ielfa,is,ie,ifabou,ipompc,vfa,coefmpc,
     &  nodempc,ipnei,neifa,labmpc,xbounact,nactdoh)
!
!     applies MPC's to the faces
!
      implicit none
!
      character*20 labmpc(*)
!
      integer i,j,nface,ielfa(4,*),ipointer,is,ie,ifabou(*),mpc,
     &  ipompc(*),index,iel,iface,nodempc(3,*),ipnei(*),neifa(*),
     &  nactdoh(*),ielorig
!
      real*8 coefmpc(*),denominator,vfa(0:5,*),sum,xbounact(*)
!
      do i=1,nface
         if(ielfa(2,i).ge.0) cycle
         ipointer=-ielfa(2,i)
         do j=is,ie
            if(ifabou(ipointer+j).ge.0) cycle
            mpc=-ifabou(ipointer+j)
            index=ipompc(mpc)
            denominator=coefmpc(index)
            sum=0.d0
            index=nodempc(3,index)
            do
               if(index.eq.0) exit
               if(nodempc(1,index).lt.0) then
c               if((nodempc(3,index).eq.0).and.
c     &            (labmpc(mpc)(6:8).eq.'SPC')) then
!
!                 a negative number refers to a boundary
!                 condition (fields nodeboun, ndirboun..)
!                 resulting from a SPC in local coordinates
!                  
                  sum=sum+coefmpc(index)*xbounact(-nodempc(1,index))
               else
!
!                 face term
!
                  ielorig=int(nodempc(1,index)/10.d0)
                  iel=nactdoh(ielorig)
c                  write(*,*) 'applympc ',mpc,nodempc(1,index),iel
                  iface=nodempc(1,index)-10*ielorig
c                  write(*,*) index
c                  write(*,*) coefmpc(index)
c                  write(*,*) iface
c                  write(*,*) ipnei(iel)
c                  write(*,*) neifa(ipnei(iel)+iface)
c                  write(*,*) nodempc(2,index)
c                  write(*,*) vfa(nodempc(2,index),neifa(ipnei(iel)+iface))
                  sum=sum+coefmpc(index)
     &                 *vfa(nodempc(2,index),neifa(ipnei(iel)+iface))
               endif
               index=nodempc(3,index)
            enddo
            vfa(j,i)=-sum/denominator
c            write(*,*) 'applympc ',ielfa(1,i),ielfa(4,i),j,vfa(j,i)
         enddo
      enddo
!     
      return
      end
