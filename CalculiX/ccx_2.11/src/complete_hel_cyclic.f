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
!     completing hel:
!
!     at the start of the subroutine: rhs of conservation of momentum
!                                     without pressure contribution
!     at the end of the subroutine: neighboring velocity terms subtracted
!
      subroutine complete_hel_cyclic(neq,bv,hel,adv,auv,jq,irow,
     &  ipnei,neiel,ifatie,c,lakonf,neifa,nzs)
!
      implicit none
!
      character*8 lakonf(*)
!
      integer irow(*),neq,nzs,j,k,l,jdof1,jq(*),ifa,neifa(*),numfaces,
     &  indexf,ipnei(*),neiel(1),ifatie(*)
!
      real*8 hel(3,*),bv(neq,3),auv(*),adv(*),c(3,3)
!
!     off-diagonal terms
!
c      write(*,*) 'complete_hel_cyclic',nzs
c      do j=1,8
c         write(*,*) j,lakonf(j)
c      enddo
c      do j=1,48
c         write(*,*) j,neifa(j)
c      enddo
c      do j=1,8
c         write(*,*) j,ipnei(j)
c      enddo
c      do j=1,48
c         write(*,*) j,neiel(j)
c      enddo
c      do j=1,38
c         write(*,*) j,ifatie(j)
c      enddo
      do j=1,neq
         do l=jq(j),jq(j+1)-1
            jdof1=irow(l)
!
!           determining the face between elements jdof1 and j
!
            indexf=ipnei(jdof1)
            if(lakonf(jdof1)(4:4).eq.'8') then
               numfaces=6
            elseif(lakonf(jdof1)(4:4).eq.'6') then
               numfaces=5
            else
               numfaces=4
            endif
            do k=1,numfaces
               indexf=indexf+1
               if(neiel(indexf).eq.j) exit
            enddo
            ifa=neifa(indexf)
            if(ifatie(ifa).eq.0) then
!
!              no cyclic symmetry face
!
!              subdiagonal terms
!
               do k=1,3
                  hel(k,jdof1)=hel(k,jdof1)-auv(l)*bv(j,k)
               enddo
!     
!              superdiagonal terms
!
               do k=1,3
                  hel(k,j)=hel(k,j)-auv(l+nzs)*bv(jdof1,k)
               enddo
!
            elseif(ifatie(ifa).gt.0) then
!
!              subdiagonal terms
!
               do k=1,3
                  hel(k,jdof1)=hel(k,jdof1)-auv(l)*
     &                 (c(k,1)*bv(j,1)+c(k,2)*bv(j,2)+c(k,3)*bv(j,3))
               enddo
!     
!              superdiagonal terms
!
               do k=1,3
                  hel(k,j)=hel(k,j)-auv(l+nzs)*
     &                 (c(1,k)*bv(jdof1,1)+c(2,k)*bv(jdof1,2)
     &                 +c(3,k)*bv(jdof1,3))
               enddo
            else
!
!              subdiagonal terms
!
               do k=1,3
                  hel(k,jdof1)=hel(k,jdof1)-auv(l)*
     &                 (c(1,k)*bv(j,1)+c(2,k)*bv(j,2)+c(3,k)*bv(j,3))
               enddo
!     
!              superdiagonal terms
!
               do k=1,3
                  hel(k,j)=hel(k,j)-auv(l+nzs)*
     &                 (c(k,1)*bv(jdof1,1)+c(k,2)*bv(jdof1,2)
     &                 +c(k,3)*bv(jdof1,3))
               enddo
            endif
         enddo
      enddo
!
      return
      end
