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
      subroutine complete_hel(neq,bv,hel,adv,auv,jq,irow,nzs)
!
      implicit none
!
      integer irow(*),neq,nzs,j,l,jdof1,jq(*)
      real*8 hel(3,*),bv(neq,3),auv(*),adv(*)
c      write(*,*) 'complete_hel1 ',hel(1,1),hel(2,1),hel(3,1)
c      do j=1,neq
c      write(*,*) 'complete_hel_1 ',j,hel(1,j),hel(2,j),hel(3,j)
c      enddo
!
!     off-diagonal terms
!
      do j=1,neq
         do l=jq(j),jq(j+1)-1
            jdof1=irow(l)
!
!           subdiagonal terms
!
            hel(1,jdof1)=hel(1,jdof1)-auv(l)*bv(j,1)
            hel(2,jdof1)=hel(2,jdof1)-auv(l)*bv(j,2)
            hel(3,jdof1)=hel(3,jdof1)-auv(l)*bv(j,3)
!
!           superdiagonal terms
!
            hel(1,j)=hel(1,j)-auv(l+nzs)*bv(jdof1,1)
            hel(2,j)=hel(2,j)-auv(l+nzs)*bv(jdof1,2)
            hel(3,j)=hel(3,j)-auv(l+nzs)*bv(jdof1,3)
         enddo
      enddo
c      write(*,*) 'complete_hel2 ',hel(1,1),hel(2,1),hel(3,1)
c      do j=1,neq
c         write(*,*) 'complete_hel ',j,bv(j,1),bv(j,2),bv(j,3)
c      enddo
!
      return
      end
