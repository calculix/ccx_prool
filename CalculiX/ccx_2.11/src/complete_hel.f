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
      integer irow(*),neq,nzs,j,k,l,jdof1,jq(*)
      real*8 hel(3,*),bv(neq,3),auv(*),adv(*)
!
!     off-diagonal terms
!
      do j=1,neq
         do l=jq(j),jq(j+1)-1
            jdof1=irow(l)
!
!           subdiagonal terms
!
            do k=1,3
               hel(k,jdof1)=hel(k,jdof1)-auv(l)*bv(j,k)
            enddo
!
!           superdiagonal terms
!
            do k=1,3
               hel(k,j)=hel(k,j)-auv(l+nzs)*bv(jdof1,k)
            enddo
         enddo
      enddo
!
      return
      end
