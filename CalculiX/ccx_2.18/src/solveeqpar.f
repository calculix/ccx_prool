!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine solveeqpar(aub,adl,b,sol,aux,irow,jq,neqa,neqb)
!     
!     solving a system of equations by iteratively solving the
!     lumped version
!     Ref: The Finite Element Method for Fluid Dynamics,
!     O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!     6th edition (2006) ISBN 0 7506 6322 7
!     p. 61
!     
      implicit none
!     
      integer irow(*),jq(*),i,j,neqa,neqb
!     
      real*8 aub(*),adl(*),b(*),sol(*),aux(*)
!     
!     multiplying the difference of the original matrix
!     with the lumped matrix with the actual solution 
!     
      do i=neqa,neqb
        aux(i)=0.d0
        do j=jq(i),jq(i+1)-1
          aux(i)=aux(i)+aub(j)*sol(irow(j))
        enddo
        sol(i)=(b(i)-aux(i))*adl(i)
      enddo
!     
      return
      end
