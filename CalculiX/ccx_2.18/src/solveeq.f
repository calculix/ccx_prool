!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2021 Guido Dhondt
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
      subroutine solveeq(adb,aub,adl,b,sol,aux,irow,jq,
     &  neq,maxit)
!
!     solving a system of equations by iteratively solving the
!     lumped version
!     The diagonal terms f the original system are stored in adb,
!     the off-diagonal terms in aub
!     Ref: The Finite Element Method for Fluid Dynamics,
!          O.C. Zienkiewicz, R.L. Taylor & P. Nithiarasu
!          6th edition (2006) ISBN 0 7506 6322 7
!          p. 61
!
      implicit none
!
      integer irow(*),jq(*),neq,i,k,maxit
!
      real*8 adb(*),aub(*),adl(*),b(*),sol(*),aux(*)
!
!     first iteration
!
      do i=1,neq
         sol(i)=b(i)*adl(i)
      enddo
      if(maxit.eq.1) return
!
!     iterating maxit times
!
      do k=2,maxit
!
!        multiplying the difference of the original matrix
!        with the lumped matrix with the actual solution 
!
         call op(neq,sol,aux,adb,aub,jq,irow)
!
         do i=1,neq
            sol(i)=(b(i)-aux(i))*adl(i)
         enddo
!
      enddo
!
      return
      end
