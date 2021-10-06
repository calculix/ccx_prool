!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2021 Guido Dhondt
!
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!$
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
      subroutine detectactivecont1(vold, nk, mi, aubi, irowbi, jqbi, 
     &      neqtot, fext, aubb, adbb, ltot,
     &      irowbb, jqbb, auw, iroww, jqw, neqslav, gapdof,
     &      auib, irowib, jqib,  icolbb, nactdof, qtmp, neq)
!
      implicit none
!
      integer irowbi(*),jqbi(*),nk, neqtot, irowbb(*), jqbb(*),
     & iroww(*), jqw(*), neqslav,  irowib(*), jqib(*),
     & icolbb(*),mi(*), nactdof(0:mi(2),*), ltot(*),
     & i,j,icol,idof,irow,jdof, neq(*)
!
      real*8 vold(0:mi(2),*), fext(*), aubb(*), adbb(*), 
     & auw(*),  auib(*),aubi(*), gapdof(*),node,value,qtmp(*)
!    
!     create field qtmp, from vold
!      from sorting nodes to DOF
      do i=1,nk
         do j=1,3
!           write(*,*)'DEBUG i j dof',i,j, nactdof(j,i)
            if(nactdof(j,i).gt.0) then
C               idof = nactdof(j,i)
              qtmp(nactdof(j,i)) = vold(j,i) ! TODO: where is vold coming from ? prediction.c? initial force?
            endif
         enddo
      enddo
!
!      We compute g as qtmp = (Kbi * qtmp) + (Kib*qtmp)
!      to account for the missing terms due to the low triangle structure
!       of the matrices
!
!      calculate Kbi * qtmp
      do i=1,neq(1)
         do j=jqbi(i), jqbi(i+1)-1
!           write(*,*)'DEBUG i j dof',i,j
            value   = aubi(j)
            irow    = irowbi(j)
            gapdof(irow) = gapdof(irow)+value*qtmp(i)
         enddo
      enddo
!
!
!     calculate Kib'*qtmp and add to g.
!     transposed multiplication 
      do i=1,neqtot
         do j=jqib(i),jqib(i+1)-1
            value = auib(j)
            icol  = irowib(j)
            gapdof(i)  = gapdof(i)+value*qtmp(icol)
         enddo
      enddo
!
!     add external force
!
      do i=1,neqtot
         jdof = ltot(i)
         node = int(jdof/10.d0)
         idof = jdof - 10*node
C          value = fext(idof) ! DEBUG
         gapdof(i) = fext(idof) - gapdof(i)
      enddo
!
      return
      end
