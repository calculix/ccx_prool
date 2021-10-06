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
!     MASSLESS DYNAMIC CONTACT

C>    *MASSLESS DYNAMIC CONTACT*: Extracting the submatrices (i/b) from the global stiffness matrix

C> Submatrices:
C> + the bb matrix (neqtot x neqtot)
C> + the bi matrix (neqtot x neq(1))
C> + the ib matrix (neq(1) x neqtot)
C> + the ii matrix (neq(1) x neq(1))
C>    - this matrix has the same
C>     size of the global matrix, in which the bi-entries,
C>     ib-entries, bb-entries (off-diagonal) are set to zero
C>     the diagonal bb-entries are set to one

C>     @param au        LOWER triangle of STIFFNESS matrix of size: neq[0]
C>     @param ad        DIAGONAL       of STIFFNESS matrix of size: neq[0]
C>     @param jq        Location in field **irow** of the first subdiagonal nonzero in column i (only for symmetric matrices)
C>     @param irow      Row of element i in field au (i.e. au(i)) 
C>     @param  neq        NTOT of equations: 
C>                            + neq[0]: mechanical TOTDOF, 
C>                            + neq[1]: TOTDOF_mech + TOTDOF_thermal
C>                            + neq[2]: neq[1] + # of single point constraints (only for modal calculations)
C>     @param nzs       see **massless.c**
C>     @param aubb      the param
C>     @param adbb      the param
C>     @param jqbb      the param
C>     @param irowbb    the param
C>     @param neqtot    the param
C>     @param nzsbb     the param
C>     @param islavnode see **massless.c**
C>     @param nslavs    see **massless.c**
C>     @param imastnode see **massless.c**
C>     @param nmasts    NTOT of master nodes
C>     @param aubi      the param
C>     @param jqbi      the param
C>     @param irowbi    the param
C>     @param nzsbi     the param
C>     @param auib      the param
C>     @param jqib      the param
C>     @param irowib    the param
C>     @param nzsib     the param
C>     @param kslav     slave contact dofs list
C>     @param ktot       slave and master dofs list

C>    @see massless.c
C>    @author Guido Dhondt

      subroutine extract_matrices(au,ad,jq,irow,neq,nzs,aubb,adbb,jqbb,
     &     irowbb,neqtot,nzsbb,islavnode,nslavs,imastnode,nmasts,aubi,
     &     jqbi,irowbi,nzsbi,auib,jqib,irowib,nzsib,kslav,ktot, icolbb)
      implicit none
!
      integer jq(*),irow(*),neq(*),nzs(*),jqbb(*),irowbb(*),neqtot,
     &     nzsbb,islavnode(*),nslavs,imastnode(*),nmasts,jqbi(*),
     &     irowbi(*),nzsbi,jqib(*),irowib(*),nzsib,kbb,i,j,kslav(*),
     &     ktot(*),id,lbb, icolbb(*)
!
      real*8 au(*),ad(*),aubb(*),adbb(*),aubi(*),auib(*)
!
!     number of potential nonzero subdiagonal entries
!
      nzsbb=0
      nzsbi=0
      nzsib=0
!
!     treating column by column in the global stiffness matrix
!
      kbb=1
      do i=1,neq(1)
C       write(*,*)'DEBUG i',i
!
!        filling submatrix bi
!
         if(i.ne.ktot(kbb)) then !does not correspond to MAS or SLV node therefore potentially to bi?
            jqbi(i)=nzsbi+1

!           write(*,*)'DEBUG jq(i+1)-1',jq(i+1)-1
            do j=jq(i),jq(i+1)-1 ! all items in column i in  big K
!              write(*,*)'DEBUG j',j
               call nident(ktot,irow(j),neqtot,id) ! if belongs to ktot (then is b) or not (then is i) entry
               if(id.gt.0) then
                  if(ktot(id).eq.irow(j)) then ! ths its bi entry
!                    write(*,*)'(its bi) i,j,ktot(id),: ',i,j,ktot(id)
                     nzsbi         = nzsbi+1
                     aubi(nzsbi)   = au(j)
                     au(j)         = 0.d0
                     irowbi(nzsbi) = id
                  endif
               endif
            enddo
            cycle ! was in bi, therefore we cycle
         endif
!       then its a b column
!        filling submatrix ib and bb
!         
         lbb=1
         jqbi(i)=nzsbi+1! check if workssssss
         jqbb(kbb)=nzsbb+1
         jqib(kbb)=nzsib+1 !
!
!        LOOP over all NONDIAG entries of current columni of i of K matrix
         loop: do j=jq(i),jq(i+1)-1 !QUESTION why the -1 ?
!           write(*,*)'DEBUG j',j
            if(irow(j).lt.ktot(lbb)) then ! its an i row, not a b row.   SMALLER
!               write(*,*)'(its ib (<) )', irow(j),ktot(lbb)
               nzsib=nzsib+1
               auib(nzsib)=au(j)
               au(j)=0.d0
               irowib(nzsib)=irow(j) ! TODO: off by +18 or more
               cycle
            elseif(irow(j).eq.ktot(lbb)) then ! b row its only if its equal EQUAL
!              write(*,*)'(its bb (=))', irow(j),ktot(lbb)
               nzsbb=nzsbb+1
               aubb(nzsbb)=au(j)
               au(j)=0.d0
               irowbb(nzsbb)=lbb
               lbb=lbb+1
               cycle
            else ! BIGGER ( prob wont happen but just in case) --> happens all the time D:
!               write(*,*) 'bigger >! ',irow(j),ktot(lbb)
               do !QUESTION is this a loop?
                  lbb=lbb+1
!                 write(*,*) 'lbb++: ' , lbb !  debug
                  if(irow(j).lt.ktot(lbb)) then
!                    write(*,*)'(its ib (<))', irow(j),ktot(lbb)
                     nzsib=nzsib+1
                     auib(nzsib)=au(j)
                     au(j)=0.d0 ! QUESTION why???
                     irowib(nzsib)=irow(j) ! TODO: off by +18 or more
                     cycle loop ! TODO not sure if w/wo loop ??
                  elseif(irow(j).eq.ktot(lbb)) then
!                    write(*,*)'(its bb (=))', irow(j),ktot(lbb)
                     nzsbb=nzsbb+1
                     aubb(nzsbb)=au(j)
                     au(j)=0.d0 ! QUESTION why??? --> erasing from Kfe which is now Kii
                     irowbb(nzsbb)=lbb
                     lbb=lbb+1
                     cycle loop
                  endif
               enddo
            endif
         enddo loop
         adbb(kbb)=ad(i) ! copy diag at bb
         ad(i)=1.d0 ! setting as 1 in Kfe, (to use at Kii)
         kbb=kbb+1
      enddo
      jqbb(neqtot+1)=nzsbb+1
      jqbi(neq(1)+1)=nzsbi+1
      jqib(neqtot+1)=nzsib+1 ! must be ++1 jq: tot subdiags +1
!

      do i=1,neqtot
         icolbb(i) = jqbb(i+1) - jqbb(i)
      enddo

      return
      end
