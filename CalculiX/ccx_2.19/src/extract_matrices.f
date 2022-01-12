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
!     MASSLESS DYNAMIC CONTACT

C     >    *MASSLESS DYNAMIC CONTACT*: Extracting the submatrices (i/b) from the global stiffness matrix

C     > Submatrices:
C     > + the bb matrix (neqtot x neqtot)
C     > + the bi matrix (neqtot x neq(1))
C     > + the ib matrix (neq(1) x neqtot)
C     > + the ii matrix (neq(1) x neq(1))
C     >    - this matrix has the same
C     >     size of the global matrix, in which the bi-entries,
C     >     ib-entries, bb-entries (off-diagonal) are set to zero
C     >     the diagonal bb-entries are set to one

C     >     @param au        LOWER triangle of STIFFNESS matrix of size: neq[0]
C     >     @param ad        DIAGONAL       of STIFFNESS matrix of size: neq[0]
C     >     @param jq        Location in field **irow** of the first subdiagonal nonzero in column i (only for symmetric matrices)
C     >     @param irow      Row of element i in field au (i.e. au(i)) 
C     >     @param  neq        NTOT of equations: 
C     >                            + neq[0]: mechanical TOTDOF, 
C     >                            + neq[1]: TOTDOF_mech + TOTDOF_thermal
C     >                            + neq[2]: neq[1] + # of single point constraints (only for modal calculations)
C     >     @param aubb      the param
C     >     @param adbb      the param
C     >     @param jqbb      the param
C     >     @param irowbb    the param
C     >     @param neqtot    the param
C     >     @param nzsbb     the param
C     >     @param aubi      the param
C     >     @param jqbi      the param
C     >     @param irowbi    the param
C     >     @param nzsbi     the param
C     >     @param auib      the param
C     >     @param jqib      the param
C     >     @param irowib    the param
C     >     @param nzsib     the param
C     >     @param ktot       slave and master dofs list

C     >    @see massless.c
C     >    @author Guido Dhondt
!
!     extract matrices bb, bi and ib from the stiffness matrix
!     setting those entries to zero (off-diagonal terms) and one
!     (diagonal terms) to get ii
!
!     all matrices only contain the subdiagonal terms
!
!     to obtain the complete bi-matrix it has to be complement by
!     the transpose of the ib-matrix
!
      subroutine extract_matrices(au,ad,jq,irow,neq,aubb,adbb,jqbb,
     &     irowbb,neqtot,nzsbb,aubi,
     &     jqbi,irowbi,nzsbi,auib,jqib,irowib,nzsib,ktot,icolbb)
!      
      implicit none
!     
      integer jq(*),irow(*),neq(*),jqbb(*),irowbb(*),neqtot,
     &     nzsbb,jqbi(*),irowbi(*),nzsbi,jqib(*),irowib(*),nzsib,
     &     kbb,i,j,ktot(*),id,lbb, icolbb(*)
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
!     column counter for the contact dofs
!
      kbb=1
!      
      do i=1,neq(1)
!     
!     filling submatrix bi
!     
        if(i.ne.ktot(kbb)) then
!
!         i-column
!
          jqbi(i)=nzsbi+1
!
          do j=jq(i),jq(i+1)-1
            call nident(ktot,irow(j),neqtot,id)
            if(id.gt.0) then
              if(ktot(id).eq.irow(j)) then
!
!                b-row: copy entry to bi-matrix; remove from ii-matrix.
!
                nzsbi=nzsbi+1
                aubi(nzsbi)=au(j)
                au(j)=0.d0
                irowbi(nzsbi)=id
              endif
            endif
          enddo
          cycle
        endif
!        
!       b-column
!
!       row counter for the contact dofs
!
        lbb=1
!
        jqbi(i)=nzsbi+1
        jqbb(kbb)=nzsbb+1
        jqib(kbb)=nzsib+1
!     
        loop: do j=jq(i),jq(i+1)-1
          if(lbb.gt.neqtot) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!
            nzsib=nzsib+1
            auib(nzsib)=au(j)
            au(j)=0.d0
            irowib(nzsib)=irow(j)
            cycle
          elseif(irow(j).lt.ktot(lbb)) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!           
            nzsib=nzsib+1
            auib(nzsib)=au(j)
            au(j)=0.d0
            irowib(nzsib)=irow(j)
            cycle
          elseif(irow(j).eq.ktot(lbb)) then
!
!           b-row: copy to bb-matrix and remove from ii-matrix          
!           
            nzsbb=nzsbb+1
            aubb(nzsbb)=au(j)
            au(j)=0.d0
            irowbb(nzsbb)=lbb
            lbb=lbb+1
            cycle
          else
!     
!     actual row entry irow(j) exceeds actual contact dof ktot(lbb).
!     Look for the next contact dof exceeding the actual row entry
!     (the procedure with kbb could be replaced by a call to nident;
!     however, it is assumed that the kbb-procedure is faster)
!     
            do
              lbb=lbb+1
              if(lbb.gt.neqtot) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!
                nzsib=nzsib+1
                auib(nzsib)=au(j)
                au(j)=0.d0
                irowib(nzsib)=irow(j)
                cycle loop
              elseif(irow(j).lt.ktot(lbb)) then
!
!           i-row: copy to ib-matrix and remove from ii-matrix          
!
                nzsib=nzsib+1
                auib(nzsib)=au(j)
                au(j)=0.d0
                irowib(nzsib)=irow(j)
                cycle loop
              elseif(irow(j).eq.ktot(lbb)) then
!
!           b-row: copy to bb-matrix and remove from ii-matrix          
!
                nzsbb=nzsbb+1
                aubb(nzsbb)=au(j)
                au(j)=0.d0
                irowbb(nzsbb)=lbb
                lbb=lbb+1
                cycle loop
              endif
            enddo
          endif
        enddo loop
!
!       b-diagonal row: copy to bb-matrix and set entry in
!       ii-matrix to 1
!
        adbb(kbb)=ad(i)
        ad(i)=1.d0
!        
        if(kbb.lt.neqtot) kbb=kbb+1
      enddo
!
!     finalizing the jq-entries
!
      jqbb(neqtot+1)=nzsbb+1
      jqbi(neq(1)+1)=nzsbi+1
      jqib(neqtot+1)=nzsib+1
!     
      do i=1,neqtot
        icolbb(i)=jqbb(i+1)-jqbb(i)
      enddo
!
      return
      end
