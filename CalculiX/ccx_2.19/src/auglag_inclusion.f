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
      subroutine auglag_inclusion(gmatrix,cvec,iacti,nacti,
     &     ncdim,mufric,atol,rtol,pkglob,kitermax,
     &     auw,jqw,iroww,nslavs,al,alnew,r)
!     
      implicit none
      include 'omp_lib.h'
!
      logical iscvg
!
      character*1 uplo
!     
      integer i,j,k,jj,nacti,ncdim,kitermax,iacti(*),incx,incy,
     &     icont,nsize,inorm,irow,jqw(*),iroww(*),nslavs
!     
!     Al IS in local contact coordinates!
!     pkglob is in global coordinates
!     
      real*8 mufric,atol,rtol,pkglob(*),alsize,al(*),err,
     &     alnew(*),r(*),cvec(*),gmatrix(nacti*ncdim,*),
     &     omega,value,auw(*),alpha,beta
!
      nsize=nacti*ncdim
      iscvg=.false.
!     
!     arguments for dsymv
!     
      uplo='U'
      alpha=1.d0
      incx=1
      beta=0.d0
      incy=1
!
      err=1.d30
      omega=1.d0
      icont=0
!
!     determine the relaxation parameter
!
      call relaxval_al(r,gmatrix,nacti,ncdim)
!
      do while((icont.le.kitermax).and.(.not.(iscvg)))
!        
!     G*lam via BLAS symmetric matrix vector multiplication
!     BLAS is part of the ARPACK library
!
        call dsymv(uplo,nsize,alpha,gmatrix,
     &       nsize,al,incx,beta,alnew,incy)
!
!     al-omega*r*(alnew+c)
!
        do i=1,nsize
          alnew(i)=al(i)-omega*r(i)*(alnew(i)+cvec(i))
        enddo
!
!     only NORMAL CONTACT PROJECTION p = max(p,0)
!
        if(ncdim.eq.1) then
          do i=1,nsize
            alnew(i)=max(alnew(i),0.d0)
          enddo
        else
          write(*,*) '*ERROR in auglag_inclusion.f: contact type'
          write(*,*) '      other than 1 (NORMAL) is not yet '
          write(*,*) '      implemented.'
          call exit(201)
        endif
!
!       determining the change in solution
!
        err=0.d0
        alsize=0.d0
!        
        do i=1,nsize
          err=err+(alnew(i)-al(i))**2
          alsize=alsize+al(i)**2
        enddo
!        
        err=dsqrt(err)
        alsize=dsqrt(alsize)
!
        do i=1,nsize
          al(i)=alnew(i)
        enddo
!
!       check for convergence          
!          
        if(err.le.(alsize*rtol+atol)) then
          iscvg=.true.
        endif
!
        icont=icont+1
      enddo
!
      if(icont.gt.kitermax)then
        write(*,*) '*WARNING: maximum iterations for massless'
        write(*,*) '          contact solution reached:' ,kitermax
        call exit(201)
      else
        write(*,*) 'AL converged! it,err:',icont,err
      endif
!
!     Expansion of pk to global cys: pkglob = Wb*al
!
      do i=1,nslavs
        if(iacti(i).ne.0)then
          inorm=3*(i-1) + 1
          do j=0,(ncdim-1)    !N or NTT
            jj=inorm+j
            do k=jqw(jj),jqw(jj+1)-1 ! each row
              value=auw(k)
              irow=iroww(k)
              pkglob(irow)=pkglob(irow)+value*al(iacti(i)+j)
            enddo
          enddo
        endif
      enddo
!
      return
      end
