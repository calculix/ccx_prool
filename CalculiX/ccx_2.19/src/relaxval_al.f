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
      subroutine relaxval_al(r,gmatrix,nacti,ncdim)
!     
      implicit none
!     
      integer i,j,nacti,ncdim,nsize,kcol(ncdim)
!     
      real*8 dssd2,r(*),vssd1(ncdim),gmatrix((nacti*ncdim),*)
!
      logical isDiagDomi
!
      nsize=nacti*ncdim
!
      isDiagDomi=.true.
!
!     check for diagonal dominance
!
      do i=1,nsize
        dssd2=sum(gmatrix(:,i))
        if(2.0d0*gmatrix(i,i).lt.dssd2) then
          isDiagDomi=.false.
          exit
        endif
      enddo
!
      do i=1,nacti
!
!       correct column
!
        do j=1,ncdim
          kcol(j)=ncdim*(i-1)+j
        enddo
!        
        if(isDiagDomi) then
          do j=1,ncdim
            vssd1(j)=gmatrix(kcol(j),kcol(j))
          enddo
        else
          do j=1,ncdim
            vssd1(j)=sum(abs(gmatrix(:,kcol(j))))
          enddo
        endif
!        
        dssd2=1.0d0/maxval(vssd1)
!        
        do j=1,ncdim
          r(kcol(j))=dssd2
        enddo
      enddo
!
      return
      end
