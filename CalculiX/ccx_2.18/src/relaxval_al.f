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
      subroutine relaxval_al(gcontfull, nacti, ncdim)

!
      implicit none
C !
      integer i,j, nacti, ncdim, nsize
      integer ksv(ncdim)
C !
      real*8 gcontfull(*), dssd2, r_al(nacti*ncdim) ! ?? TODO: can I declear fixed size if I know it already?
      real*8 gmatrix(nacti*ncdim,nacti*ncdim), vssd1(ncdim)

      logical isDiagDomi ! is strictly diagionally dominant

      nsize = nacti*ncdim

      do i=1,nsize
        do j=1,nsize
            gmatrix(j,i) = gcontfull(  (i-1)*nsize + j  )
        enddo
      enddo

      isDiagDomi=.true.

! WE CHECK IF ITS "STRICTLY DIAGONALLY DOMINANT"
!_strictly_diagonal_dominant = ...
!   all((2*abs(diag(G))) >= sum(abs(G),2));

      do i=1,nsize
!       dssd2=0.0d0
!       do j=1,nsize
!           dssd2 = dssd2 + abs(gmatrix(j,i))
!       enddo
        dssd2 = sum(gmatrix(:,1))

        if ( 2.0d0*gmatrix(i,i) < dssd2 ) then
            isDiagDomi=.false.
            exit
        endif
      enddo

      do i=1,ncdim
        ksv(i) = i
      enddo

      do i=1,nacti
C         II = nCdim*(ii-1)+(1:nCdim);
        do j = 1,ncdim
            ksv(j) = nCdim*(i-1) + j
        enddo

        if(isDiagDomi) then
C             r_al(ii) = 1/max(diag(G(II,II)));
            do j=1,ncdim
                vssd1(j)=gmatrix(ksv(j),ksv(j))
            enddo
        else
C             r_al(ii) = 1/max(sum(abs(G(II,:)),2));
            do j=1,ncdim
                vssd1(j)=sum(abs(gmatrix(:,ksv(j))))!F90
            enddo
        endif

        dssd2 = 1.0d0/maxval(vssd1)!F90

        do j=1,ncdim
            r_al(ksv(j)) = dssd2
        enddo

      enddo

      return
      end
