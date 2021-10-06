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
      subroutine detectactivecont2(gapnorm, gapdof, auw, iroww, jqw,
     &            neqtot, nslavs, springarea , iacti ,nacti)
C       WE compute g_Npre = g0 +Wb[:,1:3:end]'.gapdof
!
      implicit none
C !
      integer i, inorm, j,icol,iroww(*), jqw(*),neqtot,
     & nslavs, iacti(*), nacti
C !
      real*8 gapnorm(*), gapdof(*), auw(*), springarea(2,*),value

      do i=1,nslavs
            inorm= 3*(i-1) + 1 !normal contact index
            write(*,*) i, inorm
            do j=jqw(inorm),jqw(inorm+1)-1
                value = auw(j)
                icol  = iroww(j)
                gapnorm(i) = gapnorm(i)+value*gapdof(icol)
            enddo
      enddo

      nacti=0
      do i=1,nslavs
!       add initial clearance at time step 0
        gapnorm(i) = gapnorm(i)+springarea(2,i)
!       contact evaluation
        if(gapnorm(i)<= 0.d0)then
            iacti(i) = i
            nacti = nacti+1
        endif
      enddo

      return
      end
