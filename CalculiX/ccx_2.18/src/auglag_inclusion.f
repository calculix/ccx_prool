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
      subroutine auglag_inclusion(conttype, gcontfull,nacti, 
     &     ncdim, mufric, atol,rtol, pkvec, kitermax, timek  )

!
      implicit none
C !
      integer i,j, conttype, nacti, ncdim,kitermax

      real*8 gcontfull(*), mufric, atol,rtol, pkvec(*), timek


      return
      end
