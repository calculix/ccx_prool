!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine utemp(temp,msecpt,kstep,kinc,time,node,coords,vold,
     &  mi)
!
!     user subroutine utemp
!
!
!     INPUT:
!
!     msecpt             number of temperature values (for volume elements:1)
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     node               node number
!     coords(1..3)       global coordinates of the node
!     vold(0..4,1..nk)   solution field in all nodes 
!                        (not available for CFD-calculations)
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: not used
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!
!     OUTPUT:
!
!     temp(1..msecpt)    temperature in the node
!           
      implicit none
!
      integer msecpt,kstep,kinc,node,mi(*) 
      real*8 temp(msecpt),time(2),coords(3),vold(0:mi(2),*)
!
!     In order to use this user subroutine:
!       1. delete the next call to utemp_ccxlib
!       2. write your own code replacing the line "temp(1)=293.d0"
!
      call utemp_ccxlib(temp,msecpt,kstep,kinc,time,node,coords,vold,
     &  mi)
!
!     Start here your own code. The next line is an example of how your
!     code could look like.
!
c      temp(1)=293.d0
!
      return
      end

