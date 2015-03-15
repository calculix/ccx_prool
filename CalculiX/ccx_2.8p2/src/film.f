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
      subroutine film(h,sink,temp,kstep,kinc,time,noel,npt,
     &  coords,jltyp,field,nfield,loadtype,node,area,vold,mi)
!
!     user subroutine film
!
!
!     INPUT:
!
!     sink               most recent sink temperature
!     temp               current temperature value
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        11 = face 1 
!                        12 = face 2 
!                        13 = face 3 
!                        14 = face 4 
!                        15 = face 5 
!                        16 = face 6
!     field              currently not used
!     nfield             currently not used (value = 1)
!     loadtype           load type label
!     node               network node (only for forced convection)
!     area               area covered by the integration point
!     vold(0..4,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!
!     OUTPUT:
!
!     h(1)               magnitude of the film coefficient
!     h(2)               not used; please do NOT assign any value
!     sink               (updated) sink temperature (need not be
!                        defined for forced convection)
!           
      implicit none
!
      character*20 loadtype
      integer kstep,kinc,noel,npt,jltyp,nfield,node,mi(*)
      real*8 h(2),sink,time(2),coords(3),temp,field(nfield),area,
     &  vold(0:mi(2),*)
!
c      h(1)=10.d0*coords(2)
      h(1)=200.d0
c      sink=-12.d0
!
      return
      end

