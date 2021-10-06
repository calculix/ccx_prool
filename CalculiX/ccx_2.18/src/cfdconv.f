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
      subroutine cfdconv(vold,vcon,v,nk,nmethod,iconvergence,ithermal,
     &     iit,iturbulent,mi,dtimef,vconmax,iexplicit)
!     
!     calculates the change in solution
!     
      implicit none
!     
      integer iconvergence,iit,iturbulent,mi(*),nk,ithermal(*),i,j,
     &     nmethod,iexplicit
!     
      real*8 v(0:mi(2),*),vold(0:mi(2),*),vcon(0:mi(2),*),vmax(0:6),
     &     vconmax(0:6),ratio(0:6),dtimef
!     
!     first subiteration: calculate the size of the conservative
!     fields
!     
      do j=0,6
        vconmax(j)=0.d0
      enddo
!     
      if(iexplicit.eq.1) then
!     
!     for incompressible fluids the density is stored
!     in vcon(4,*), the change in density in v(4,*)
!     
        do i=1,nk
          do j=0,4
            vconmax(j)=vconmax(j)+vcon(j,i)**2
          enddo
        enddo
      else
        do i=1,nk
          do j=0,3
            vconmax(j)=vconmax(j)+vcon(j,i)**2
          enddo
!     
!     for incompressible fluids the pressure is stored
!     in vold(4,*), the change in pressure in v(4,*)
!     
          vconmax(4)=vconmax(4)+vold(4,i)**2
        enddo
      endif
      if(iturbulent.ne.0) then
        do i=1,nk
          do j=5,6
            vconmax(j)=vconmax(j)+vcon(j,i)**2
          enddo
        enddo
      endif
!     
!     all subiterations: calculate the size of the change of
!     the conservative variables
!     
      do j=0,6
        vmax(j)=0.d0
      enddo
!     
      if(iexplicit.eq.1) then
!     
!     for incompressible fluids the density is stored
!     in vcon(4,*), the change in density in v(4,*)
!     
        do i=1,nk
          do j=0,4
            vmax(j)=vmax(j)+v(j,i)**2
          enddo
        enddo
      else
        do i=1,nk
          do j=0,3
            vmax(j)=vmax(j)+v(j,i)**2
          enddo
!     
!     for incompressible fluids the pressure is stored
!     in vold(4,*), the change in pressure in v(4,*)
!     
          vmax(4)=vmax(4)+v(4,i)**2
        enddo
      endif
      if(iturbulent.ne.0) then
        do i=1,nk
          do j=5,6
            vmax(j)=vmax(j)+v(j,i)**2
          enddo
        enddo
      endif
!     
!     check convergence
!     
      if(iturbulent.eq.0) then
!     
!     laminar
!     
        do i=0,4
          vconmax(i)=dsqrt(vconmax(i))
        enddo
        do i=0,4
          vmax(i)=dsqrt(vmax(i))
          if(vconmax(i).lt.1.d-10) then
            ratio(i)=0.d0
            vmax(i)=0.d0
          else
            ratio(i)=vmax(i)/vconmax(i)
          endif
        enddo
        if(nmethod.eq.1) then
          if(((vmax(0).lt.1.d-8*vconmax(0)).or.
     &         (vconmax(0).lt.1.d-10)).and.
     &         ((vmax(1).lt.1.d-8*vconmax(1)).or.
     &         (vconmax(1).lt.1.d-10)).and.
     &         ((vmax(2).lt.1.d-8*vconmax(2)).or.
     &         (vconmax(2).lt.1.d-10)).and.
     &         ((vmax(3).lt.1.d-8*vconmax(3)).or.
     &         (vconmax(3).lt.1.d-10)).and.
     &         ((vmax(4).lt.1.d-8*vconmax(4)).or.
     &         (vconmax(4).lt.1.d-10)).and.
     &         (iit.gt.1)) iconvergence=1
        endif
        if(iit.gt.1)
     &       write(12,'(i7,15(1x,e10.3))') iit,ratio(0),
     &       ratio(1),ratio(2),
     &       ratio(3),ratio(4),
     &       dtimef
      else
!     
!     turbulent
!     
        do i=0,6
          vconmax(i)=dsqrt(vconmax(i))
        enddo
        do i=0,6
          vmax(i)=dsqrt(vmax(i))
          if(vconmax(i).lt.1.d-10) then
            ratio(i)=0.d0
            vmax(i)=0.d0
          else
            ratio(i)=vmax(i)/vconmax(i)
          endif
        enddo
        if(ithermal(1).eq.0) vconmax(0)=1.d0
        if(nmethod.eq.1) then
          if(((vmax(0).lt.1.d-8*vconmax(0)).or.
     &         (vconmax(0).lt.1.d-10)).and.
     &         ((vmax(1).lt.1.d-8*vconmax(1)).or.
     &         (vconmax(1).lt.1.d-10)).and.
     &         ((vmax(2).lt.1.d-8*vconmax(2)).or.
     &         (vconmax(2).lt.1.d-10)).and.
     &         ((vmax(3).lt.1.d-8*vconmax(3)).or.
     &         (vconmax(3).lt.1.d-10)).and.
     &         ((vmax(4).lt.1.d-8*vconmax(4)).or.
     &         (vconmax(4).lt.1.d-10)).and.
     &         ((vmax(5).lt.1.d-8*vconmax(5)).or.
     &         (vconmax(5).lt.1.d-10)).and.
     &         ((vmax(6).lt.1.d-8*vconmax(6)).or.
     &         (vconmax(6).lt.1.d-10)).and.
     &         (iit.gt.1)) iconvergence=1
        endif
        if(iit.gt.1)
     &       write(12,'(i7,15(1x,e10.3))') iit,ratio(0),
     &       ratio(1),ratio(2),
     &       ratio(3),ratio(4),
     &       ratio(5),ratio(6),dtimef
      endif
!     
      return
      end
      
