!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998 Guido Dhondt
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
      subroutine approxplane(col,straight,xn,np)
!
!     calculate the equation of the planes through the
!     edges of a quadrilateral and parallel to the vector xn together
!     with a plane perpendicular to xn and through the center of gravity 
!     of the four corner nodes of the quadrilateral  
!     (so-called mean quadrilateral plane) with 
!     (col(1,1),col(2,1),col(3,1)),(col(1,2),col(2,2),col(3,2)),
!     (col(1,3),col(2,3),col(3,3)),(col(1,4),col(2,4),col(3,4))
!     as vertices. The equation of the planes through the first edge 
!     (connecting the first and the second node) is of the form
!     straight(1)*x+straight(2)*y+straight(3)*z+straight(4)=0, such that the
!     vector (straight(1),straight(2),straight(3)) points outwards (replace 
!     (1) by (5),(9) and (13) for the second, third and fourth edge,
!     similar offset for (2),(3) and (4); 
!     The equation of the mean quadrilateral plane is
!     straight(17)*x+straight(18)*y+straight(19)*z+straight(20)=0 such
!     that the quadrilateral is numbered clockwise when looking in the 
!     direction of vector (straight(17),straight(18),straight(19)).
!
!     adapted for quadratic elements hex20, tet10 
!     Author: Saskia Sitzmann    
!
      implicit none
!     
      integer i,j,modf,np
!
      real*8 col(3,np),colmean(3),straight(36),ps(8,3),dd,xn(3)
!     
!     sides of the quadrilateral
!
      do i=1,3
         do j=1, np
            ps(j,i)=col(i,modf(np,j+1))-col(i,modf(np,j))
         enddo
      enddo
!     
!     mean normal to the quadrilateral (given)
!     
      do i=1,3
         straight(4*np+i)=xn(i)
      enddo
!     
!     ps(j,:) x xn
!     
      do j=1,np
         straight((j-1)*4+1)=ps(j,2)*xn(3)-ps(j,3)*xn(2)
         straight((j-1)*4+2)=ps(j,3)*xn(1)-ps(j,1)*xn(3)
         straight((j-1)*4+3)=ps(j,1)*xn(2)-ps(j,2)*xn(1)
         dd=dsqrt(straight((j-1)*4+1)*straight((j-1)*4+1)
     &        +straight((j-1)*4+2)*straight((j-1)*4+2)
     &        +straight((j-1)*4+3)*straight((j-1)*4+3))
         do i=1,3
            straight((j-1)*4+i)=straight((j-1)*4+i)/dd
         enddo
      enddo
!     
!     determining the inhomogeneous terms
!     
      do i=1,3
         colmean(i)=0.0
      enddo
      do j=1,np
         straight((j-1)*4+4)=-straight((j-1)*4+1)*col(1,j)
     &        -straight((j-1)*4+2)*col(2,j)
     &        -straight((j-1)*4+3)*col(3,j)
         do i=1,3
            colmean(i)=colmean(i)+col(i,j)
         enddo
      enddo
      straight(4*np+4)=-xn(1)*colmean(1)/(real(np)) 
     &     -xn(2)*colmean(2)/(real(np))
     &     -xn(3)*colmean(3)/(real(np))
!     
      return
      end
      
