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
      subroutine distattachline(xig,etg,pneigh,pnode,a,p,
     & ratio,nterms,xn)
!
!     calculates the distance between a straight line through the node 
!     with coordinates in "pnode" and direction vector "xn" and 
!     the node with local coordinates xig and etg
!     in a face described by "nterms" nodes with coordinates
!     in "pneigh"
!
      implicit none
!
      integer nterms,i,j
!
      real*8 ratio(8),pneigh(3,*),pnode(3),a,xi,et,xig,etg,p(3),
     &  dummy,xn(3),coeff
!
      if(nterms.eq.3) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
         ratio(1)=1.d0-xi-et
         ratio(2)=xi
         ratio(3)=et
      elseif(nterms.eq.4) then
         xi=xig
         et=etg
         ratio(1)=(1.d0-xi)*(1.d0-et)/4.d0
         ratio(2)=(1.d0+xi)*(1.d0-et)/4.d0
         ratio(3)=(1.d0+xi)*(1.d0+et)/4.d0
         ratio(4)=(1.d0-xi)*(1.d0+et)/4.d0
      elseif(nterms.eq.6) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
         ratio(1)=2.d0*(0.5d0-xi-et)*(1.d0-xi-et)
         ratio(2)=xi*(2.d0*xi-1.d0)
         ratio(3)=et*(2.d0*et-1.d0)
         ratio(4)=4.d0*xi*(1.d0-xi-et)
         ratio(5)=4.d0*xi*et
         ratio(6)=4.d0*et*(1.d0-xi-et)  
      elseif(nterms.eq.8) then
         xi=xig
         et=etg
         ratio(1)=(1.d0-xi)*(1.d0-et)*(-xi-et-1.d0)/4.d0
         ratio(2)=(1.d0+xi)*(1.d0-et)*(xi-et-1.d0)/4.d0
         ratio(3)=(1.d0+xi)*(1.d0+et)*(xi+et-1.d0)/4.d0
         ratio(4)=(1.d0-xi)*(1.d0+et)*(-xi+et-1.d0)/4.d0
         ratio(5)=(1.d0-xi*xi)*(1.d0-et)/2.d0
         ratio(6)=(1.d0+xi)*(1.d0-et*et)/2.d0
         ratio(7)=(1.d0-xi*xi)*(1.d0+et)/2.d0
         ratio(8)=(1.d0-xi)*(1.d0-et*et)/2.d0
      else
         write(*,*) '*ERROR in distattach: case with ',nterms
         write(*,*) '       terms is not covered'
         call exit(201)
      endif
!
!     calculating the position in the face
!      
      do i=1,3
         p(i)=0.d0
         do j=1,nterms
            p(i)=p(i)+ratio(j)*pneigh(i,j)
         enddo
      enddo
!
!     calculating the distance
!
c      a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+(pnode(3)-p(3))**2
      coeff=0.0
      do i=1,3
         coeff=coeff+xn(i)*(p(i)-pnode(i))
      enddo
      a=(p(1)-pnode(1)-coeff*xn(1))**2+(p(2)-pnode(2)-
     &     coeff*xn(2))**2+(p(3)-pnode(3)-coeff*xn(3))**2
!     
      return
      end
      
