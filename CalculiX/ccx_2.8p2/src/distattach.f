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
      subroutine distattach(xig,etg,pneigh,pnode,a,p,ratio,nterms)
!
!     calculates the distance between the node with coordinates
!     in "pnode" and the node with local coordinates xig and etg
!     in a face described by "nterms" nodes with coordinates
!     in pneigh
!
      implicit none
!
      integer nterms,i,j
!
      real*8 ratio(9),pneigh(3,*),pnode(3),a,xi,et,xig,etg,p(3),
     &  dummy,fxi1,fxi2,fxi3,fet1,fet2,fet3,b,xip,xim,etp,etm,
     &  xim2,etm2
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
c         xi=xig
c         et=etg
c         ratio(1)=(1.d0-xi)*(1.d0-et)/4.d0
c         ratio(2)=(1.d0+xi)*(1.d0-et)/4.d0
c         ratio(3)=(1.d0+xi)*(1.d0+et)/4.d0
c         ratio(4)=(1.d0-xi)*(1.d0+et)/4.d0
         xip=1.d0+xig
         xim=1.d0-xig
         etp=1.d0+etg
         etm=1.d0-etg
         ratio(1)=xim*etm/4.d0
         ratio(2)=xip*etm/4.d0
         ratio(3)=xip*etp/4.d0
         ratio(4)=xim*etp/4.d0
      elseif(nterms.eq.6) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
         a=1.d0-xi-et
         ratio(1)=a*(2.d0*a-1.d0)
         ratio(2)=xi*(2.d0*xi-1.d0)
         ratio(3)=et*(2.d0*et-1.d0)
         ratio(4)=4.d0*xi*a
         ratio(5)=4.d0*xi*et
         ratio(6)=4.d0*et*a 
      elseif(nterms.eq.7) then
         xi=(xig+1.d0)/2.d0
         et=(etg+1.d0)/2.d0
         if(xi+et.gt.1.d0) then
            dummy=xi
            xi=1.d0-et
            et=1.d0-dummy
         endif
         a=1.d0-xi-et
         b=a*xi*et
         ratio(1)=a*(2.d0*a-1.d0)+3.d0*b
         ratio(2)=xi*(2.d0*xi-1.d0)+3.d0*b
         ratio(3)=et*(2.d0*et-1.d0)+3.d0*b
         ratio(4)=4.d0*xi*a-12.d0*b
         ratio(5)=4.d0*xi*et-12.d0*b
         ratio(6)=4.d0*et*a-12.d0*b
         ratio(7)=27.d0*b
      elseif(nterms.eq.8) then
c         xi=xig
c         et=etg
c         ratio(1)=(1.d0-xi)*(1.d0-et)*(-xi-et-1.d0)/4.d0
c         ratio(2)=(1.d0+xi)*(1.d0-et)*(xi-et-1.d0)/4.d0
c         ratio(3)=(1.d0+xi)*(1.d0+et)*(xi+et-1.d0)/4.d0
c         ratio(4)=(1.d0-xi)*(1.d0+et)*(-xi+et-1.d0)/4.d0
c         ratio(5)=(1.d0-xi*xi)*(1.d0-et)/2.d0
c         ratio(6)=(1.d0+xi)*(1.d0-et*et)/2.d0
c         ratio(7)=(1.d0-xi*xi)*(1.d0+et)/2.d0
c         ratio(8)=(1.d0-xi)*(1.d0-et*et)/2.d0
         xip=1.d0+xig
         xim=1.d0-xig
         xim2=xip*xim
         etp=1.d0+etg
         etm=1.d0-etg
         etm2=etp*etm
         ratio(1)=xim*etm*(-xig-etp)/4.d0
         ratio(2)=xip*etm*(xig-etp)/4.d0
         ratio(3)=xip*etp*(xig-etm)/4.d0
         ratio(4)=xim*etp*(-xig-etm)/4.d0
         ratio(5)=xim2*etm/2.d0
         ratio(6)=xip*etm2/2.d0
         ratio(7)=xim2*etp/2.d0
         ratio(8)=xim*etm2/2.d0
      elseif(nterms.eq.9) then
         xi=xig
         et=etg
!
         fxi1=xi*(xi-1.d0)/2.d0
         fxi2=(1.d0-xi)*(1.d0+xi)
         fxi3=xi*(xi+1.d0)/2.d0
         fet1=et*(et-1.d0)/2.d0
         fet2=(1.d0-et)*(1.d0+et)
         fet3=et*(et+1.d0)/2.d0
         ratio(1)=fxi1*fet1
         ratio(2)=fxi3*fet1
         ratio(3)=fxi3*fet3
         ratio(4)=fxi1*fet3
         ratio(5)=fxi2*fet1
         ratio(6)=fxi3*fet2
         ratio(7)=fxi2*fet3
         ratio(8)=fxi1*fet2
         ratio(9)=fxi2*fet2
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
      a=(pnode(1)-p(1))**2+(pnode(2)-p(2))**2+(pnode(3)-p(3))**2
!
      return
      end
      
