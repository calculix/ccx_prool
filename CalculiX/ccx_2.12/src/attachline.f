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
!
!     finds the intersection of a straight line through the node
!     with coordinates in pnode and direction vector xn with
!     the face containing 
!     "nterms" nodes with coordinates in field "pneigh" (nterms < 9).
!     cave: the coordinates are stored in pneigh(1..3,*)
!
      subroutine attachline(pneigh,pnode,nterms,ratio,dist,xil,etl,xn)
!
!     Author: Saskia Sitzmann
!
      implicit none
!
      integer nterms,i,j,imin,jmin,idepth,k
!
      real*8 ratio(8),pneigh(3,8),pnode(3),dummy,
     &  a(-1:1,-1:1),xi(-1:1,-1:1),et(-1:1,-1:1),p(3),aold(-1:1,-1:1),
     &  xiold(-1:1,-1:1),etold(-1:1,-1:1),distmin,xiopt,etopt,
     &  dist,xil,etl,xn(3),d(5),th,coeff
!
      idepth=3
      d=(/1.d-2,1.d-4,1.d-6,1.d-9,1.d-14/)
      th=1.d-21
!     
      do k=1,idepth   
!     
!     initialisation
!     
         if(k.eq.1)then
            xiopt=0.0
            etopt=0.0
         else
            xiopt=xi(0,0)
            etopt=et(0,0)
         endif
         do i=-1,1
            do j=-1,1
               xi(i,j)=xiopt+i*d(k)
               et(i,j)=etopt+j*d(k)
               xi(i,j)=min(xi(i,j),1.d0)
               xi(i,j)=max(xi(i,j),-1.d0)
               et(i,j)=min(et(i,j),1.d0)
               et(i,j)=max(et(i,j),-1.d0)
               call distattachline(xi(i,j),et(i,j),pneigh,pnode,a(i,j),
     &              p,ratio,nterms,xn)
            enddo
         enddo
!     
!     minimizing the distance from the face to the node
!     
         do
            distmin=a(0,0)
            imin=0
            jmin=0
            do i=-1,1
               do j=-1,1
                  if(a(i,j).lt.distmin) then
                     distmin=a(i,j)
                     imin=i
                     jmin=j
                  endif
               enddo
            enddo
!     
!     exit if minimum found
!     
            if((imin.eq.0).and.(jmin.eq.0)) exit
!     
            do i=-1,1
               do j=-1,1
                  aold(i,j)=a(i,j)
                  xiold(i,j)=xi(i,j)
                  etold(i,j)=et(i,j)
               enddo
            enddo
!     
            do i=-1,1
               do j=-1,1
                  if((i+imin.ge.-1).and.(i+imin.le.1).and.
     &                 (j+jmin.ge.-1).and.(j+jmin.le.1)) then
                     a(i,j)=aold(i+imin,j+jmin)
                     xi(i,j)=xiold(i+imin,j+jmin)
                     et(i,j)=etold(i+imin,j+jmin)
                  else
                     xi(i,j)=xi(i,j)+imin*d(k)
                     et(i,j)=et(i,j)+jmin*d(k)
!     
                     xi(i,j)=min(xi(i,j),1.d0)
                     xi(i,j)=max(xi(i,j),-1.d0)
                     et(i,j)=min(et(i,j),1.d0)
                     et(i,j)=max(et(i,j),-1.d0)
!     
                     call distattachline(xi(i,j),et(i,j),pneigh,
     &          pnode,a(i,j),p,ratio,nterms,xn)
            endif
          enddo
        enddo
       enddo
       if(distmin.lt.th) exit
      enddo

      dist=a(0,0)
!
      if(nterms.eq.3) then
         xil=(xi(0,0)+1.d0)/2.d0
         etl=(et(0,0)+1.d0)/2.d0
         if(xil+etl.gt.1.d0) then
            dummy=xil
            xil=1.d0-etl
            etl=1.d0-dummy
         endif
      elseif(nterms.eq.4) then
         xil=xi(0,0)
         etl=et(0,0)
      elseif(nterms.eq.6) then
         xil=(xi(0,0)+1.d0)/2.d0
         etl=(et(0,0)+1.d0)/2.d0
         if(xil+etl.gt.1.d0) then
            dummy=xil
            xil=1.d0-etl
            etl=1.d0-dummy
         endif
      elseif(nterms.eq.8) then
         xil=xi(0,0)
         etl=et(0,0)
      endif
!      
      call distattachline(xi(0,0),et(0,0),pneigh,
     &     pnode,a(0,0),p,ratio,nterms,xn)
      coeff= (p(1)-pnode(1))*xn(1)+
     &     (p(2)-pnode(2))*xn(2)+
     &     (p(3)-pnode(3))*xn(3)
!     
      return
      end
      
      
