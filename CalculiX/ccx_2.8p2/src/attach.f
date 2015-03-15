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
      subroutine attach(pneigh,pnode,nterms,ratio,dist,xil,etl)
!
!     ataches node with coordinates in "pnode" to the face containing 
!     "nterms" nodes with coordinates in field "pneigh" (nterms < 9).
!     cave: the coordinates are stored in pneigh(1..3,*)
!
      implicit none
!
      integer nterms,i,j,imin,jmin
!
      real*8 ratio(9),pneigh(3,9),pnode(3),dummy,
     &  a(-1:1,-1:1),xi(-1:1,-1:1),et(-1:1,-1:1),p(3),aold(-1:1,-1:1),
     &  xiold(-1:1,-1:1),etold(-1:1,-1:1),distmin,xiopt,etopt,
     &  d1,d2,d3,d4,dist,xil,etl
!
c      d1=0.25d0
c      d2=3.125d-2
c      d3=3.9063d-3
c      d4=1.d-3
      d1=1.d-2
      d2=1.d-4
      d3=1.d-6
      d4=1.d-8
!
!     initialisation
!
      do i=-1,1
        do j=-1,1
          xi(i,j)=i*d1
          et(i,j)=j*d1
          call distattach(xi(i,j),et(i,j),pneigh,pnode,a(i,j),p,
     &      ratio,nterms)
c             write(*,*) i,j,a(i,j)
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
!       exit if minimum found
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
     &         (j+jmin.ge.-1).and.(j+jmin.le.1)) then
              a(i,j)=aold(i+imin,j+jmin)
              xi(i,j)=xiold(i+imin,j+jmin)
              et(i,j)=etold(i+imin,j+jmin)
            else
              xi(i,j)=xi(i,j)+imin*d1
              et(i,j)=et(i,j)+jmin*d1
!
              xi(i,j)=min(xi(i,j),1.d0)
              xi(i,j)=max(xi(i,j),-1.d0)
              et(i,j)=min(et(i,j),1.d0)
              et(i,j)=max(et(i,j),-1.d0)
!
              call distattach(xi(i,j),et(i,j),pneigh,
     &          pnode,a(i,j),p,ratio,nterms)
!              write(*,*) a(i,j)
            endif
          enddo
        enddo
      enddo
!
!     2nd run
!     initialisation
!
      xiopt=xi(0,0)
      etopt=et(0,0)
      do i=-1,1
        do j=-1,1
          xi(i,j)=xiopt+i*d2
          et(i,j)=etopt+j*d2
          xi(i,j)=min(xi(i,j),1.d0)
          xi(i,j)=max(xi(i,j),-1.d0)
          et(i,j)=min(et(i,j),1.d0)
          et(i,j)=max(et(i,j),-1.d0)
          call distattach(xi(i,j),et(i,j),pneigh,pnode,a(i,j),p,
     &       ratio,nterms)
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
!       exit if minimum found
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
     &         (j+jmin.ge.-1).and.(j+jmin.le.1)) then
              a(i,j)=aold(i+imin,j+jmin)
              xi(i,j)=xiold(i+imin,j+jmin)
              et(i,j)=etold(i+imin,j+jmin)
            else
              xi(i,j)=xi(i,j)+imin*d2
              et(i,j)=et(i,j)+jmin*d2
!
              xi(i,j)=min(xi(i,j),1.d0)
              xi(i,j)=max(xi(i,j),-1.d0)
              et(i,j)=min(et(i,j),1.d0)
              et(i,j)=max(et(i,j),-1.d0)
!
              call distattach(xi(i,j),et(i,j),pneigh,
     &          pnode,a(i,j),p,ratio,nterms)
!              write(*,*) a(i,j)
            endif
          enddo
        enddo
      enddo
!
!     3rd run
!     initialisation
!
      xiopt=xi(0,0)
      etopt=et(0,0)
      do i=-1,1
        do j=-1,1
          xi(i,j)=xiopt+i*d3
          et(i,j)=etopt+j*d3
          xi(i,j)=min(xi(i,j),1.d0)
          xi(i,j)=max(xi(i,j),-1.d0)
          et(i,j)=min(et(i,j),1.d0)
          et(i,j)=max(et(i,j),-1.d0)
          call distattach(xi(i,j),et(i,j),pneigh,pnode,a(i,j),p,
     &       ratio,nterms)
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
!       exit if minimum found
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
     &         (j+jmin.ge.-1).and.(j+jmin.le.1)) then
              a(i,j)=aold(i+imin,j+jmin)
              xi(i,j)=xiold(i+imin,j+jmin)
              et(i,j)=etold(i+imin,j+jmin)
            else
              xi(i,j)=xi(i,j)+imin*d3
              et(i,j)=et(i,j)+jmin*d3
!
              xi(i,j)=min(xi(i,j),1.d0)
              xi(i,j)=max(xi(i,j),-1.d0)
              et(i,j)=min(et(i,j),1.d0)
              et(i,j)=max(et(i,j),-1.d0)
!
              call distattach(xi(i,j),et(i,j),pneigh,
     &          pnode,a(i,j),p,ratio,nterms)
!              write(*,*) a(i,j)
            endif
          enddo
        enddo
      enddo
!
!     4th run
!     initialisation
!
      xiopt=xi(0,0)
      etopt=et(0,0)
      do i=-1,1
        do j=-1,1
          xi(i,j)=xiopt+i*d4
          et(i,j)=etopt+j*d4
          xi(i,j)=min(xi(i,j),1.d0)
          xi(i,j)=max(xi(i,j),-1.d0)
          et(i,j)=min(et(i,j),1.d0)
          et(i,j)=max(et(i,j),-1.d0)
          call distattach(xi(i,j),et(i,j),pneigh,pnode,a(i,j),p,
     &        ratio,nterms)
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
!       exit if minimum found
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
     &         (j+jmin.ge.-1).and.(j+jmin.le.1)) then
              a(i,j)=aold(i+imin,j+jmin)
              xi(i,j)=xiold(i+imin,j+jmin)
              et(i,j)=etold(i+imin,j+jmin)
            else
              xi(i,j)=xi(i,j)+imin*d4
              et(i,j)=et(i,j)+jmin*d4
!
              xi(i,j)=min(xi(i,j),1.d0)
              xi(i,j)=max(xi(i,j),-1.d0)
              et(i,j)=min(et(i,j),1.d0)
              et(i,j)=max(et(i,j),-1.d0)
!
              call distattach(xi(i,j),et(i,j),pneigh,
     &          pnode,a(i,j),p,ratio,nterms)
!              write(*,*) a(i,j)
            endif
          enddo
        enddo
      enddo
!
      call distattach(xi(0,0),et(0,0),pneigh,pnode,a(0,0),p,
     &     ratio,nterms)
!
      do i=1,3
        pnode(i)=p(i)
      enddo
!
      dist=dsqrt(a(0,0))
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
      elseif((nterms.eq.6).or.(nterms.eq.7)) then
         xil=(xi(0,0)+1.d0)/2.d0
         etl=(et(0,0)+1.d0)/2.d0
         if(xil+etl.gt.1.d0) then
            dummy=xil
            xil=1.d0-etl
            etl=1.d0-dummy
         endif
      elseif((nterms.eq.8).or.(nterms.eq.9)) then
         xil=xi(0,0)
         etl=et(0,0)
      endif
!     
      return
      end
      
