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
      subroutine crackplane(ncrack,istartcrackfro,iendcrackfro,ifront,
     &     co,xa,xplanecrack,istartcrackbou,iendcrackbou,
     &     ibounnod,xn,xt,ifrontrel,iedno,ibounedg,ieled,kontri,
     &     nfront,costruc,xnplane,xaplane,nstep,xm,xb)
!     
!     determine for each front a plane approximating mode-I
!     crack propagation
!     
      implicit none
!     
      integer i,j,k,l,istartcrackfro(*),iendcrackfro(*),nodestart,
     &     nodeend,ifront(*),node,ncrack,istartcrackbou(*),noderel,
     &     iendcrackbou(*),ibounnod(*),iedge,iedgerel,ielem,nfront,
     &     ifrontrel(*),iedno(2,*),ibounedg(*),ieled(2,*),kontri(3,*),
     &     n1,n2,n3,nstep,m
!     
      real*8 a(3,3),b(3),co(3,*),xlength,x,y,z,xa(3,nstep,*),det,
     &     dd,xplanecrack(4,nstep,*),xn(3,nstep,*),xt(3,*),cg(3),
     &     costruc(3,*),p12(3),p23(3),xn1(3),xn2(3),xnplane(3,nstep,*),
     &     xaplane(3,nstep,*),ca,cb,cc,xm(3,*),xb(3,*)
!     
!     loop over all cracks (each crack has a unique mean plane,
!     no matter how many fronts belong to the crack)
!     
      do i=1,ncrack
!     
!     setting the LHS (a; only upper triangle is used since the matrix
!     is symmetric) and the RHS (b) to zero
!     
        do m=1,nstep
!     
          do k=1,3
            b(k)=0.d0
            do l=k,3
              a(k,l)=0.d0
            enddo
          enddo
!     
!     calculate the "length" of the crack front
!     
          nodestart=ifront(istartcrackfro(i))
          nodeend=ifront(iendcrackfro(i))
          xlength=dsqrt((co(1,nodestart)-co(1,nodeend))**2+
     &         (co(2,nodestart)-co(2,nodeend))**2+
     &         (co(3,nodestart)-co(3,nodeend))**2)
!     
!     loop over all nodes belonging to the front
!     
          do j=istartcrackfro(i),iendcrackfro(i)
!     
!     taking into account the front node
!     
            node=ifront(j)
            x=co(1,node)
            y=co(2,node)
            z=co(3,node)
            a(1,1)=a(1,1)+x*x
            a(1,2)=a(1,2)+x*y
            a(1,3)=a(1,3)+x*z
            a(2,2)=a(2,2)+y*y
            a(2,3)=a(2,3)+y*z
            a(3,3)=a(3,3)+z*z
            b(1)=b(1)+x
            b(2)=b(2)+y
            b(3)=b(3)+z
!     
!     taking into account the fictitiously propagated front node
!     along the largest principal stress plane (mode-I)
!     
            x=x+xlength*xa(1,m,j)
            y=y+xlength*xa(2,m,j)
            z=z+xlength*xa(3,m,j)
            a(1,1)=a(1,1)+x*x
            a(1,2)=a(1,2)+x*y
            a(1,3)=a(1,3)+x*z
            a(2,2)=a(2,2)+y*y
            a(2,3)=a(2,3)+y*z
            a(3,3)=a(3,3)+z*z
            b(1)=b(1)+x
            b(2)=b(2)+y
            b(3)=b(3)+z
          enddo
!     
!     solving the system of equations (3 x 3 system)
!     
          det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(2,3))
     &         -a(1,2)*(a(1,2)*a(3,3)-a(1,3)*a(2,3))
     &         +a(1,3)*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
          ca=b(1)*(a(2,2)*a(3,3)-a(2,3)*a(2,3))
     &         -b(2)*(a(1,2)*a(3,3)-a(2,3)*a(1,3))
     &         +b(3)*(a(1,2)*a(2,3)-a(2,2)*a(1,3))
          cb=-b(1)*(a(1,2)*a(3,3)-a(1,3)*a(2,3))
     &         +b(2)*(a(1,1)*a(3,3)-a(1,3)*a(1,3))
     &         -b(3)*(a(1,1)*a(2,3)-a(1,2)*a(1,3))
          cc=b(1)*(a(1,2)*a(2,3)-a(1,3)*a(2,2))
     &         -b(2)*(a(1,1)*a(2,3)-a(1,3)*a(1,2))
     &         +b(3)*(a(1,1)*a(2,2)-a(1,2)*a(1,2))
          ca=ca/det
          cb=cb/det
          cc=cc/det
!     
!     scaling the equation of the plane such that ca*ca+cb*cb+cc*cc=1
!     
          dd=dsqrt(ca*ca+cb*cb+cc*cc)
          xplanecrack(1,m,i)=ca/dd
          xplanecrack(2,m,i)=cb/dd
          xplanecrack(3,m,i)=cc/dd
          xplanecrack(4,m,i)=1.d0/dd
!     
!     note: the direction of vector xplanecrack is hereto unrelated
!     to the direction of the tangential vectors along the crack front
!     
c          write(*,*)
c          write(*,*) 'crackplane'
c          write(*,*) 'equation of the crack plane for crack ',i
c          write(*,*) xplanecrack(1,m,i),' * x + ',xplanecrack(2,m,i),
c     &         ' * y + ',xplanecrack(3,m,i),' * z = ',xplanecrack(4,m,i)
c          write(*,*)
!     
        enddo
      enddo
!     
!     creating a new local system based on
!     - the tangential vector t
!     - the normal vector on the adjacent triangles of the crack mesh
!     - a = t x n
!     the direction of n corresponds according to the corkscrew rule
!     with the node numbering of the crack elements      
!     
      do i=1,ncrack
        do j=istartcrackfro(i),iendcrackfro(i)
!     
!     normal on the triangle belonging to the one adjacent edge
!     
          noderel=ifrontrel(j)
          iedgerel=iedno(1,noderel)
          iedge=ibounedg(iedgerel)
          ielem=ieled(1,iedge)
!     
          n1=kontri(1,ielem)
          n2=kontri(2,ielem)
          n3=kontri(3,ielem)
          do k=1,3
            p12(k)=co(k,n2)-co(k,n1)
            p23(k)=co(k,n3)-co(k,n2)
          enddo
          xn1(1)=p12(2)*p23(3)-p12(3)*p23(2)
          xn1(2)=p12(3)*p23(1)-p12(1)*p23(3)
          xn1(3)=p12(1)*p23(2)-p12(2)*p23(1)
!     
!     normal on the triangle belonging to the other adjacent edge
!     
          iedgerel=iedno(2,noderel)
          iedge=ibounedg(iedgerel)
          ielem=ieled(1,iedge)
!     
          n1=kontri(1,ielem)
          n2=kontri(2,ielem)
          n3=kontri(3,ielem)
          do k=1,3
            p12(k)=co(k,n2)-co(k,n1)
            p23(k)=co(k,n3)-co(k,n2)
          enddo
          xn2(1)=p12(2)*p23(3)-p12(3)*p23(2)
          xn2(2)=p12(3)*p23(1)-p12(1)*p23(3)
          xn2(3)=p12(1)*p23(2)-p12(2)*p23(1)
!     
!     taking the mean (factor of 2 is not important due to
!     subsequent normalization)
!     
          do k=1,3
            xm(k,j)=xn1(k)+xn2(k)
          enddo
!     
!     projection on a plane orthogonal to the local tangent vector
!     xm.xt=0 must apply
!     
          dd=xm(1,j)*xt(1,j)+xm(2,j)*xt(2,j)+xm(3,j)*xt(3,j)
          do k=1,3
            xm(k,j)=xm(k,j)-dd*xt(k,j)
          enddo
!     
!     normalizing vector xm
!     
          dd=dsqrt(xm(1,j)*xm(1,j)+xm(2,j)*xm(2,j)+xm(3,j)*xm(3,j))
          do k=1,3
            xm(k,j)=xm(k,j)/dd
          enddo
!     
!     propagation direction a=t x n 
!     
          xb(1,j)=xt(2,j)*xm(3,j)-xt(3,j)*xm(2,j)
          xb(2,j)=xt(3,j)*xm(1,j)-xt(1,j)*xm(3,j)
          xb(3,j)=xt(1,j)*xm(2,j)-xt(2,j)*xm(1,j)
        enddo
      enddo
!     
!     calculating the plane propagation direction(xaplane) based on 
!     the normal of crack plane (xnplane) and the tangential vector;
!     
!     xaplane is needed because the crack length is calculated
!     (in cracklength.f) as the distance between a front node and
!     the intersection of the vector xaplane in that node with
!     the free surface
!     
      do i=1,ncrack
        do j=istartcrackfro(i),iendcrackfro(i)
          do m=1,nstep
            xnplane(1,m,j)=xplanecrack(1,m,i)
            xnplane(2,m,j)=xplanecrack(2,m,i)
            xnplane(3,m,j)=xplanecrack(3,m,i)
!     
!     projection on a plane orthogonal to the local tangent vector
!     xnplane.xt=0 must apply
!     
            dd=xnplane(1,m,j)*xt(1,j)+xnplane(2,m,j)*xt(2,j)+
     &           xnplane(3,m,j)*xt(3,j)
            do k=1,3
              xnplane(k,m,j)=xnplane(k,m,j)-dd*xt(k,j)
            enddo
!     
!     normalizing vector xnplane
!     
            dd=dsqrt(xnplane(1,m,j)*xnplane(1,m,j)
     &           +xnplane(2,m,j)*xnplane(2,m,j)
     &           +xnplane(3,m,j)*xnplane(3,m,j))
            do k=1,3
              xnplane(k,m,j)=xnplane(k,m,j)/dd
            enddo
!     
!     propagation direction a_plane=+t x n_plane or a_plane=-t x n_plane 
!     
            xaplane(1,m,j)=xt(2,j)*xnplane(3,m,j)-xt(3,j)*xnplane(2,m,j)
            xaplane(2,m,j)=xt(3,j)*xnplane(1,m,j)-xt(1,j)*xnplane(3,m,j)
            xaplane(3,m,j)=xt(1,j)*xnplane(2,m,j)-xt(2,j)*xnplane(1,m,j)
!     
!     check if the propagation direction is correct: xnplane.xn>0 must apply
!     
            if(xnplane(1,m,j)*xm(1,j)+xnplane(2,m,j)*xm(2,j)+
     &           xnplane(3,m,j)*xm(3,j).lt.0.d0) then
!     
!     inverting the propagation direction       
!     
              do k=1,3
                xaplane(k,m,j)=-xaplane(k,m,j)
              enddo
c              write(*,*) '*WARNING: inverting the normal for plane',i,j
            endif
          enddo
        enddo
      enddo
!     
c     write(*,*) 'crackplane xaplane'
c     do i=1,nfront
c     write(*,100) 'xaplane ',i,xaplane(1,i),xaplane(2,i),xaplane(3,i)
c     enddo
c      write(*,*) 'crackplane xt'
c      do i=1,nfront
c        write(*,101) i,xt(1,i),xt(2,i),xt(3,i)
c      enddo
c      write(*,*) 'crackplane xm'
c      do i=1,nfront
c        write(*,102) i,xm(1,i),xm(2,i),xm(3,i)
c      enddo
c      write(*,*) 'crackplane xa'
c      do i=1,nfront
c        write(*,103) i,xb(1,i),xb(2,i),xb(3,i)
c      enddo
c      write(*,*)
c 100  format('xaplane',i5,3(1x,e11.4))
c 101  format('xt',i5,3(1x,e11.4))
c 102  format('xm',i5,3(1x,e11.4))
c 103  format('xa',i5,3(1x,e11.4))
!     
      return
      end

