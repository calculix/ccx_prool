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
      subroutine createlocalsys(nnfront,istartfront,iendfront,ifront,
     &     co,xt,xn,xa,nfront,ifrontrel,stress,iedno,ibounedg,ieled,
     &     kontri,isubsurffront,istartcrackfro,iendcrackfro,ncrack,
     &     angle,nstep,ier)
!     
!     sorting the crack boundary nodes according to adjacency
!     
      implicit none
!     
      integer nnfront,i,j,k,j1,j2,istartfront(*),iendfront(*),ifront(*),
     &     node,nodenext,nodelast,noderel,n,ier,iedge,iedgerel,ielem,
     &     matz,nfront,ifrontrel(*),iedno(2,*),ibounedg(*),ieled(2,*),
     &     kontri(3,*),isubsurffront(*),istartcrackfro(*),ncrack,
     &     iendcrackfro(*),nstep,m
!     
      real*8 co(3,*),xt(3,*),xn(3,nstep,*),xa(3,nstep,*),dd,s(3,3),w(3),
     &     fv1(3),fv2(3),cg(3),stress(6,nstep,*),xtj2m1(3),angle(*),
     &     z(3,3)
!     
      do i=1,nnfront
        j1=istartfront(i)
        j2=iendfront(i)
!     
!     calculate the tangent vector t
!     at the end nodes: tangent to the adjacent front edge        
!     in all other nodes: mean of the tangent to the adjacent front edges
!     
        do j=j1,j2
          node=ifront(j)
          if(j.lt.j2) then
            nodenext=ifront(j+1)
            do k=1,3
              xt(k,j)=co(k,nodenext)-co(k,node)
            enddo
          elseif(j.eq.j2) then
            if(isubsurffront(i).eq.1) then
              nodenext=ifront(j1)
              do k=1,3
                xt(k,j)=co(k,nodenext)-co(k,node)
              enddo
            else
              nodelast=ifront(j-1)
              do k=1,3
                xt(k,j)=co(k,node)-co(k,nodelast)
              enddo
            endif
          endif
!     
!     normalizing
!     
          dd=dsqrt(xt(1,j)*xt(1,j)+xt(2,j)*xt(2,j)+xt(3,j)*xt(3,j))
          do k=1,3
            xt(k,j)=xt(k,j)/dd
          enddo
        enddo
!     
!     taking the mean and normalizing (due to normalizing the factor of
!     2 is not important)
!     
        if(isubsurffront(i).eq.1) then
          do k=1,3
            xtj2m1(k)=xt(k,j2-1)
          enddo
        endif
!     
        do j=j2-1,j1+1,-1
          do k=1,3
            xt(k,j)=(xt(k,j-1)+xt(k,j))
          enddo
          dd=dsqrt(xt(1,j)*xt(1,j)+xt(2,j)*xt(2,j)+xt(3,j)*xt(3,j))
          do k=1,3
            xt(k,j)=xt(k,j)/dd
          enddo
        enddo
!     
!     for subsurface cracks: adjust the tangent at the starting and
!     end node
!     
        if(isubsurffront(i).eq.1) then
          do k=1,3
            xt(k,j1)=xt(k,j1)+xt(k,j2)
          enddo
          dd=dsqrt(xt(1,j1)*xt(1,j1)+xt(2,j1)*xt(2,j1)+
     &         xt(3,j1)*xt(3,j1))
          do k=1,3
            xt(k,j1)=xt(k,j1)/dd
          enddo
          do k=1,3
            xt(k,j2)=xt(k,j2)+xtj2m1(k)
          enddo
          dd=dsqrt(xt(1,j2)*xt(1,j2)+xt(2,j2)*xt(2,j2)+
     &         xt(3,j2)*xt(3,j2))
          do k=1,3
            xt(k,j2)=xt(k,j2)/dd
          enddo
        else
!     
!     for surface cracks: calculate the angle between the tangents
!     at the crossing points of the crack fronts with the free
!     surface; can be used to determine an appropriate shape factor
!     in shapefactor.f
!     
          angle(i)=xt(1,j1)*xt(1,j2)+xt(2,j1)*xt(2,j2)+xt(3,j1)*xt(3,j2)
          angle(i)=dacos(angle(i))
        endif
!     
!     calculating the normal vector n to the largest principal stress
!     plane
!     
!     stress order in cgx: xx,yy,zz,xy,yz,xz
!     
        do j=j1,j2
          noderel=ifrontrel(j)
          do m=1,nstep
            s(1,1)=stress(1,m,noderel)
            s(1,2)=stress(4,m,noderel)
            s(1,3)=stress(6,m,noderel)
            s(2,1)=s(1,2)
            s(2,2)=stress(2,m,noderel)
            s(2,3)=stress(5,m,noderel)
            s(3,1)=s(1,3)
            s(3,2)=s(2,3)
            s(3,3)=stress(3,m,noderel)
!     
!     determining the eigenvalues
!     
            n=3
            matz=1
            ier=0
            call rs(n,n,s,w,matz,z,fv1,fv2,ier)
            if(ier.ne.0) then
              write(*,*)
     &             '*ERROR in createlocalsys while calculating the'
              write(*,*) '       eigenvalues/eigenvectors'
              ier=1
c              call exit(201)
            endif
!     
!     normal to largest principal stress plane
!     
            do k=1,3
              xn(k,m,j)=z(k,3)
            enddo
!     
!     projection on a plane orthogonal to the local tangent vector
!     (necessary since factory roofing effect is not taken into account)
!     
            dd=xn(1,m,j)*xt(1,j)+xn(2,m,j)*xt(2,j)+xn(3,m,j)*xt(3,j)
            do k=1,3
              xn(k,m,j)=xn(k,m,j)-dd*xt(k,j)
            enddo
!     
!     normalizing vector xn
!     
            dd=dsqrt(xn(1,m,j)*xn(1,m,j)+xn(2,m,j)*xn(2,m,j)
     &                                  +xn(3,m,j)*xn(3,m,j))
            do k=1,3
              xn(k,m,j)=xn(k,m,j)/dd
            enddo
          enddo
        enddo
!     
!     propagation direction a=+t x n or a=-t x n 
!     
        do j=j1,j2
          node=ifront(j)
          do m=1,nstep
            xa(1,m,j)=xt(2,j)*xn(3,m,j)-xt(3,j)*xn(2,m,j)
            xa(2,m,j)=xt(3,j)*xn(1,m,j)-xt(1,j)*xn(3,m,j)
            xa(3,m,j)=xt(1,j)*xn(2,m,j)-xt(2,j)*xn(1,m,j)
            dd=-xa(1,m,j)*co(1,node)-xa(2,m,j)*co(2,node)
     &                              -xa(3,m,j)*co(3,node)
!     
!     the equation of the plane orthogonal to xa and through
!     node is xa(1,j)*x+xa(2,j)*y+xa(3,j)*z+dd=0
!     
!     the direction of xa is away from the crack plane is the center
!     of gravity cg of an adjacent triangle is such that
!     xa(1,j)*cg(1)+xa(2,j)*cg(2)+xa(3,j)*cg(3)+dd<0
!     
            noderel=ifrontrel(j)
            iedgerel=iedno(1,noderel)
            iedge=ibounedg(iedgerel)
            ielem=ieled(1,iedge)
            do k=1,3
              cg(k)=(co(k,kontri(1,ielem))+
     &             co(k,kontri(2,ielem))+
     &             co(k,kontri(3,ielem)))/3.d0
            enddo
!     
            if(xa(1,m,j)*cg(1)+xa(2,m,j)*cg(2)+xa(3,m,j)*cg(3)+dd.gt.0)
     &           then
            do k=1,3
              xa(k,m,j)=-xa(k,m,j)
            enddo
          endif
        enddo
      enddo
      enddo
!     
c     do i=1,nfront
c     write(*,*) 'createlocalsys xt'
c     write(*,*) 'xt ',i,xt(1,i),xt(2,i),xt(3,i)
c     enddo
c     do i=1,nfront
c     write(*,*) 'createlocalsys xn'
c     write(*,*) 'xn ',i,xn(1,i),xn(2,i),xn(3,i)
c     enddo
c     do i=1,nfront
c     write(*,*) 'createlocalsys xa'
c     write(*,*) 'xa ',i,xa(1,i),xa(2,i),xa(3,i)
c     enddo
c     write(*,*)
!     
      return
      end

