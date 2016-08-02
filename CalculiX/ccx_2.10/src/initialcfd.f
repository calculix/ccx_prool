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
      subroutine initialcfd(nef,ipkonf,kon,lakonf,co,coel,cofa,nface,
     &  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa,
     &  volume,neifa,xxj,cosb,vel,dmin)
!
!     calculating geometric variables of the cells and their faces
!
      implicit none
!
      character*8 lakonf(*)
!
      integer nef,ipkonf(*),kon(*),nface,ielfa(4,*),ipnei(*),neiel(*),
     &  ifaceq(8,6),i,j,k,indexe,kflag,index1,index2,j1,j2,nope,
     &  nodes(4),iel1,iel2,iel3,iface,indexf,neifa(*),nf(5),ifacet(7,4),
     &  ifacew(8,5),numfaces,ied4(2,6),ied6(2,9),ied8(2,12)
!
      real*8 co(3,*),coel(3,*),cofa(3,*),area(*),xxn(3,*),xxi(3,*),
     &  xle(*),xlen(*),xlet(*),xrlfa(3,*),cosa(*),xsj2(3),xi,et,
     &  shp2(7,4),xs2(3,7),xl2(3,8),xl13,volume(*),dxsj2,xl(3,8),
     &  xxj(3,*),cosb(*),dmin,vel(nef,0:5)
!
!     nodes belonging to the cell faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data nf /3,3,4,4,4/
      data ied4 /1,2,2,3,3,1,1,4,2,4,3,4/
      data ied6 /1,2,2,3,3,1,4,5,5,6,6,4,1,4,2,5,3,6/
      data ied8 /1,2,2,3,3,4,4,1,5,6,6,7,7,8,8,1,1,5,2,6,3,7,4,8/
!
!     coordinates of the center of the cells
!
      do i=1,nef
         if(ipkonf(i).lt.0) cycle
         if(lakonf(i)(1:1).ne.'F') cycle
         indexe=ipkonf(i)
         if(lakonf(i)(4:4).eq.'8') then
            nope=8
         else if(lakonf(i)(4:4).eq.'6') then
            nope=6
         else
            nope=4
         endif
         do j=1,3
            do k=1,nope
               coel(j,i)=coel(j,i)+co(j,kon(indexe+k))
            enddo
            coel(j,i)=coel(j,i)/nope
         enddo
      enddo
!
      kflag=2
!
!     loop over all faces
!
      do i=1,nface
         iel1=ielfa(1,i)
         indexe=ipkonf(iel1)
         j1=ielfa(4,i)
         if(lakonf(iel1)(4:4).eq.'8') then
!
!           hexahedral element
!
!           coordinates of the face centers
!
            do j=1,4
               nodes(j)=kon(indexe+ifaceq(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/4.d0
            enddo
!
            xi=0.d0
            et=0.d0
            call shape4q(xi,et,xl2,xsj2,xs2,shp2,kflag)
!
!           area of the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
            area(i)=4.d0*dxsj2
!
!           normal and xi-vector on face viewed from cell 1
!
            index1=ipnei(iel1)+j1
            do k=1,3
               xxn(k,index1)=xsj2(k)/dxsj2
               xxi(k,index1)=cofa(k,i)-coel(k,iel1)
            enddo
!
!           distance from face center to the center of cell 1
!
            xle(index1)=dsqrt(xxi(1,index1)**2+xxi(2,index1)**2+
     &                         xxi(3,index1)**2)
            do k=1,3
               xxi(k,index1)=xxi(k,index1)/xle(index1)
            enddo
!
!           angle between the normal and the xi-vector
!
            cosa(index1)=xxn(1,index1)*xxi(1,index1)+
     &                   xxn(2,index1)*xxi(2,index1)+
     &                   xxn(3,index1)*xxi(3,index1)
!
            iel2=ielfa(2,i)
!
!           check whether there is an adjacent cell
!
            if(iel2.ne.0) then
               index2=ipnei(iel2)
               do j2=1,6
                  index2=index2+1
                  if(neiel(index2).eq.iel1) exit
               enddo
!
!              normal and xi-vector on face viewed from cell 2
!
               do k=1,3
                  xxi(k,index2)=cofa(k,i)-coel(k,iel2)
                  xxn(k,index2)=-xxn(k,index1)
               enddo
               xle(index2)=dsqrt(xxi(1,index2)**2+xxi(2,index2)**2+
     &              xxi(3,index2)**2)
               do k=1,3
                  xxi(k,index2)=xxi(k,index2)/xle(index2)
               enddo
!
!              angle between the normal and the xi-vector: xxn.xxi
!
               cosa(index2)=xxn(1,index2)*xxi(1,index2)+
     &                      xxn(2,index2)*xxi(2,index2)+
     &                      xxn(3,index2)*xxi(3,index2)
!
!              distance from the face center to the center of the
!              adjacent cell
!
               xlen(index1)=xle(index2)
               xlen(index2)=xle(index1)
!
               do k=1,3
                  xxj(k,index2)=coel(k,iel1)-coel(k,iel2)
               enddo
!
!              distance between the cell center and the center of the
!              adjacent cell
!
               xlet(index1)=dsqrt(xxj(1,index2)**2+xxj(2,index2)**2
     &                                            +xxj(3,index2)**2)
               xlet(index2)=xlet(index1)
!
!              xxj is the unit vector connecting neighboring cell centers
!
               do k=1,3
                  xxj(k,index2)=xxj(k,index2)/xlet(index2)
                  xxj(k,index1)=-xxj(k,index2)
               enddo
!
!              xxn.xxj
!
               cosb(index2)=xxn(1,index1)*xxj(1,index1)+
     &                      xxn(2,index1)*xxj(2,index1)+
     &                      xxn(3,index1)*xxj(3,index1)
               cosb(index1)=cosb(index2)
!
               xrlfa(1,i)=xle(index2)/(xle(index1)+xle(index2))
               xrlfa(2,i)=xle(index1)/(xle(index1)+xle(index2))
            else
!
!              xxi and xxj coincide
!
               do k=1,3
                  xxj(k,index1)=xxi(k,index1)
               enddo
               cosb(index1)=cosa(index1)
!
!              external face: determining the cell next to the
!              adjacent cell
!
               iel3=ielfa(3,i)
               if(iel3.eq.0) cycle
               xl13=dsqrt((coel(1,iel1)-coel(1,iel3))**2+
     &                    (coel(2,iel1)-coel(2,iel3))**2+
     &                    (coel(3,iel1)-coel(3,iel3))**2)
               xrlfa(1,i)=(xl13+xle(index1))/xl13
               xrlfa(3,i)=1.d0-xrlfa(1,i)
            endif
         else if(lakonf(iel1)(4:4).eq.'6') then
!
!           wedge element
!
!           coordinates of the face centers
!
            do j=1,nf(j1)
               nodes(j)=kon(indexe+ifacew(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/nf(j1)
            enddo
!
            xi=0.d0
            et=0.d0
            if(nf(j1).eq.3) then
               call shape3tri(xi,et,xl2,xsj2,xs2,shp2,kflag)
            else
               call shape4q(xi,et,xl2,xsj2,xs2,shp2,kflag)
            endif
!
!           area of the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
            if(nf(j1).eq.3) then
               area(i)=dxsj2/2.d0
            else
               area(i)=4.d0*dxsj2
            endif
!
!           normal and xi-vector on face viewed from cell 1
!
            index1=ipnei(iel1)+j1
            do k=1,3
               xxn(k,index1)=xsj2(k)/dxsj2
               xxi(k,index1)=cofa(k,i)-coel(k,iel1)
            enddo
!
!           distance from face center to the center of cell 1
!
            xle(index1)=dsqrt(xxi(1,index1)**2+xxi(2,index1)**2+
     &                         xxi(3,index1)**2)
            do k=1,3
               xxi(k,index1)=xxi(k,index1)/xle(index1)
            enddo
!
!           angle between the normal and the xi-vector
!
            cosa(index1)=xxn(1,index1)*xxi(1,index1)+
     &                   xxn(2,index1)*xxi(2,index1)+
     &                   xxn(3,index1)*xxi(3,index1)
!
            iel2=ielfa(2,i)
!
!           check whether there is an adjacent cell
!
            if(iel2.ne.0) then
               index2=ipnei(iel2)
               do j2=1,5
                  index2=index2+1
                  if(neiel(index2).eq.iel1) exit
               enddo
!
!              normal and xi-vector on face viewed from cell 2
!
               do k=1,3
                  xxi(k,index2)=cofa(k,i)-coel(k,iel2)
                  xxn(k,index2)=-xxn(k,index1)
               enddo
               xle(index2)=dsqrt(xxi(1,index2)**2+xxi(2,index2)**2+
     &              xxi(3,index2)**2)
               do k=1,3
                  xxi(k,index2)=xxi(k,index2)/xle(index2)
               enddo
!
!              angle between the normal and the xi-vector: xxn.xxi
!
               cosa(index2)=xxn(1,index2)*xxi(1,index2)+
     &                      xxn(2,index2)*xxi(2,index2)+
     &                      xxn(3,index2)*xxi(3,index2)
!
!              distance from the face center to the center of the
!              adjacent cell
!
               xlen(index1)=xle(index2)
               xlen(index2)=xle(index1)
!
               do k=1,3
                  xxj(k,index2)=coel(k,iel1)-coel(k,iel2)
               enddo
!
!              distance between the cell center and the center of the
!              adjacent cell
!
               xlet(index1)=dsqrt(xxj(1,index2)**2+xxj(2,index2)**2
     &                                            +xxj(3,index2)**2)
               xlet(index2)=xlet(index1)
!
!              xxj is the unit vector connecting neighboring cell centers
!
               do k=1,3
                  xxj(k,index2)=xxj(k,index2)/xlet(index2)
                  xxj(k,index1)=-xxj(k,index2)
               enddo
!
!              xxn.xxj
!
               cosb(index2)=xxn(1,index1)*xxj(1,index1)+
     &                      xxn(2,index1)*xxj(2,index1)+
     &                      xxn(3,index1)*xxj(3,index1)
               cosb(index1)=cosb(index2)
!
               xrlfa(1,i)=xle(index2)/(xle(index1)+xle(index2))
               xrlfa(2,i)=xle(index1)/(xle(index1)+xle(index2))
            else
!
!              xxi and xxj coincide
!
               do k=1,3
                  xxj(k,index1)=xxi(k,index1)
               enddo
               cosb(index1)=cosa(index1)
!
!              external face: determining the cell next to the
!              adjacent cell
!
               iel3=ielfa(3,i)
               if(iel3.eq.0) cycle
               xl13=dsqrt((coel(1,iel1)-coel(1,iel3))**2+
     &                    (coel(2,iel1)-coel(2,iel3))**2+
     &                    (coel(3,iel1)-coel(3,iel3))**2)
               xrlfa(1,i)=(xl13+xle(index1))/xl13
               xrlfa(3,i)=1.d0-xrlfa(1,i)
            endif
         else
!
!           tetrahedral element
!
!           coordinates of the face centers
!
            do j=1,3
               nodes(j)=kon(indexe+ifacet(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
                  cofa(k,i)=cofa(k,i)+xl2(k,j)
               enddo
            enddo
            do k=1,3
               cofa(k,i)=cofa(k,i)/3.d0
            enddo
!
            xi=0.d0
            et=0.d0
            call shape3tri(xi,et,xl2,xsj2,xs2,shp2,kflag)
!
!           area of the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &                    xsj2(3)*xsj2(3))
            area(i)=dxsj2/2.d0
!
!           normal and xi-vector on face viewed from cell 1
!
            index1=ipnei(iel1)+j1
            do k=1,3
               xxn(k,index1)=xsj2(k)/dxsj2
               xxi(k,index1)=cofa(k,i)-coel(k,iel1)
            enddo
!
!           distance from face center to the center of cell 1
!
            xle(index1)=dsqrt(xxi(1,index1)**2+xxi(2,index1)**2+
     &                         xxi(3,index1)**2)
            do k=1,3
               xxi(k,index1)=xxi(k,index1)/xle(index1)
            enddo
!
!           angle between the normal and the xi-vector
!
            cosa(index1)=xxn(1,index1)*xxi(1,index1)+
     &                   xxn(2,index1)*xxi(2,index1)+
     &                   xxn(3,index1)*xxi(3,index1)
!
            iel2=ielfa(2,i)
!
!           check whether there is an adjacent cell
!
            if(iel2.ne.0) then
               index2=ipnei(iel2)
               do j2=1,4
                  index2=index2+1
                  if(neiel(index2).eq.iel1) exit
               enddo
!
!              normal and xi-vector on face viewed from cell 2
!
               do k=1,3
                  xxi(k,index2)=cofa(k,i)-coel(k,iel2)
                  xxn(k,index2)=-xxn(k,index1)
               enddo
               xle(index2)=dsqrt(xxi(1,index2)**2+xxi(2,index2)**2+
     &              xxi(3,index2)**2)
               do k=1,3
                  xxi(k,index2)=xxi(k,index2)/xle(index2)
               enddo
!
!              angle between the normal and the xi-vector: xxn.xxi
!
               cosa(index2)=xxn(1,index2)*xxi(1,index2)+
     &                      xxn(2,index2)*xxi(2,index2)+
     &                      xxn(3,index2)*xxi(3,index2)
!
!              distance from the face center to the center of the
!              adjacent cell
!
               xlen(index1)=xle(index2)
               xlen(index2)=xle(index1)
!
               do k=1,3
                  xxj(k,index2)=coel(k,iel1)-coel(k,iel2)
               enddo
!
!              distance between the cell center and the center of the
!              adjacent cell
!
               xlet(index1)=dsqrt(xxj(1,index2)**2+xxj(2,index2)**2
     &                                            +xxj(3,index2)**2)
               xlet(index2)=xlet(index1)
!
!              xxj is the unit vector connecting neighboring cell centers
!
               do k=1,3
                  xxj(k,index2)=xxj(k,index2)/xlet(index2)
                  xxj(k,index1)=-xxj(k,index2)
               enddo
!
!              xxn.xxj
!
               cosb(index2)=xxn(1,index1)*xxj(1,index1)+
     &                      xxn(2,index1)*xxj(2,index1)+
     &                      xxn(3,index1)*xxj(3,index1)
               cosb(index1)=cosb(index2)
!
               xrlfa(1,i)=xle(index2)/(xle(index1)+xle(index2))
               xrlfa(2,i)=xle(index1)/(xle(index1)+xle(index2))
            else
!
!              xxi and xxj coincide
!
               do k=1,3
                  xxj(k,index1)=xxi(k,index1)
               enddo
               cosb(index1)=cosa(index1)
!
!              external face: determining the cell next to the
!              adjacent cell
!
               iel3=ielfa(3,i)
               if(iel3.eq.0) cycle
               xl13=dsqrt((coel(1,iel1)-coel(1,iel3))**2+
     &                    (coel(2,iel1)-coel(2,iel3))**2+
     &                    (coel(3,iel1)-coel(3,iel3))**2)
               xrlfa(1,i)=(xl13+xle(index1))/xl13
               xrlfa(3,i)=1.d0-xrlfa(1,i)
            endif
         endif
c         write(*,*)
      enddo
!
!     calculation of the volume of the elements
!
      do i=1,nef
         if(ipkonf(i).lt.0) cycle
         if(lakonf(i)(1:1).ne.'F') cycle
         indexf=ipnei(i)
         volume(i)=0.d0
         if(lakonf(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakonf(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
            iface=neifa(indexf+j)
            volume(i)=volume(i)+
     &            area(iface)*cofa(1,iface)*xxn(1,indexf+j)
         enddo
c         write(*,*) 'initialcfd volume ',i,volume(i)
      enddo
!     
!     calculation of the minimum length within the cells
!
      dmin=1.d30
      do i=1,nef
         indexe=ipkonf(i)
         read(lakonf(i)(4:4),'(i1)') nope
         do j=1,nope
            do k=1,3
               xl(k,j)=co(k,kon(indexe+j))
            enddo
         enddo
         if(nope.eq.4) then
            do j=1,6
               dmin=min(dmin,(xl(1,ied4(1,j))-xl(1,ied4(2,j)))**2+
     &                       (xl(2,ied4(1,j))-xl(2,ied4(2,j)))**2+
     &                       (xl(3,ied4(1,j))-xl(3,ied4(2,j)))**2)
            enddo
         elseif(nope.eq.6) then
            do j=1,9
               dmin=min(dmin,(xl(1,ied6(1,j))-xl(1,ied6(2,j)))**2+
     &                       (xl(2,ied6(1,j))-xl(2,ied6(2,j)))**2+
     &                       (xl(3,ied6(1,j))-xl(3,ied6(2,j)))**2)
            enddo
         else
            do j=1,12
               dmin=min(dmin,(xl(1,ied8(1,j))-xl(1,ied8(2,j)))**2+
     &                       (xl(2,ied8(1,j))-xl(2,ied8(2,j)))**2+
     &                       (xl(3,ied8(1,j))-xl(3,ied8(2,j)))**2)
            enddo
         endif
      enddo
      dmin=dsqrt(dmin)
!
      return
      end
