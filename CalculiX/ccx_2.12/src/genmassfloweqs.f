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
      subroutine genmassfloweqs(ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,co,ikmpc,ilmpc,labmpc,
     &  lakon,nload,sideload,ipkon,kon,nelemload)
!
!     generations fluid equations for prescribed zero mass
!     flow in CFD-calculations
!
      implicit none
!
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,
     &  i,j,number,indexe,iel1,nodes(4),
     &  mpcfreeold,ikmpc(*),ilmpc(*),id,idof,nf(5),kflag,j1,
     &  k,iface,ifaceq(8,6),ifacet(6,4),ifacew(8,5),
     &  ipkon(*),nload,nelemload(2,*),kon(*)
!
      real*8 coefmpc(*),co(3,*),xl2(3,8),
     &  xn(3),xsj2(3),dxsj2,shp2(7,4),xs2(3,7),xi,et
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,
     &             1,2,4,5,9,8,
     &             2,3,4,6,10,9,
     &             1,4,3,8,10,7/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data nf /3,3,4,4,4/
!
      kflag=2
!
      do i=1,nload
         if(sideload(i)(1:1).ne.'M') cycle
         iel1=nelemload(1,i)
         read(sideload(i)(2:2),'(i1)') j1
         iface=10*iel1+j1
!
         nmpc=nmpc+1
         if(nmpc.gt.nmpc_) then
            write(*,*) '*ERROR in genmassfloweqs: increase nmpc_'
            call exit(201)
         endif
!
         labmpc(nmpc)='FLUID               '
         ipompc(nmpc)=mpcfree
!
         indexe=ipkon(iel1)
         if(lakon(iel1)(4:4).eq.'8') then
!
!           hexahedral element
!     
            do j=1,4
               nodes(j)=kon(indexe+ifaceq(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
               enddo
            enddo
!     
            xi=0.d0
            et=0.d0
            call shape4q(xi,et,xl2,xsj2,xs2,shp2,kflag)
!     
!     normal to the face
!     
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
            do k=1,3
               xn(k)=xsj2(k)/dxsj2
            enddo
         else if(lakon(iel1)(4:4).eq.'6') then
!     
!           wedge element
!
            do j=1,nf(j1)
               nodes(j)=kon(indexe+ifacew(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
               enddo
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
!           normal to the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
            do k=1,3
               xn(k)=xsj2(k)/dxsj2
            enddo
         else
!
!           tetrahedral element
!
            do j=1,3
               nodes(j)=kon(indexe+ifacet(j,j1))
               do k=1,3
                  xl2(k,j)=co(k,nodes(j))
               enddo
            enddo
!
            xi=0.d0
            et=0.d0
            call shape3tri(xi,et,xl2,xsj2,xs2,shp2,kflag)
!     
!           normal to the face
!
            dxsj2=dsqrt(xsj2(1)*xsj2(1)+xsj2(2)*xsj2(2)+
     &           xsj2(3)*xsj2(3))
            do k=1,3
               xn(k)=xsj2(k)/dxsj2
            enddo
         endif
!
!
!     determining which direction to use for the
!     dependent side: should not occur on the dependent
!     side in another MPC and should have a nonzero
!     coefficient
! 
         number=0
         do j=1,3
            number=number+1
c            if(number.gt.3) number=1
            idof=-(8*(iface-1)+number)
            call nident(ikmpc,idof,nmpc-1,id)
            if(id.gt.0) then
               if(ikmpc(id).eq.idof) then
                  cycle
               endif
            endif
            if(dabs(xn(number)).lt.1.d-5) cycle
            exit
         enddo
         if(j.gt.3) then
            write(*,*) 
     &           '*ERROR reading *EQUATIONF: MPC on face'
            write(*,*) iface,' of element',iel1
            write(*,*) ' in transformed coordinates'
            write(*,*) ' cannot be converted in MPC: all'
            write(*,*) ' DOFs in the node are used as'
            write(*,*) ' dependent nodes in other MPCs'
            call exit(201)
         endif
         number=number-1
!     
!     updating ikmpc and ilmpc
!     
         do j=nmpc,id+2,-1
            ikmpc(j)=ikmpc(j-1)
            ilmpc(j)=ilmpc(j-1)
         enddo
         ikmpc(id+1)=idof
         ilmpc(id+1)=nmpc
!     
         do j=1,3
            number=number+1
            if(number.gt.3) number=1
            if(dabs(xn(number)).lt.1.d-5) cycle
            nodempc(1,mpcfree)=iface
            nodempc(2,mpcfree)=number
            coefmpc(mpcfree)=xn(number)
            mpcfreeold=mpcfree
            mpcfree=nodempc(3,mpcfree)
            if(mpcfree.eq.0) then
               write(*,*) 
     &              '*ERROR in genmassfloweqs: increase nmpc_'
               call exit(201)
            endif
         enddo
         nodempc(3,mpcfreeold)=0
      enddo
!
      return
      end

