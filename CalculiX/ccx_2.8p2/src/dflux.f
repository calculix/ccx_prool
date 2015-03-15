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
      subroutine dflux(flux,sol,kstep,kinc,time,noel,npt,coords,
     &     jltyp,temp,press,loadtype,area,vold,co,lakonl,konl,
     &     ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,iscale,mi)
!
!     user subroutine dflux
!
!
!     INPUT:
!
!     sol                current temperature value
!     kstep              step number
!     kinc               increment number
!     time(1)            current step time
!     time(2)            current total time
!     noel               element number
!     npt                integration point number
!     coords(1..3)       global coordinates of the integration point
!     jltyp              loading face kode:
!                        1  = body flux
!                        11 = face 1 
!                        12 = face 2 
!                        13 = face 3 
!                        14 = face 4 
!                        15 = face 5 
!                        16 = face 6
!     temp               currently not used
!     press              currently not used
!     loadtype           load type label
!     area               for surface flux: area covered by the
!                            integration point
!                        for body flux: volume covered by the
!                            integration point
!     vold(0..4,1..nk)   solution field in all nodes
!                        0: temperature
!                        1: displacement in global x-direction
!                        2: displacement in global y-direction
!                        3: displacement in global z-direction
!                        4: static pressure
!     co(3,1..nk)        coordinates of all nodes
!                        1: coordinate in global x-direction
!                        2: coordinate in global y-direction
!                        3: coordinate in global z-direction
!     lakonl             element label
!     konl(1..20)        nodes belonging to the element
!     ipompc(1..nmpc))   ipompc(i) points to the first term of
!                        MPC i in field nodempc
!     nodempc(1,*)       node number of a MPC term
!     nodempc(2,*)       coordinate direction of a MPC term
!     nodempc(3,*)       if not 0: points towards the next term
!                                  of the MPC in field nodempc
!                        if 0: MPC definition is finished
!     coefmpc(*)         coefficient of a MPC term
!     nmpc               number of MPC's
!     ikmpc(1..nmpc)     ordered global degrees of freedom of the MPC's
!                        the global degree of freedom is
!                        8*(node-1)+direction of the dependent term of
!                        the MPC (direction = 0: temperature;
!                        1-3: displacements; 4: static pressure;
!                        5-7: rotations)
!     ilmpc(1..nmpc)     ilmpc(i) is the MPC number corresponding
!                        to the reference number in ikmpc(i)   
!     mi(1)              max # of integration points per element (max
!                        over all elements)
!     mi(2)              max degree of freedomm per node (max over all
!                        nodes) in fields like v(0:mi(2))...
!
!     OUTPUT:
!
!     flux(1)            magnitude of the flux
!     flux(2)            not used; please do NOT assign any value
!     iscale             determines whether the flux has to be
!                        scaled for increments smaller than the 
!                        step time in static calculations
!                        0: no scaling
!                        1: scaling (default)
!           
      implicit none
!
      character*8 lakonl
      character*20 loadtype
!
      integer kstep,kinc,noel,npt,jltyp,konl(20),ipompc(*),
     &  nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),node,idof,id,iscale,mi(*)
!
      real*8 flux(2),time(2),coords(3),sol,temp,press,vold(0:mi(2),*),
     &  area,co(3,*),coefmpc(*)
!
!     the code starting here up to the end of the file serves as
!     an example for combined mechanical-lubrication problems. 
!     Please replace it by your own code for your concrete application.
!
      integer ifaceq(9,6),ifacet(7,4),ifacew(8,5),ig,nelem,nopes,
     &  iflag,i,j,k,nope
!
      real*8 xl21(3,9),xi,et,al,rho,um,h,pnode1(3),pnode2(3),
     &  ratio(9),dist,xl22(3,9),dpnode1(3,3),dpnode2(3,3),v1(3),
     &  v2(3),dh(3),xsj2(3),xs2(3,7),shp2(7,8)
!
      data ifaceq /4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data iflag /3/
!
      include "gauss.f"
!
      nelem=noel
!
      if(lakonl(4:4).eq.'2') then
         nope=20
         nopes=8
      elseif(lakonl(4:4).eq.'8') then
         nope=8
         nopes=4
      elseif(lakonl(4:5).eq.'10') then
         nope=10
         nopes=6
      elseif(lakonl(4:4).eq.'4') then
         nope=4
         nopes=3
      elseif(lakonl(4:5).eq.'15') then
         nope=15
      elseif(lakonl(4:4).eq.'6') then
         nope=6
      endif
!     
!     treatment of wedge faces
!     
      if(lakonl(4:4).eq.'6') then
         if(ig.le.2) then
            nopes=3
         else
            nopes=4
         endif
      endif
      if(lakonl(4:5).eq.'15') then
         if(ig.le.2) then
            nopes=6
         else
            nopes=8
         endif
      endif
!
!     first side of the oil film
!     
      ig=1
!
      if((nope.eq.20).or.(nope.eq.8)) then
         do i=1,nopes
            node=konl(ifaceq(i,ig))
            idof=8*(node-1)+4
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dflux: node ',node
               write(*,*) '       is not connected to the structure'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            do j=1,3
               xl21(j,i)=co(j,node)+
     &              vold(j,node)
            enddo
         enddo
      elseif((nope.eq.10).or.(nope.eq.4)) then
         write(*,*) '*ERROR in dload: tetrahedral elements'
         write(*,*) '       are not allowed'
         call exit(201)
      else
         do i=1,nopes
            node=konl(ifacew(i,ig))
            idof=8*(node-1)+4
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dflux: node ',node
               write(*,*) '       is not connected to the structure'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            do j=1,3
               xl21(j,i)=co(j,node)+
     &              vold(j,node)
            enddo
         enddo
      endif
!
!     projecting the integration point on the first side of the
!     oil film
!
      do j=1,3
         pnode1(j)=coords(j)
      enddo
!
      call attach(xl21,pnode1,nopes,ratio,dist,xi,et)
! 
!     derivative of the shape functions in (xi,et)
!    
      if(nopes.eq.8) then
         call shape8q(xi,et,xl21,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.4) then
         call shape4q(xi,et,xl21,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.6) then
         call shape6tri(xi,et,xl21,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xi,et,xl21,xsj2,xs2,shp2,iflag)
      endif
!
!     the gradient of pnode1 
!     dpnode1(j,k)=dpnode1(j)/dx(k)
!
      do i=1,3
         do j=1,3
            dpnode1(i,j)=0.d0
            do k=1,nopes
               dpnode1(i,j)=dpnode1(i,j)+shp2(j,k)*xl21(i,k)
            enddo
         enddo
      enddo
!
!     second side of the oil film
!     
      ig=2
!
      if((nope.eq.20).or.(nope.eq.8)) then
         do i=1,nopes
            node=konl(ifaceq(i,ig))
            idof=8*(node-1)+4
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dflux: node ',node
               write(*,*) '       is not connected to the structure'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            do j=1,3
               xl22(j,i)=co(j,node)+
     &              vold(j,node)
            enddo
         enddo
      elseif((nope.eq.10).or.(nope.eq.4)) then
         write(*,*) '*ERROR in dload: tetrahedral elements'
         write(*,*) '       are not allowed'
         call exit(201)
      else
         do i=1,nopes
            node=konl(ifacew(i,ig))
            idof=8*(node-1)+4
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in dflux: node ',node
               write(*,*) '       is not connected to the structure'
               call exit(201)
            endif
            node=nodempc(1,nodempc(3,ipompc(ilmpc(id))))
            do j=1,3
               xl22(j,i)=co(j,node)+
     &              vold(j,node)
            enddo
         enddo
      endif
!
!     projecting the integration point on the second side of the
!     oil film
!
      do j=1,3
         pnode2(j)=coords(j)
      enddo
!
      call attach(xl22,pnode2,nopes,ratio,dist,xi,et)
! 
!     derivative of the shape functions in (xi,et)
!    
      if(nopes.eq.8) then
         call shape8q(xi,et,xl22,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.4) then
         call shape4q(xi,et,xl22,xsj2,xs2,shp2,iflag)
      elseif(nopes.eq.6) then
         call shape6tri(xi,et,xl22,xsj2,xs2,shp2,iflag)
      else
         call shape3tri(xi,et,xl22,xsj2,xs2,shp2,iflag)
      endif
!
!     the gradient of pnode1 
!     dpnode2(j,k)=dpnode2(j)/dx(k)
!
      do i=1,3
         do j=1,3
            dpnode2(i,j)=0.d0
            do k=1,nopes
               dpnode2(i,j)=dpnode2(i,j)+shp2(j,k)*xl22(i,k)
            enddo
         enddo
      enddo
!
!     calculating the thickness of the oil film
!
      h=dsqrt((pnode1(1)-pnode2(1))**2+
     &        (pnode1(2)-pnode2(2))**2+
     &        (pnode1(3)-pnode2(3))**2)
!
!     calculating the gradient of the oil film thickness
!
      do i=1,3
         dh(i)=((pnode1(1)-pnode2(1))*(dpnode1(1,i)-dpnode2(1,i))
     &         +(pnode1(2)-pnode2(2))*(dpnode1(2,i)-dpnode2(2,i))
     &         +(pnode1(3)-pnode2(3))*(dpnode1(3,i)-dpnode2(3,i)))/h
      enddo
!
!     velocity of the parts adjoining the film
!     the axis or rotation is assumed to be the x-axis
!
      do i=1,3
         v1(i)=0.d0
      enddo
      v2(1)=0.d0
      v2(2)=-26000.d0*coords(3)/dsqrt(coords(2)**2+coords(3)**2)
      v2(3)=26000.d0*coords(2)/dsqrt(coords(2)**2+coords(3)**2)
!
!     density (oil, N-mm-s-K system)
!
      rho=890.d-9
!
!     body flux
!
      flux(1)=-rho*((v1(1)+v2(1))*dh(1)+
     &              (v1(2)+v2(2))*dh(2)+
     &              (v1(3)+v2(3))*dh(3))/2.d0
!
      iscale=0
!
      return
      end

