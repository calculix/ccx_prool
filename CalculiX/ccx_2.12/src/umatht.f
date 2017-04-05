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
      subroutine umatht(u,dudt,dudg,flux,dfdt,dfdg,
     &  statev,temp,dtemp,dtemdx,time,dtime,predef,dpred,
     &  cmname,ntgrd,nstatv,props,nprops,coords,pnewdt,
     &  noel,npt,layer,kspt,kstep,kinc,vold,co,lakonl,konl,
     &  ipompc,nodempc,coefmpc,nmpc,ikmpc,ilmpc,mi)
!
!     heat transfer material subroutine
!
!     INPUT:
!
!     statev(nstatv)      internal state variables at the start
!                         of the increment
!     temp                temperature at the start of the increment
!     dtemp               increment of temperature
!     dtemdx(ntgrd)       current values of the spatial gradients of the 
!                         temperature
!     time(1)             step time at the beginning of the increment
!     time(2)             total time at the beginning of the increment
!     dtime               time increment
!     predef              not used
!     dpred               not used
!     cmname              material name
!     ntgrd               number of spatial gradients of temperature
!     nstatv              number of internal state variables as defined
!                         on the *DEPVAR card
!     props(nprops)       user defined constants defined by the keyword
!                         card *USER MATERIAL,TYPE=THERMAL
!     nprops              number of user defined constants, as specified
!                         on the *USER MATERIAL,TYPE=THERMAL card
!     coords              global coordinates of the integration point
!     pnewd               not used
!     noel                element number
!     npt                 integration point number
!     layer               not used
!     kspt                not used
!     kstep               not used
!     kinc                not used
!     vold(0..4,1..nk)    solution field in all nodes
!                         0: temperature
!                         1: displacement in global x-direction
!                         2: displacement in global y-direction
!                         3: displacement in global z-direction
!                         4: static pressure
!     co(3,1..nk)         coordinates of all nodes
!                         1: coordinate in global x-direction
!                         2: coordinate in global y-direction
!                         3: coordinate in global z-direction
!     lakonl              element label
!     konl(1..20)         nodes belonging to the element
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
!     u                   not used
!     dudt                not used
!     dudg(ntgrd)         not used
!     flux(ntgrd)         heat flux at the end of the increment
!     dfdt(ntgrd)         not used
!     dfdg(ntgrd,ntgrd)   variation of the heat flux with respect to the 
!                         spatial temperature gradient
!     statev(nstatv)      internal state variables at the end of the
!                         increment
!
      implicit none
!
      character*8 lakonl
      character*80 cmname
!
      integer ntgrd,nstatv,nprops,noel,npt,layer,kspt,kstep,kinc,
     &  konl(20),ipompc(*),nodempc(3,*),nmpc,ikmpc(*),ilmpc(*),mi(*)
!
      real*8 u,dudt,dudg(ntgrd),flux(ntgrd),dfdt(ntgrd),
     &  statev(nstatv),pnewdt,temp,dtemp,dtemdx(ntgrd),time(2),dtime,
     &  predef,dpred,props(nprops),coords(3),dfdg(ntgrd,ntgrd),
     &  vold(0:mi(2),*),co(3,*),coefmpc(*)
!
!     the code starting here up to the end of the file serves as
!     an example for combined mechanical-lubrication problems. 
!     Please replace it by your own code for your concrete application.
!
      integer ifaceq(8,6),ifacet(6,4),ifacew(8,5),ig,nelem,nopes,
     &  iflag,i,j,nope,node,idof,id
!
      real*8 xl21(3,9),xi,et,al,rho,um,h,pnode1(3),pnode2(3),
     &  ratio(9),dist,xl22(3,9)
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
      data iflag /2/
!
      nelem=noel
      i=npt
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
               write(*,*) '*ERROR in umatht: node ',node
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
         write(*,*) '*ERROR in umatht: tetrahedral elements'
         write(*,*) '       are not allowed'
         call exit(201)
      else
         do i=1,nopes
            node=konl(ifacew(i,ig))
            idof=8*(node-1)+4
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in umatht: node ',node
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
               write(*,*) '*ERROR in umatht: node ',node
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
         write(*,*) '*ERROR in umatht: tetrahedral elements'
         write(*,*) '       are not allowed'
         call exit(201)
      else
         do i=1,nopes
            node=konl(ifacew(i,ig))
            idof=8*(node-1)+4
            call nident(ikmpc,idof,nmpc,id)
            if((id.eq.0).or.(ikmpc(id).ne.idof)) then
               write(*,*) '*ERROR in umatht: node ',node
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
!     calculating the thickness of the oil film
!
      h=dsqrt((pnode1(1)-pnode2(1))**2+
     &        (pnode1(2)-pnode2(2))**2+
     &        (pnode1(3)-pnode2(3))**2)
!
!     density, viscosity (oil, SI units, 290 K)
!
      rho=890.d-9
      um=1.d-6
!
      al=(h**3)*rho/(12.d0*um)
!
!     filling the tangent matrix
!
      do i=1,3
         do j=1,3
            dfdg(i,j)=0.d0
         enddo
         dfdg(i,i)=al
      enddo
!
!     determining the equivalent flux
!
      do j=1,ntgrd
         flux(j)=-al*dtemdx(j)
      enddo
!
      return
      end
