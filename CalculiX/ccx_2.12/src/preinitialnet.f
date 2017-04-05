!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine preinitialnet(ieg,lakon,v,ipkon,kon,nflow,prop,ielprop,
     &     ielmat,ntmat_,shcon,nshcon,rhcon,nrhcon,mi,iponoel,inoel,
     &     itg,ntg)
!
!     this routine only applies to compressible networks
!
!     determination of initial values based on the boundary conditions
!     and the initial values given by the user by propagating these
!     through the network (using information on the mass flow direction
!     derived from unidirectional network elements or mass flow given
!     by the user (boundary conditions or initial conditions)) 
!
!     a mass flow with size 1.d-30 is used to propagate the sign
!     of the mass flow (in case only the sign and not the size is known)
!
!     it is assumed that pressures, temperatures and mass flows cannot
!     be identically zero (a zero pressure does not make sense,
!     temperature has to be the absolute temperature,
!     a zero mass flow leads to convergence problems).
!     
      implicit none
!
      character*8 lakon(*)
!           
      integer mi(*),ieg(*),nflow,i,ielmat(mi(3),*),ntmat_,node1,node2,
     &     nelem,index,nshcon(*),ipkon(*),kon(*),nodem,imat,ielprop(*),
     &     nrhcon(*),neighbor,ichange,iponoel(*),inoel(2,*),indexe,
     &     itg(*),ntg,node,imin,imax
!     
      real*8 prop(*),shcon(0:3,ntmat_,*),xflow,v(0:mi(2),*),cp,r,
     &     dvi,rho,rhcon(0:1,ntmat_,*),kappa,cti,Ti,ri,ro,p1zp2,omega,
     &     p2zp1,xmin,xmax
!
c      write(*,*) 'preinitialnet '
c      do i=1,ntg
c         write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c      enddo
!
      do
         ichange=0
!
!        propagation of the pressure through the network
!
         do i=1,nflow
            nelem=ieg(i)
            indexe=ipkon(nelem)
!
            node1=kon(indexe+1)
            if(node1.eq.0) cycle
            nodem=kon(indexe+2)
            node2=kon(indexe+3)
            if(node2.eq.0) cycle
!
!           at least one total pressure value unknown in the element
!
            if(((v(2,node1).ne.0.d0).and.(v(2,node2).eq.0.d0)).or.
     &         ((v(2,node1).eq.0.d0).and.(v(2,node2).ne.0.d0))) then
!
               if(lakon(nelem)(2:3).eq.'OR') then
!
!               directional elements: small slope               
!
                  if(v(2,node1).eq.0.d0) then
                     v(2,node1)=v(2,node2)*1.01d0
                  else
                     v(2,node2)=v(2,node1)*0.99d0
                  endif
                  ichange=1
                  if(v(1,nodem).eq.0.d0) v(1,nodem)=1.d-30
               elseif(lakon(nelem)(2:3).eq.'VO') then
!
!                 vortex: pressure ratio can be determined
!                 from geometry
!
                  index=ielprop(nelem)
!
                  if(prop(index+1).lt.prop(index+2)) then
!
!                    r2 < r1
!
                     if(v(0,node2).ne.0.d0) then
                        ri=prop(index+1)
                        ro=prop(index+2)
                        Ti=v(0,node2)
                        if(lakon(nelem)(4:5).eq.'FO') then
                           omega=prop(index+5)
                        else
                           omega=prop(index+7)
                        endif
!
                        imat=ielmat(1,nelem)
                        call materialdata_tg(imat,ntmat_,Ti,shcon,
     &                       nshcon,cp,
     &                       r,dvi,rhcon,nrhcon,rho)
!
                        kappa=cp/(cp-r)
                        cti=omega*ri
!
                        if(lakon(nelem)(4:5).eq.'FO') then
!
!                          forced vortex
!
                           p1zp2=(1.d0+cti**2*((ro/ri)**2-1.d0)/
     &                           (2.d0*cp*Ti))**(kappa/(kappa-1.d0))
                        else
!
!                          free vortex
!
                           p1zp2=(1.d0+cti**2*(1.d0-(ri/ro)**2)/
     &                           (2.d0*cp*Ti))**(kappa/(kappa-1.d0))
                        endif
!
                        if(v(2,node1).eq.0.d0) then
                           v(2,node1)=v(2,node2)*p1zp2
                        else
                           v(2,node2)=v(2,node1)/p1zp2
                        endif
                        ichange=1
                     endif
                  else
!
!                    r1 <= r2
!
                     if(v(0,node1).ne.0.d0) then
                        ri=prop(index+2)
                        ro=prop(index+1)
                        Ti=v(0,node1)
                        if(lakon(nelem)(4:5).eq.'FO') then
                           omega=prop(index+5)
                        else
                           omega=prop(index+7)
                        endif
!
                        imat=ielmat(1,nelem)
                        call materialdata_tg(imat,ntmat_,Ti,shcon,
     &                       nshcon,cp,
     &                       r,dvi,rhcon,nrhcon,rho)
!
                        kappa=cp/(cp-r)
                        cti=omega*ri
!
                        if(lakon(nelem)(4:5).eq.'FO') then
!
!                          forced vortex
!
                           p2zp1=(1.d0+cti**2*((ro/ri)**2-1.d0)/
     &                           (2.d0*cp*Ti))**(kappa/(kappa-1.d0))
                        else
!
!                          free vortex
!
                           p2zp1=(1.d0+cti**2*(1.d0-(ri/ro)**2)/
     &                           (2.d0*cp*Ti))**(kappa/(kappa-1.d0))
                        endif
!
                        if(v(2,node1).eq.0.d0) then
                           v(2,node1)=v(2,node2)/p2zp1
                        else
                           v(2,node2)=v(2,node1)*p2zp1
                        endif
                        ichange=1
                     endif
                  endif
               elseif(v(1,nodem).ne.0.d0) then
!
!                 mass flow is given (either size or just the
!                 direction): small slope
!
                  if(v(1,nodem).gt.0.d0) then
                     if(v(2,node1).eq.0.d0) then
                        v(2,node1)=v(2,node2)*1.01d0
                     else
                        v(2,node2)=v(2,node1)*0.99d0
                     endif
                  else
                     if(v(2,node1).eq.0.d0) then
                        v(2,node1)=v(2,node2)*0.99d0
                     else
                        v(2,node2)=v(2,node1)*1.01d0
                     endif
                  endif
                  ichange=1
               endif
            endif
         enddo
!
!        propagation of the mass flow through the network
!
         do i=1,nflow
            nelem=ieg(i)
            indexe=ipkon(nelem)
            nodem=kon(indexe+2)
!
            if(dabs(v(1,nodem)).le.1.d-30) then
!
!              no initial mass flow given yet
!              check neighbors for mass flow (only if not
!              branch nor joint)
!
!              first end node
!
               node1=kon(indexe+1)
!
               if(node1.ne.0) then
                  index=iponoel(node1)
!
                  if(inoel(2,inoel(2,index)).eq.0) then
!
!                 no branch nor joint; determine neighboring element
!
                     if(inoel(1,index).eq.nelem) then
                        neighbor=inoel(1,inoel(2,index))
                     else
                        neighbor=inoel(1,index)
                     endif
!
!                 initial mass flow in neighboring element
!
                     xflow=v(1,kon(ipkon(neighbor)+2))
!
                     if(dabs(v(1,nodem)).gt.0.d0) then
!
!                    propagate initial mass flow
!
                        if(dabs(xflow).gt.1.d-30) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     else
!
!                    propagate only the sign of the mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     endif
                  endif
               endif
!
!              second end node
!
               node2=kon(indexe+3)
!
               if(node2.ne.0) then
                  index=iponoel(node2)
!
                  if(inoel(2,inoel(2,index)).eq.0) then
!
!                 no branch nor joint; determine neighboring element
!
                     if(inoel(1,index).eq.nelem) then
                        neighbor=inoel(1,inoel(2,index))
                     else
                        neighbor=inoel(1,index)
                     endif
!
!                 initial mass flow in neighboring element
!
                     xflow=v(1,kon(ipkon(neighbor)+2))
!
                     if(dabs(v(1,nodem)).gt.0.d0) then
!
!                    propagate initial mass flow
!
                        if(dabs(xflow).gt.1.d-30) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     else
!
!                    propagate only the sign of the mass flow
!
                        if(dabs(xflow).gt.0.d0) then
                           v(1,nodem)=xflow
                           ichange=1
                           cycle
                        endif
                     endif
                  endif
               endif
            endif
         enddo
!
!        propagation of the temperature
!
         do i=1,nflow
            nelem=ieg(i)
            indexe=ipkon(nelem)
            node1=kon(indexe+1)
            if(node1.eq.0) cycle
            node2=kon(indexe+3)
            if(node2.eq.0) cycle
!
            if(((v(0,node1).ne.0.d0).and.(v(0,node2).ne.0.d0)).or.
     &         ((v(0,node1).eq.0.d0).and.(v(0,node2).eq.0.d0))) cycle
!
!           If the element is a adiabatic gas pipe the
!           total temperature at both ends is equal
!
            if(lakon(nelem)(2:6).eq.'GAPFA') then
               if(v(0,node1).eq.0.d0) then
                  v(0,node1)=v(0,node2)
               else
                  v(0,node2)=v(0,node1)
               endif
               ichange=1
               cycle
            endif
!
            nodem=kon(indexe+2)
!
            if(v(1,nodem).eq.0.d0) then
!
!              direction of mass flow unknown in the element
!
               cycle
            elseif(v(1,nodem).gt.0.d0) then
!
!              positive mass flow (i.e. going from node1 to node2)
!
               if(v(0,node1).eq.0.d0) cycle
!
!              propagating the temperature to node2
!
               v(0,node2)=v(0,node1)
               ichange=1
               cycle
            else
!
!              negative mass flow (i.e. going from node2 to node1)
!
               if(v(0,node2).eq.0.d0) cycle
!
!              propagating the temperature to node1
!
               v(0,node1)=v(0,node2)
               ichange=1
               cycle
            endif
         enddo
c         write(*,*) 'preinitialnet '
c         do i=1,ntg
c            write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c         enddo
         if(ichange.eq.0) exit
      enddo
!
!     set of mass flow of +-1.d-30 to zero
!
      do i=1,nflow
         nelem=ieg(i)
         indexe=ipkon(nelem)
         nodem=kon(indexe+2)
         if(dabs(v(1,nodem)).eq.1.d-30) v(1,nodem)=0.d0
      enddo
!
!     check the pressures: set pressures to zero (i.e. no initial condtion)
!     which lie in between the neighboring pressures => Laplace method
!     is applied in initialnet
!
      loop: do i=1,ntg
         node=itg(i)
!
!        neighboring elements (excluding nodes which do not belong to
!        any network element or just to one network element (middle nodes)
!
         index=iponoel(node)
         if((index.eq.0).or.(inoel(2,index).eq.0)) cycle
!
         imin=node
         imax=node
         xmin=v(2,node)
         xmax=v(2,node)
!
         do
            nelem=inoel(1,index)
            if(lakon(nelem)(2:3).eq.'VO') cycle loop
            indexe=ipkon(nelem)
!
!           neighboring vertex node
!
            if(kon(indexe+1).ne.node) then
               neighbor=kon(indexe+1)
            else
               neighbor=kon(indexe+3)
            endif
            if(neighbor.eq.0) cycle loop
!
!           check its value
!
            if(dabs(v(2,neighbor)).lt.xmin) then
               xmin=dabs(v(2,neighbor))
               imin=neighbor
            elseif(dabs(v(2,neighbor)).gt.xmax) then
               xmax=dabs(v(2,neighbor))
               imax=neighbor
            endif
!
            index=inoel(2,index)
            if(index.eq.0) exit
         enddo
!
!        if value lies in between the neighboring values => assign a
!        negative sign as marker
!
         if((imin.ne.node).and.(imax.ne.node)) then
            v(2,node)=-v(2,node)
         endif
      enddo loop
!
!     set marked values to zero => Laplace equation will be used
!     in initialnet.f
!
      do i=1,ntg
         node=itg(i)
         index=iponoel(node)
         if((index.eq.0).or.(inoel(2,index).eq.0)) cycle
c         if(v(2,node).lt.0.d0) v(2,node)=0.d0
         if(v(2,node).lt.0.d0) v(2,node)=-v(2,node)
      enddo
c         write(*,*) 'preinitialnet end '
c         do i=1,ntg
c            write(*,'(i10,3(1x,e11.4))') itg(i),(v(j,itg(i)),j=0,2)
c         enddo
!
!     same for temperatures?
!     
      return
      end
      
      
