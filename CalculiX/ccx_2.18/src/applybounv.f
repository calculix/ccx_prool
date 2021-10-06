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
      subroutine applybounv(nodeboun,ndirboun,nboun,v,
     &     compressible,nmpc,nodempc,ipompc,coefmpc,inomat,
     &     mi,coefmodmpc)
!     
!     1) applies velocity SPC's for rho*v (rho is assumed constant in
!        step 1 of the algorithm)
!     2) applies MPC's for rho*v [rho is assumed constant in step 1 of
!     the algorithm and it has the same value for all terms in the
!     MPC since the only MPC's allowed are:
!     a) SPC's in nonglobal coordinates (density at only one location)
!     b) cyclic symmetry MPC's (density on both sides equal)]
!     
      implicit none
!     
      integer compressible,mi(*),nodeboun(*),ndirboun(*),i,nboun,node,
     &     index,nodei,nmpc,nodempc(3,*),ipompc(*),ist,ndir,ndiri,
     &     inomat(*)
!     
      real*8 xnorm,coefmodmpc(*),sum,v(0:mi(2),*),coefmpc(*),residu,
     &     correction
!     
!     SPC's: velocity
!     
      do i=1,nboun
!     
!     only velocity boundary conditions
!     
        ndir=ndirboun(i)
        if((ndir.lt.1).or.(ndir.gt.3)) cycle
!     
!     check whether fluid node
!     
        node=nodeboun(i)
        if(inomat(node).eq.0) cycle
!     
!     calculating the physical variables for the node at stake
!     
        v(ndir,node)=0.d0
      enddo
!     
!     MPC's: velocity
!     
      if(compressible.eq.0) then
!     
!     incompressible fluids
!     
        do i=1,nmpc
          ist=ipompc(i)
!     
          ndir=nodempc(2,ist)
          if((ndir.lt.1).or.(ndir.gt.3)) cycle
!     
!     check whether fluid node
!     
          node=nodempc(1,ist)
          if(inomat(node).eq.0) cycle
!     
!     calculating the value of the dependent DOF of the MPC
!     
          index=nodempc(3,ist)
          if(index.eq.0) cycle
          sum=0.d0
          do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
            sum=sum+coefmpc(index)*v(ndiri,nodei)
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
!     
!     calculating the physical variables for the node at stake
!     
          v(ndir,node)=-sum/coefmpc(ist)
        enddo
      else
!     
!     compressible fluids
!     
!     MPC's are treated by distributing the residual proportional to
!     the coefficients
!     
        do i=1,nmpc
          index=ipompc(i)
!     
!     only velocity dofs
!     
          ndir=nodempc(2,index)
          if((ndir.lt.1).or.(ndir.gt.3)) cycle
!     
!     check whether fluid node
!     
          node=nodempc(1,index)
          if(inomat(node).eq.0) cycle
!     
!     calculating the value of the dependent DOF of the MPC
!     
          residu=coefmpc(index)*v(ndir,node)
          xnorm=1.d0
          if(index.eq.0) cycle
          do
            index=nodempc(3,index)
            if(index.eq.0) exit
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
            residu=residu+coefmpc(index)*v(ndiri,nodei)
          enddo
!     
!     correcting all terms of the MPC
!     
          index=ipompc(i)
          do
            nodei=nodempc(1,index)
            ndiri=nodempc(2,index)
!     
            correction=-residu*coefmodmpc(index)
            v(ndiri,nodei)=v(ndiri,nodei)+correction
            index=nodempc(3,index)
            if(index.eq.0) exit
          enddo
        enddo
!     
      endif
!     
      return
      end
      
