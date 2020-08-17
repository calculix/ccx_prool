!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2020 Guido Dhondt
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
      subroutine phys2con(inomat,node,vold,ntmat_,shcon,nshcon,physcon,
     &  compressible,vcon,rhcon,nrhcon,ithermal,mi)
!
!     applies temperature and velocity boundary conditions for
!     incompressible fluids (liquids)
!
      implicit none
!
      integer inomat(*),node,ntmat_,nshcon(*),compressible,nrhcon(*),
     &  ithermal(*),mi(*),imat,k
!
      real*8 rhcon(0:1,ntmat_,*),shcon(0:3,ntmat_,*),vold(0:mi(2),*),
     &  vcon(0:4,*),physcon(*),temp,cp,r,rho
!
      imat=inomat(node)
      temp=vold(0,node)
      call materialdata_cp_sec(imat,ntmat_,temp,shcon,nshcon,
     &     cp,physcon)
!     
      if(compressible.eq.1) then
!     
!        gas
!     
         r=shcon(3,1,imat)
         rho=vold(4,node)/(r*(vold(0,node)-physcon(1)))
         vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &        (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &        /2.d0)-vold(4,node)
         vcon(4,node)=rho
      else
!     
!        liquid 
!     
         call materialdata_rho(rhcon,nrhcon,imat,rho,
     &        temp,ntmat_,ithermal)
         if(ithermal(1).gt.1)
     &        vcon(0,node)=rho*(cp*(temp-physcon(1))+
     &        (vold(1,node)**2+vold(2,node)**2+vold(3,node)**2)
     &        /2.d0)
      endif
!     
      do k=1,3
         vcon(k,node)=rho*vold(k,node)
      enddo
!     
      return
      end
      
