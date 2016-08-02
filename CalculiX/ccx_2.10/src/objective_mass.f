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
      subroutine objective_mass(co,kon,ipkon,lakon,v,nelcon,rhcon,
     &  ielmat,ielorien,norien,ntmat_,vold,matname,mi,
     &  nal,thicke,mortar,nea,neb,ielprop,prop,distmin,
     &  ndesi,nodedesi,ndirdesi,nobject,g0,dgdxtot,numobject)
!
!     calculates the total differential of the objective function MASS
!
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*)
!
      integer kon(*),nea,neb,mi(*),nelcon(2,*),ielmat(mi(3),*),
     &  ielorien(mi(3),*),ntmat_,ipkon(*),norien,nal,ielprop(*),
     &  mortar,ndesi,nodedesi(*),ndirdesi(*),nobject,numobject
!
      real*8 co(3,*),v(0:mi(2),*),prop(*),rhcon(0:1,ntmat_,*),
     &  vold(0:mi(2),*),thicke(mi(3),*),distmin,g0(nobject),
     &  dgdxtot(ndesi,nobject)
!
!     -------------------------------------------------------------
!     Calculation of the derivative w.r.t to the designvaribales
!     (For the objective function mass this is equal to the
!     total difference
!     -------------------------------------------------------------
!
      call objective_mass_dx(co,kon,ipkon,lakon,v,nelcon,rhcon,
     &  ielmat,ielorien,norien,ntmat_,vold,matname,mi,nal,thicke,
     &  mortar,nea,neb,ielprop,prop,distmin,ndesi,nodedesi,ndirdesi,
     &  nobject,g0,dgdxtot,numobject)
!
      return
      end
