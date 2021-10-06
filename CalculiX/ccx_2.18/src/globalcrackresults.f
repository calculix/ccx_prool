      
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
      subroutine globalcrackresults(nfront,ifront,wk1,wk2,wk3,dkeq,
     &     domphi,dadn,ncyc,wk1glob,wk2glob,wk3glob,dkeqglob,phiglob,
     &     dadnglob,dnglob,acrack,acrackglob,nstep,xkeqmin,xkeqmax,
     &     xkeqminglob,xkeqmaxglob,iinc,iincglob,domstep,domstepglob)
!     
!     calculate the stress intensity factors along the crack fronts
!     
      implicit none
!     
      integer nfront,ifront(*),i,node,ncyc,nstep,iinc,iincglob(*)
!     
      real*8 wk1(*),wk2(*),wk3(*),dkeq(*),domphi(*),dadn(*),
     &     wk1glob(*),wk2glob(*),wk3glob(*),dkeqglob(*),phiglob(*),
     &     dadnglob(*),dnglob(*),acrackglob(*),acrack(nstep,*),
     &     xkeqmax(*),xkeqminglob(*),xkeqmaxglob(*),domstep(*),
     &     xkeqmin(*),domstepglob(*)
!
      do i=1,nfront
!
!     loop over all nodes belonging to the crack front(s)
!     the crack length is the one corresponding to step 1;
!     although acrackglob(node) already has values from the previous
!     increment (except in the first increment), the value of
!     acrack(1,i) may be different due to cracklength_smoothing;
!
        node=ifront(i)
        wk1glob(node)=wk1(i)
        wk2glob(node)=wk2(i)
        wk3glob(node)=wk3(i)
        xkeqminglob(node)=xkeqmin(i)
        xkeqmaxglob(node)=xkeqmax(i)
        dkeqglob(node)=dkeq(i)
        phiglob(node)=domphi(i)
        dadnglob(node)=dadn(i)
c        dnglob(node)=1.d0*ncyc
        acrackglob(node)=acrack(1,i)
        iincglob(node)=iinc
        domstepglob(node)=domstep(i)
      enddo
!      
      return
      end

