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
      subroutine crackrate(nfront,ifrontrel,xkeq,phi,ifront,
     &     dadn,ncyc,icritic,datarget,crcon,temp,ncrtem,
     &     crconloc,ncrconst,xk1,xk2,xk3,nstep,acrack,
     &     wk1,wk2,wk3,xkeqmin,xkeqmax,dkeq,domstep,domphi,
     &     param,nparam,law,ier)
!
!     User Subroutine
!     
!     calculate the crack propagation rate and the number of
!     cycles for this increment
!     
!     INPUT:
!
!     nfront             total number of crack front-nodes
!     ifrontrel(i)       entry of front-node i in field ibounnod
!     xkeq(m,i)          equivalent K-factor for step m in front-node i
!     phi(m,i)           deflection angle for step m in front-node i
!     ifront(i)          (global) node number of front-node i
!     datarget           target crack propagation increment size
!     crcon(0,1..ncrtem) temperature data points for the crack propagation
!                        data
!     crcon(1..ncrconst,j)
!                        crack propagation constants for temperature data
!                        point j
!     temp(j,i)          temperature in step j at boundary node ibounnod(i)
!     ncrtem             number of temperature data points in the material
!                        law
!     crconloc(i)        interpolated crack propagation constant i for a
!                        concrete temperature
!     ncrconst           number of constants in the crack propagation law
!                        for a given temperature
!     xk1(m,i)           K1-factor for step m in front-node i     
!     xk2(m,i)           K2-factor for step m in front-node i     
!     xk3(m,i)           K3-factor for step m in front-node i     
!     nstep              number of steps for which stresses (and      
!                        optionally temperatures) are available
!     acrack(m,i)        crack length for step m in front-node i;  
!                        in some methods the crack length is dependent
!                        on the load, and therefore on the step number;
!     
!     
!     OUTPUT
!
!     dadn(i)            crack propagation rate in front-node i
!     ncyc               number of cycles in this increment
!     icritic            0: Kc is nowhere reached
!                        <0: nowhere propagation      
!                        >0: Kc is reached in at least one front-node
!     wk1(i)             worst K1-factor in front-node i; 
!                        wk1(i) = xk1(m,i) if dabs(xk1(m,i)) is maximal
!                        over all steps
!     wk2(i)             worst K2-factor in front-node i; 
!     wk3(i)             worst K3-factor in front-node i; 
!     xkeqmin(i)         minimum value of xkeq(m,i) over all steps      
!     xkeqmax(i)         maximum value of xkeq(m,i) over all steps      
!     dkeq(i)            range of the equivalent stress intensity factor 
!                        of the main cycle at front-node i
!     domstep(i)         the step dictating the deflection angle at
!                        front-node i (dominant step)
!     domphi(i)          deflection angle for the crack propagation in
!                        front-node i
!
      implicit none
!
      character*132 param(*)
!     
      integer i,nfront,ifrontrel(*),noderel,icritic,ncyc,ncrconst,
     &     ncrtem,nstep,m,ifront(*),nparam,law,ier
!     
      real*8 datarget,damax,xkeq(nstep,*),phi(nstep,*),
     &     acrack(nstep,*),dadn(*),crcon(0:ncrconst,*),t1l,
     &     crconloc(*),dadnref,dkref,xm,epsilon,dkth,delta,dkc,
     &     xk1(nstep,*),xk2(nstep,*),xk3(nstep,*),temp(nstep,*),
     &     wk1(*),wk2(*),wk3(*),
     &     xkeqmin(*),xkeqmax(*),dkeq(*),domstep(*),domphi(*)
!     
!     determine the maximal crack propagation in one iteration
!     
      damax=0.d0
      do i=1,nfront
!
!     look for maximum keq across all steps in the mission; this
!     step is considered to be dominant and the crack propagation
!     in the mission is defined as the crack propagation of this
!     step (same applies to the deflection angle: the deflection
!     angle of the mission is the deflection angle of this step)
!
        dkeq(i)=0.d0
        xkeqmin(i)=1.d30
        xkeqmax(i)=-1.d30
        wk1(i)=0.d0
        wk2(i)=0.d0
        wk3(i)=0.d0
!        
        do m=1,nstep
          if(xkeq(m,i).gt.dkeq(i)) then
            dkeq(i)=xkeq(m,i)
            domphi(i)=phi(m,i)
            domstep(i)=1.d0*m
          endif
          xkeqmin(i)=min(xkeqmin(i),xkeq(m,i))
          xkeqmax(i)=max(xkeqmax(i),xkeq(m,i))
          if(dabs(xk1(m,i)).gt.dabs(wk1(i))) wk1(i)=xk1(m,i)
          if(dabs(xk2(m,i)).gt.dabs(wk2(i))) wk2(i)=xk2(m,i)
          if(dabs(xk3(m,i)).gt.dabs(wk3(i))) wk3(i)=xk3(m,i)
        enddo
!     
!     determine the crack propagation data for the local temperature
!
        noderel=ifrontrel(i)
        t1l=temp(1,noderel)
        call materialdata_crack(crcon,ncrconst,ncrtem,t1l,crconloc)
!
!     constants for a simple material law        
!
        dadnref=crconloc(1)
        dkref=crconloc(2)
        xm=crconloc(3)
        epsilon=crconloc(4)
        dkth=crconloc(5)
        delta=crconloc(6)
        dkc=crconloc(7)
!
        if(dkeq(i).ge.dkc) then
          write(*,*) '*WARNING in crackrate: Kc is reached'
          write(*,*) '         original K-value: ',dkeq(i)
          dkeq(i)=dkth+(dkc-dkth)*0.999
          write(*,*) '         dk is reduced to ',dkeq(i)
          icritic=icritic+1
        endif
!     
        if(dkeq(i).le.dkth) then
          dadn(i)=0.d0
        else
          dadn(i)=dadnref*(dkeq(i)/dkref)**xm*
     &         (1.d0-dexp(epsilon*(1.d0-dkeq(i)/dkth)))/
     &         (1.d0-dexp(delta*(dkeq(i)/dkc-1.d0)))
        endif
!        
        damax=max(dadn(i),damax)
      enddo
!     
!     determine the number of cycles
!
      if(damax.gt.0.d0) then
c        ncyc=nint(max(datarget/damax,1.d0))
        ncyc=nint(datarget/damax)
        if(ncyc.eq.0)then
          ier=1
          return
        endif
      else
!
!       no propagation
!
        icritic=-1
      endif
!
      if(icritic.ne.0) ncyc=1
!     
      return
      end

      
