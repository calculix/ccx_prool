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
      subroutine flux(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,rho,physcon,g,co,dvi,numf,
     &     vold,set,shcon,nshcon,rhcon,nrhcon,ntmat_,mi,ider,iaxial)
!
!     determine whether the flux in the element is an unknown 
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!      
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(8),idirf(8),kflag,ipkon(*),kon(*),
     &     nshcon(*), nrhcon(*),ntmat_,mi(*),ider,iaxial
!      
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(8),R,cp,physcon(*),rho,
     &     g(3),co(3,*),dvi,vold(0:mi(2),*),shcon(0:3,ntmat_,*),
     &     rhcon(0:1,ntmat_,*)
!
      if(lakon(nelem)(2:8).eq.'ACCTUBO') then 
!         
         call acctube_one(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,mi,ider,
     &        iaxial)
!
      elseif(lakon(nelem)(2:8).eq.'ACCTUBE') then 
!
         call acctube(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,mi,ider,
     &        iaxial)
!
      elseif(lakon(nelem)(2:6).eq.'CARBS') then  
!
         call carbon_seal(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,R,physcon,dvi,numf,set,mi,iaxial)
!     
      elseif(lakon(nelem)(2:5).eq.'CHAR') then 
!     
         call characteristic(node1,node2,nodem,nelem,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,physcon,numf,set,mi,iaxial)
!
      elseif(lakon(nelem)(2:5).eq.'CROS') then 
!         
         call cross_split(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,iaxial)
!     
      elseif(lakon(nelem)(2:8).eq.'REBRSI1') then 
!         
         call wye(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,iaxial)
!     
      elseif(lakon(nelem)(2:8).eq.'REBRSI2') then 
!         
         call tee(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,r,physcon,numf,set,mi,ider,iaxial)
!
      elseif(lakon(nelem)(2:4).eq.'LAB') then 
!         
         call labyrinth(node1,node2,nodem,nelem,lakon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,physcon,co,dvi,numf,vold,set,
     &     kon,ipkon,mi,iaxial)
! 
      elseif(lakon(nelem)(2:5).eq.'GAPI') then 
!         
         call gaspipe(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,mi,iaxial)
!
      elseif(lakon(nelem)(2:5).eq.'GAPF') then 

         call gaspipe_fanno(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,iaxial)
!
      elseif(lakon(nelem)(2:5).eq.'FDPF') then 
         call free_disc_pumping(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,iaxial)
!   
      elseif(lakon(nelem)(2:5).eq.'FCVF') then 
         call free_convection(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,iaxial)
!
      elseif(lakon(nelem)(2:5).eq.'MFPC') then 
         call massflow_percent(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,iaxial)
         
      elseif(lakon(nelem)(2:5).eq.'LIPI') then
!         
         call liquidpipe(node1,node2,nodem,nelem,lakon,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,dvi,numf,vold,mi,ipkon,kon,set,iaxial)
!
      elseif(lakon(nelem)(2:5).eq.'LICH') then
!         
         call liquidchannel(node1,node2,nodem,nelem,lakon,nactdog,
     &           identity,ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,dvi,numf,mi,ipkon,kon)
!
      elseif(lakon(nelem)(2:5).eq.'LIPU') then
!         
         call liquidpump(node1,node2,nodem,nelem,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,numf,mi,iaxial)
!
      elseif(lakon(nelem)(2:3).eq.'LP') then
!         
         call liquidpipe(node1,node2,nodem,nelem,lakon,nactdog,identity,
     &           ielprop,prop,kflag,v,xflow,f,nodef,idirf,df,
     &           rho,g,co,dvi,numf,vold,mi,ipkon,kon,set,iaxial)
!      
      elseif(lakon(nelem)(2:4).eq.'MRG') then 
!     
         call moehring(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,dvi,numf,set,mi,iaxial)
!
      elseif(lakon(nelem)(2:3).eq.'OR') then 
!         
         call orifice(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,co,vold,mi,
     &        iaxial)
!
      elseif(lakon(nelem)(2:3).eq.'RE') then 
!         
         call restrictor(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,mi,iaxial)
! 
      elseif(lakon(nelem)(2:4).eq.'RCV') then   
!
         call rcavi(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,dvi,numf,set,mi,iaxial)
!
      elseif(lakon(nelem)(2:3).eq.'RO') then   
!
         call rcavi2(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &     nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &     nodef,idirf,df,cp,R,dvi,numf,set,mi,iaxial)
!
      elseif(lakon(nelem)(2:5).eq.'RIMS') then   
!
         call rimseal(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,mi,iaxial)
!     
      elseif(lakon(nelem)(2:6).eq.'SPUMP') then 
! 
        call scavenge_pump(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,ntmat_,mi,
     &        iaxial)   
!
      elseif(lakon(nelem)(2:3).eq.'VO') then 
!     
         call vortex(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,numf,set,mi,iaxial) 
!     
      elseif((lakon(nelem)(2:4).eq.'ATR')
     &        .or.(lakon(nelem)(2:4).eq.'RTA')) then
!     
         call absolute_relative(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,kflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,numf,set,mi,iaxial) 
!
      else
         identity=.true.
!          
      endif
!    
      return
      end
      
