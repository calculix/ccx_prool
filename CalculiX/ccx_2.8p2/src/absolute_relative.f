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
      subroutine absolute_relative(node1,node2,nodem,nelem,lakon,
     &     kon,ipkon, nactdog,identity,ielprop,prop,iflag,v,
     &     xflow,f,nodef,idirf,df,cp,R,physcon,numf,set,mi,iaxial)
!     
!     orifice element
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(4),idirf(4),index,iflag,
     &     ipkon(*),kon(*),nelemswirl,mi(*),iaxial
!     
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(4),kappa,R,
     &     p1,p2,T1,T2,cp,physcon(*),km1,kp1,kdkm1,
     &     kdkp1,u,pi,Qred_crit,pt1,pt2,Tt1,Tt2,ct,fact,
     &     Cp_cor
!     
      if (iflag.eq.0) then
         identity=.true.
!     
         if(nactdog(2,node1).ne.0)then
            identity=.false.
         elseif(nactdog(2,node2).ne.0)then
            identity=.false.
         elseif(nactdog(1,nodem).ne.0)then
            identity=.false.
         endif
!     
      elseif (iflag.eq.1)then
!     
         kappa=(cp/(cp-R))
         pi=4.d0*datan(1.d0)
         index=ielprop(nelem)
         qred_crit=dsqrt(kappa/R)*
     &        (1+0.5d0*(kappa-1))**(-0.5*(kappa+1)/(kappa-1))
!     
!     Because the flow value is independant of the chosen 
!     coordinate system  initial mass flow value is set to
!     dsqrt(T1)*P1*Qred_crit with Qred_crit/2 = 0.02021518917
!     with consideration to flow direction
!     
         node1=kon(ipkon(nelem)+1)
         node2=kon(ipkon(nelem)+3)
         p1=v(2,node1)
         p2=v(2,node2)
         T1=v(0,node1)
         T2=v(0,node2)
!
        if(p1.gt.p2) then
            xflow=0.75/dsqrt(T1)*P1*qred_crit
         else
            xflow=-0.75/dsqrt(T1)*P1*qred_crit
         endif 
!     
      elseif(iflag.eq.2) then
!     
         numf=4
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!     
         index=ielprop(nelem)
!         
         u=prop(index+1)
         ct=prop(index+2)
!
        if(ct.eq.0) then
            nelemswirl=prop(index+3)
!
!     previous element is a preswirl nozzle
!
            if(lakon(nelemswirl)(2:5).eq.'ORPN') then
               ct=prop(ielprop(nelemswirl)+5)
!
!     previous element is a forced vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
               ct=prop(ielprop(nelemswirl)+7)
!
!     previous element is a free vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
               ct=prop(ielprop(nelemswirl)+9)
            endif
         endif                
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(lakon(nelem)(2:4).eq.'ATR') then
!     
            if(u/CT.ge.2) then
!     
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!     
               nodef(1)=node1
               nodef(2)=node1
               nodef(3)=nodem
               nodef(4)=node2
!     
!     in the case of a negative flow direction
!     
               if(xflow.le.0d0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=ABSOLUTE TO RELATIVE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!     
            else
               pt1=v(2,node2)
               pt2=v(2,node1)
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!
               if(xflow.le.0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=ABSOLUTE TO RELATIVE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!
               nodef(1)=node2
               nodef(2)=node2
               nodef(3)=nodem
               nodef(4)=node1
            endif
!       
         elseif(lakon(nelem)(2:4).eq.'RTA') then
!
            if(u/CT.lt.2) then
!     
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!               
               nodef(1)=node1
               nodef(2)=node1
               nodef(3)=nodem
               nodef(4)=node2
!
               if(xflow.le.0d0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=RELATIVE TO ABSOLUTE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!
            else
!
               pt1=v(2,node2)
               pt2=v(2,node1)
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!
               if(xflow.le.0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=RELATIVE TO ABSOLUTE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!
               nodef(1)=node2
               nodef(2)=node2
               nodef(3)=nodem
               nodef(4)=node1
!              
            endif
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     computing temperature corrected Cp=Cp(T) coefficient 
         call cp_corrected(cp,Tt1,Tt2,cp_cor)
!     
            if(Tt1.lt.273) then
               Tt1= Tt2
            endif
!
            if(cp_cor.eq.0) then
               cp_cor=cp
            endif
!     
!     transformation from absolute system to relative system
!     
         if(lakon(nelem)(2:4).eq.'ATR') then
!     
            fact=1+(u**2-2*u*ct)/(2*Cp_cor*Tt1)
!     
            f=Pt2-Pt1*(fact)**kdkm1
!     
!     pressure node 1
!
            df(1)=-fact**kdkm1
!
!     temperature node1
!
            df(2)=-pt1*Kdkm1*(-(u**2-2*u*ct)/(2*Cp_cor*Tt1**2))
     &           *fact**(kdkm1-1)
!
!     mass flow node m
!
            df(3)=0
!
!     pressure node 2
!     
            df(4)=1
!     
!     transformation from relative system to absolute system
!     
         elseif(lakon(nelem)(2:4).eq.'RTA') then
!     
            fact=1-(u**2-2*u*ct)/(2*Cp*Tt1)
!     
            f=Pt2-Pt1*(fact)**kdkm1
!     
            df(1)=-fact**kdkm1
!     
            df(2)=-Pt1*Kdkm1*((u**2-2*u*ct)/(2*Cp*Tt1**2))
     &           *fact**(kdkm1-1)
!     
            df(3)=0
!     
            df(4)=1
!     
         endif

      elseif(iflag.eq.3) then
            
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!     
         index=ielprop(nelem)
!         
         u=prop(index+1)
         ct=prop(index+2)
!
        if(ct.eq.0) then
            nelemswirl=prop(index+3)
!
!     previous element is a preswirl nozzle
!
            if(lakon(nelemswirl)(2:5).eq.'ORPN') then
               ct=prop(ielprop(nelemswirl)+5)
!
!     previous element is a forced vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFO') then
               ct=prop(ielprop(nelemswirl)+7)
!
!     previous element is a free vortex
!
            elseif(lakon(nelemswirl)(2:5).eq.'VOFR') then
               ct=prop(ielprop(nelemswirl)+9)
            endif
         endif                
!     
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(lakon(nelem)(2:4).eq.'ATR') then
!     
            if(u/CT.ge.2) then
!     
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!     
               nodef(1)=node1
               nodef(2)=node1
               nodef(3)=nodem
               nodef(4)=node2
!     
!     in the case of a negative flow direction
!     
               if(xflow.le.0d0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=ABSOLUTE TO RELATIVE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!     
            else
               pt1=v(2,node2)
               pt2=v(2,node1)
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!
               if(xflow.le.0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=ABSOLUTE TO RELATIVE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!
               nodef(1)=node2
               nodef(2)=node2
               nodef(3)=nodem
               nodef(4)=node1
            endif
!       
         elseif(lakon(nelem)(2:4).eq.'RTA') then
!
            if(u/CT.lt.2) then
!     
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!               
               nodef(1)=node1
               nodef(2)=node1
               nodef(3)=nodem
               nodef(4)=node2
!
               if(xflow.le.0d0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=RELATIVE TO ABSOLUTE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!
            else
!
               pt1=v(2,node2)
               pt2=v(2,node1)
               xflow=v(1,nodem)*iaxial
               Tt1=v(0,node1)+physcon(1)
               Tt2=v(0,node2)+physcon(1)
!
               if(xflow.le.0) then
                  write(*,*)''
                  write(*,*)'*WARNING:'
                  write(*,*)'in element',nelem
                  write(*,*)'TYPE=RELATIVE TO ABSOLUTE'
                  write(*,*)'mass flow negative!'
                  write(*,*)'check results and element definition'
               endif
!
               nodef(1)=node2
               nodef(2)=node2
               nodef(3)=nodem
               nodef(4)=node1
!              
            endif
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
!     
!     computing temperature corrected Cp=Cp(T) coefficient 
         call cp_corrected(cp,Tt1,Tt2,cp_cor)
!     
         if(Tt1.lt.273) then
            Tt1= Tt2
         endif
!     
         if(cp_cor.eq.0) then
            cp_cor=cp
         endif

                  write(1,*) ''
         write(1,55) 'In line',int(nodem/100),' from node',node1,
     &' to node', node2,':   air massflow rate=',xflow,'kg/s'
!     &,', oil massflow rate=',xflow_oil,'kg/s'
 55      FORMAT(1X,A,I6.3,A,I6.3,A,I6.3,A,F9.6,A,A,F9.6,A)

!         if(inv.eq.1) then
            write(1,56)'       Inlet node ',node1,':     Tt1= ',Tt1,
     &           'K, Ts1= ',Tt1,'K, Pt1= ',Pt1/1E5,
     &           'Bar'
            write(1,*)'             element T    ',set(numf)(1:20)
            write(1,57)'             u= ',u,'m/s ,Ct= ',Ct,'m/s'
            write(1,56)'       Outlet node ',node2,':    Tt2= ',T2,
     &           'K, Ts2= ',Tt2,'K, Ptt2= ',Pt2/1e5,
     &           'Bar'
!     
 56      FORMAT(1X,A,I6.3,A,f6.1,A,f6.1,A,f9.5,A,f9.5)  
 57      FORMAT(1X,A,f6.2,A,f6.2,A)

      endif
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return 
      end
