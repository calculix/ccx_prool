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
      subroutine gaspipe(node1,node2,nodem,nelem,lakon,kon,ipkon,
     &        nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,shcon,
     &        nshcon,rhcon,nrhcon,ntmat_,mi,iaxial)
!     
!     pipe with friction losses 
!
!     author: Yannick Muller
!     
      implicit none
!     
      logical identity,crit
      character*8 lakon(*)
      character*81 set(*)
!     
      integer nelem,nactdog(0:3,*),node1,node2,nodem,numf,
     &     ielprop(*),nodef(5),idirf(5),index,iflag,
     &     inv,ipkon(*),kon(*),icase,kgas,k_oil,nshcon(*),
     &     nrhcon(*),ntmat_,mi(*),iaxial
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(5),kappa,R,a,d,l,
     &     p1,p2,T1,T2,Tt1,Tt2,pt1,pt2,cp,physcon(*),p2p1,km1,dvi,
     &     kp1,kdkm1,reynolds,pi,e,lambda,lld,kdkp1,T2dTt2,
     &     T1dTt1,X_t1dTt1,X_t2dTt2,X2_den,X1_den,
     &     X1,X2,B1,B2,C1,C2,t_moy,tdkp1,ln,m2r2d2a2,
     &     pt2zpt1,ks,form_fact,Tt1dT1,Tt2dT2,M1,M2,
     &     Pt2zPt1_c,qred_crit,l_neg,Qred,Ts1,qred_max1,phi,xflow_oil,
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*)
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
         crit=.false.
!     
         index=ielprop(nelem)
         kappa=(cp/(cp-R))
         A=prop(index+1)
         d=prop(index+2)
         l=prop(index+3)
         if(l.lt.0d0) then
            l_neg=l
            l=abs(l)
         else
            l_neg=l
         endif
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPIA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPII') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=int(prop(index+7))
!
         p1=v(2,node1)
         p2=v(2,node2)
!
         if(p1.ge.p2) then
            inv=1
            T1=v(0,node1)+physcon(1)
            T2=v(0,node2)+physcon(1)
         else
            inv=-1
            p1=v(2,node2)
            p2=v(2,node1)
            T1=v(0,node2)+physcon(1)
            T2=v(0,node1)+physcon(1)
         endif
!
         p2p1=p2/p1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
!             
         if(p2p1.gt.C2) then
            xflow=inv*p1*a*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &           *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(T1)
         else
            xflow=inv*p1*a*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(T1)
         endif
!
!     calculation of the dynamic viscosity 
!
         if(dabs(dvi).lt.1E-30) then
            kgas=0
            call dynamic_viscosity(kgas,T1,dvi)
         endif  
!
         reynolds=dabs(xflow)*d/(dvi*a)
!     
         if(reynolds.lt.100) then
            reynolds = 100
         endif
!
         call friction_coefficient(l_neg,d,ks,reynolds,form_fact,lambda)
!
         call pt2zpt1_crit(p2,p1,T1,T2,lambda,kappa,r,l,d,A,iflag,
     &        inv,Pt2zPt1_c,Qred_crit,crit,qred_max1,icase)
!     
!        next location is "misused" to store the critical value
!        used in resultgas.f
!
         v(2,nodem)=Pt2zPt1_c
!
         Qred=dabs(xflow)*dsqrt(T1)/(A*P1)
         if (crit) then
            xflow=0.5*inv*Qred_crit*P1*A/dsqrt(T1)
!
            if(lakon(nelem)(2:6).eq.'GAPII') then
               call ts_calc(xflow,T1,P1,kappa,r,a,Ts1,icase)
               if (inv.eq.1) then 
                  v(3,node1)=Ts1
                  v(3,node2)=Ts1
                  v(0,node2)=Ts1*(1.d0+km1/(2*kappa))
               else
                  v(3,node2)=Ts1
                  v(3,node1)=Ts1
                  v(0,node1)=Ts1*(1.d0+km1/(2*kappa))  
               endif
            endif
         elseif(Qred.gt.Qred_crit) then
            xflow=0.5*inv*Qred_crit*P1*A/dsqrt(T1)
         else
            xflow=inv*Qred*P1*A/dsqrt(T1)
         endif
!
      elseif (iflag.eq.2)then
!
         numf=5
         crit=.false.
!
         pi=4.d0*datan(1.d0)
         e=2.7182818d0
!
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!
         index=ielprop(nelem)
         A=prop(index+1)
         d=prop(index+2)
         l=prop(index+3)
         if(l.lt.0d0) then
            l_neg=l
            l=abs(l)
         else
            l_neg=l
         endif
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPIA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPII') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=int(prop(index+7)) 
!
         pt1=v(2,node1)
         pt2=v(2,node2)
         xflow=v(1,nodem)*iaxial 
!
         if(xflow.ge.0d0) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)
!
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
!
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
!     
            nodef(1)=node1
            nodef(2)=node1
            nodef(3)=nodem
            nodef(4)=node2
            nodef(5)=node2
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            Tt1=v(0,node2)+physcon(1)
            Tt2=v(0,node1)+physcon(1)
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
!
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
!     
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
            nodef(5)=node1
         endif
!     
         idirf(1)=2
         idirf(2)=0
         idirf(3)=1
         idirf(4)=2
         idirf(5)=0
!
         pt2zpt1=pt2/pt1
!     
!     calculation of the dynamic viscosity 
!     
         if(xflow_oil.ne.0d0) then
!     
            if((k_oil.lt.0).or.(k_oil.gt.12)) then
               write(*,*) '*ERROR:in gaspipe.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the type of oil is not defined'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Current calculation stops here'
               call exit(201)
            else
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
               lambda=lambda*phi
            endif
!     
!     for pure air
!     
         else
            if(dabs(dvi).lt.1E-30) then
               kgas=0
               call dynamic_viscosity(kgas,T1,dvi)
            endif
!     
            reynolds=dabs(xflow)*d/(dvi*a)
!
            phi=1.d0
            call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &           lambda)
         endif
!     
         call pt2zpt1_crit(pt2,pt1,Tt1,Tt2,lambda,kappa,r,l,d,A,iflag,
     &     inv,pt2zpt1_c,qred_crit,crit,qred_max1,icase)
!
         if(dabs(xflow)*dsqrt(Tt1)/(A*Pt1).gt.qred_max1) then
            crit=.true.
         endif
!     
!     definition of the coefficients 
!     
         lld=lambda*l/d 
!
         if(.not.crit) then
!     
            T_moy=0.5d0*(T1+T2)
            T2dTt2=T2/Tt2
            Tt2dT2=1.d0/T2dTt2
            X_T2dTt2=T2dTt2**(2*kdkm1)
            T1dTt1=T1/Tt1
            Tt1dT1=1.d0/T1dTt1
            X_T1dTt1=T1dTt1**(2*kdkm1)
!     
            X2_den=pt2**2*X_T2dTt2
            X2=t2**2/X2_den
            X1_den=pt1**2*X_T1dTt1
            X1=T1**2/X1_den        
!     
            ln=log(Pt2zPt1*(T2dTt2/T1dTt1)**kdkm1)
!     
            m2r2d2a2=xflow**2*R**2/(2*A**2)
!     
            C1=2.d0*cp*A**2*X1_den*(1.d0-2.d0*kdkm1*(Tt1dT1-1.d0))
     &           +2.d0*xflow**2*R**2*T1
!
            C2=2.d0*cp*A**2*X2_den*(1.d0-2.d0*kdkm1*(Tt2dT2-1.d0))
     &           +2.d0*xflow**2*R**2*T2
!     
            B1=-2.d0*m2r2d2a2*(1.d0-kdkm1)*T1/X1_den*(1.d0-0.5d0*lld)
     &           +0.5d0*R*(ln-kdkm1*(T2+T1)/T1)
!
            B2=2.d0*m2r2d2a2*(1.d0-kdkm1)*T2/X2_den*(1.d0+0.5d0*lld)
     &           +0.5d0*R*(ln+kdkm1*(T2+T1)/T2)
!
!     residual
!     
            f=(m2r2d2a2*(X2*(1.d0+0.5d0*lld)-X1*(1.d0-0.5d0*lld))
     &           +R*T_moy*ln
     &           +b2/c2*(2*cp*A**2*(Tt2-T2)*X2_den-xflow**2*R**2*T2**2)
     &           +b1/c1*(2*cp*A**2*(Tt1-T1)*X1_den-xflow**2*R**2*T1**2))
!     
!     pressure node1
!     
            df(1)=(2.d0*m2r2d2a2*X1/pt1*(1.d0-0.5d0*lld)
     &           -R*T_moy/pt1
     &           +B1/C1*(4.d0*cp*A**2*(Tt1-T1)*pt1*X_T1dTt1))
!     
!     temperature node1
!     
            df(2)=(-2.d0*m2r2d2a2*(kdkm1/Tt1*X1)*(1.d0-0.5d0*lld)
     &           +r*kdkm1*T_moy/Tt1
     &          +b1/c1*(2*cp*A*A*X1_den*(1.d0-2.d0*kdkm1*(Tt1-T1)/Tt1)))
!     
!     mass flow
!     
            df(3)=(inv*xflow*R**2/a**2
     &           *(X2*(1.d0+0.5d0*lld)-X1*(1.d0-0.5d0*lld))
     &           +B2/C2*(-2.d0*inv*xflow*R*R*T2**2.d0)
     &           +B1/C1*(-2.d0*inv*xflow*R*R*T1**2.d0))
!     
!     pressure node2
!     
            df(4)=(-2*m2r2d2a2*X2/pt2*(1.d0+0.5d0*lld)
     &           +R*T_moy/pt2
     &           +B2/C2*(4.d0*cp*A*A*(Tt2-T2)*pt2*X_T2dTt2))
!     
!     temperature node2
!     
            df(5)=(2.d0*m2r2d2a2*(kdkm1/Tt2*X2)*(1.d0+0.5d0*lld)
     &           -r*kdkm1*T_moy/Tt2
     &          +b2/c2*(2*cp*A*A*X2_den*(1.d0-2.d0*kdkm1*(Tt2-T2)/Tt2)))
!
         else
!     
            pt1=pt2/pt2zpt1_c
            f=xflow*dsqrt(Tt1)/pt1-A*qred_crit
!
!     pressure node1
!     
            df(1)=-xflow*dsqrt(Tt1)/pt1**2
!     
!     temperature node1
!     
            df(2)=0.5d0*xflow/(pt1*dsqrt(Tt1))
!     
!     mass flow
!     
            df(3)=inv*dsqrt(Tt1)/pt1
!     
!     pressure node2
!     
            df(4)=0.d0
!     
!     temperature node2
!     
            df(5)=0.d0
!     
         endif
!     
      elseif(iflag.eq.3) then

         pi=4.d0*datan(1.d0)
         e=2.7182818d0
!     
         kappa=(cp/(cp-R))
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         kdkp1=kappa/kp1
!         
         index=ielprop(nelem)
         A=prop(index+1)
!
         d=prop(index+2)
!     
         l=prop(index+3)
         if(l.lt.0d0) then
            l_neg=l
            l=abs(l)
         else
            l_neg=l
         endif
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPIA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPII') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=int(prop(index+7))
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!     
         if(xflow.ge.0d0) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)
!     
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
!            
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
!     
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=-v(1,nodem)*iaxial
            Tt1=v(0,node2)+physcon(1)
            Tt2=v(0,node1)+physcon(1)
!
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
!     
            nodef(1)=node2
            nodef(2)=node2
            nodef(3)=nodem
            nodef(4)=node1
            nodef(5)=node1
         endif
!     
         pt2zpt1=pt2/pt1
!     
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1E-30) then
            kgas=0
            call dynamic_viscosity(kgas,T1,dvi)
         endif
         reynolds=dabs(xflow)*d/(dvi*a)
         if(reynolds.lt.100.d0) then
            reynolds= 100.d0
         endif
!     
!     definition of the friction coefficient for 2 phase flows and pure air
!     
!     Lockhart-Martinelli method
!
         if(lakon(nelem)(7:7).eq.'F') then
!     
            if((k_oil.lt.0).or.(k_oil.gt.12)) then
               write(*,*) '*ERROR:in gaspipe.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the type of oil is not defined'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Current calculation stops here'
               call exit(201)
            elseif(xflow_oil.eq.0) then
               write(*,*) '*WARNING:in gaspipe.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the oil mass flow rate is NULL'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Only pure air is considered'
               phi=1
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
            else
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
               lambda=lambda*phi
            endif
!     
!     for pure air
!     
         else
            phi=1.d0
            call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &           lambda)
         endif
!     
         call pt2zpt1_crit(pt2,pt1,Tt1,Tt2,lambda,kappa,r,l,d,A,iflag,
     &     inv,pt2zpt1_c,qred_crit,crit,qred_max1,icase)
!     
!     definition of the coefficients 
!     
         M1=dsqrt(2/km1*((Tt1/T1)-1))
         M2=dsqrt(2/km1*((Tt2/T2)-1))
            
         write(1,*) ''
         write(1,55) 'In line',int(nodem/100),' from node',node1,
     &' to node', node2,':   air massflow rate= ',xflow,'kg/s',
     &', oil massflow rate= ',xflow_oil,'kg/s'
 55      FORMAT(1X,A,I6.3,A,I6.3,A,I6.3,A,F9.6,A,A,F9.6,A)
! 
         if(inv.eq.1) then
            write(1,53)'       Inlet node ',node1,':    Tt1= ',Tt1,
     &           'K, Ts1= ',T1,'K, Pt1= ',Pt1/1E5,
     &           'Bar, M1= ',M1
            write(1,*)'             element W    ',set(nelem+numf)(1:20)
            write(1,57)'             eta=',dvi,'kg/(m*s), Re= '
     &           ,reynolds,', Phi= ',phi,', lambda= ',lambda,
     &           ', lambda*l/d= ',lambda*l/d,', zeta= ',phi*lambda*l/d
            write(1,53)'       Outlet node ',node2,'    Tt2= ',Tt2,
     &           'K, Ts2= ',T2,'K, Pt2= ',Pt2/1e5,
     &           'Bar, M2= ',M2
!     
         else if(inv.eq.-1) then
            write(1,53)'       Inlet node ',node2,':    Tt1= ',Tt1,
     &           'K, Ts1= ',T1,'K, Pt1= ',Pt1/1E5,
     &           'Bar, M1= ',M1
            write(1,*)'             element W    ',set(nelem+numf)(1:20)
            write(1,57)'             eta= ',dvi,'kg/(m*s), Re= '
     &           ,reynolds,' ,Phi= ',phi,', lambda= ',lambda,
     &           ', lambda*l/d= ',lambda*l/d,', zeta= ',phi*lambda*l/d
            write(1,53)'       Outlet node ',node1,'    Tt2= ',Tt2,
     &           'K, Ts2= ',T2,'K, Pt2=',Pt2/1e5,
     &           'Bar, M2= ',M2
         endif
      endif
      
 53   FORMAT(1X,A,I6.3,A,f6.1,A,f6.1,A,f9.5,A,f8.5)  
 57   FORMAT(1X,A,G11.4,A,G11.4,A,f8.5,A,f8.5,A,f8.5,A,f8.5)
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!
      return
      end
