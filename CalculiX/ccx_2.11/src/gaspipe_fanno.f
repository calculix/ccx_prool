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
      subroutine gaspipe_fanno(node1,node2,nodem,nelem,lakon,kon,
     &        ipkon,nactdog,identity,ielprop,prop,iflag,v,xflow,f,
     &        nodef,idirf,df,cp,r,physcon,dvi,numf,set,
     &        shcon,nshcon,rhcon,nrhcon,ntmat_,co,vold,mi,ttime,time,
     &        iaxial)
!     
!     pipe with friction losses (Fanno Formulas) GAPF 
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
     &     inv,ipkon(*),kon(*),icase,kgas,k_oil
     &     ,nshcon(*),nrhcon(*),ntmat_,i,mi(*),nodea,nodeb,
     &     nodec,iaxial
!
      real*8 prop(*),v(0:mi(2),*),xflow,f,df(5),kappa,R,a,d,l,
     &     p1,p2,T1,T2,Tt1,Tt2,pt1,pt2,cp,physcon(3),p2p1,km1,dvi,
     &     kp1,kdkm1,reynolds,pi,e,lambda,lld,kdkp1,T2dTt2,
     &     T1dTt1,X_t1dTt1,X_t2dTt2,X2_den,X1_den,
     &     X1,X2,B1,B2,C1,C2,tdkp1,ttime,time,
     &     pt2zpt1,ks,form_fact,xflow_oil,Tt1dT1,Tt2dT2,
     &     Pt2zPt1_c,qred_crit,l_neg,Qred,
     &     expon1,expon2,cte,term1,term2,term3,term4,term5,term6,
     &     term,phi,M1,M2,qred2,qred1,qred_max1,qred_crit_out,co(3,*),
     &     shcon(0:3,ntmat_,*),rhcon(0:1,ntmat_,*),vold(0:mi(2),*),
     &     radius,initial_radius,l_initial
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
         pi=4.d0*datan(1.d0)
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
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=nint(prop(index+7))
!
         if((lakon(nelem)(2:6).eq.'GAPFF').and.
     &        (lakon(nelem)(2:7).ne.'GAPFF2')) then
!            
            icase=0
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
c            iaxial=nint(prop(index+3))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
!
            initial_radius=dsqrt((co(1,nodeb)-co(1,nodea))**2)
!
c            if(iaxial.ne.0) then
c              A=pi*radius**2/iaxial
c            else
              A=pi*radius**2
c            endif
            d=2*radius
            l=prop(index+4)
            if(l.lt.0d0) then
               l_neg=l
               l=abs(l)
            else
               l_neg=l
            endif
            ks=prop(index+5)
            form_fact=prop(index+6)
            xflow_oil=prop(index+7)
            k_oil=nint(prop(index+8))
!
         elseif (lakon(nelem)(2:7).eq.'GAPFF2') then
            write(*,*) nelem,lakon(nelem)(1:6)
            icase=0
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            nodec=nint(prop(index+3))
c            iaxial=nint(prop(index+4))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            initial_radius=dsqrt((co(1,nodeb)-co(1,nodea))**2)
            d=2*radius
c           if(iaxial.ne.0) then
c              A=pi*radius**2/iaxial
c            else
               A=pi*radius**2
c            endif
            l_initial=dsqrt((co(2,nodec)-co(2,nodeb))**2)
            l=dsqrt((co(2,nodec)+vold(2,nodec)-
     &           co(2,nodeb)-vold(2,nodeb))**2)
            if(l.lt.0d0) then
               l_neg=l
               l=abs(l)
            else
               l_neg=l
            endif
            ks=prop(index+5)
            form_fact=prop(index+6)           
            xflow_oil=prop(index+7)
            k_oil=nint(prop(index+8))
         endif
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!
         if(pt1.ge.pt2) then
            inv=1
            Tt1=v(0,node1)+physcon(1)
            Tt2=v(0,node2)+physcon(1)
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            Tt1=v(0,node2)+physcon(1)
            Tt2=v(0,node1)+physcon(1)
         endif
!
         p2p1=pt2/pt1
         km1=kappa-1.d0
         kp1=kappa+1.d0
         kdkm1=kappa/km1
         tdkp1=2.d0/kp1
         C2=tdkp1**kdkm1
!        
!     incompressible flow
         xflow=inv*A*dsqrt(d/l*2*Pt1/(R*Tt1)*(pt1-pt2))
         if(p2p1.gt.C2) then
           xflow=inv*pt1*a*dsqrt(2.d0*kdkm1*p2p1**(2.d0/kappa)
     &           *(1.d0-p2p1**(1.d0/kdkm1))/r)/dsqrt(Tt1)
         else
            xflow=inv*pt1*a*dsqrt(kappa/r)*tdkp1**(kp1/(2.d0*km1))/
     &           dsqrt(Tt1)
         endif
!
!     calculation of the dynamic viscosity 
!     
         if(dabs(dvi).lt.1E-30) then
            kgas=0
            call dynamic_viscosity(kgas,Tt1,dvi)
         endif  
!
         reynolds=dabs(xflow)*d/(dvi*a)
!
         call friction_coefficient(l_neg,d,ks,reynolds,form_fact,lambda)
         xflow=inv*A*dsqrt(d/(lambda*l)*2*Pt1/(R*Tt1)*(pt1-pt2))
!
         call pt2zpt1_crit(pt2,pt1,Tt1,Tt2,lambda,kappa,r,l,d,A,iflag,
     &     inv,pt2zpt1_c,qred_crit,crit,qred_max1,icase)
!
         Qred=dabs(xflow)*dsqrt(Tt1)/(A*pt1)
!
         if (crit) then
            xflow=0.5*inv*Qred_crit*Pt1*A/dsqrt(Tt1)
            if(icase.eq.1) then
!     
               call ts_calc(xflow,Tt1,pt1,kappa,r,a,T1,icase)
               if (inv.eq.1) then 
                  v(3,node1)=T1
                  v(3,node2)=T1
                  if(nactdog(0,node2).eq.1) then
                     v(0,node2)=T1*(1.d0+km1/(2*kappa))
                  endif
               else
                  v(3,node2)=T1
                  v(3,node1)=T1
                  if(nactdog(0,node1).eq.1) then
                     v(0,node1)=T1*(1.d0+km1/(2*kappa)) 
                  endif
               endif
            endif
         elseif(Qred.gt.Qred_crit) then
            xflow=0.5*inv*Qred_crit*pt1*A/dsqrt(Tt1)
         else
            xflow=inv*Qred*pt1*A/dsqrt(Tt1)
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
!
         l=prop(index+3)
         if(l.lt.0d0) then
            l_neg=l
            l=abs(l)
         else
            l_neg=l
         endif
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=nint(prop(index+7))
!
                  if((lakon(nelem)(2:6).eq.'GAPFF').and.
     &        (lakon(nelem)(2:7).ne.'GAPFF2')) then
            icase=0
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
c            iaxial=nint(prop(index+3))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            initial_radius=dsqrt((co(1,nodeb)-co(1,nodea))**2)
            d=2*radius
c           if(iaxial.ne.0) then
c              A=pi*radius**2/iaxial
c            else
               A=pi*radius**2
c            endif
            l=prop(index+4)
            if(l.lt.0d0) then
               l_neg=l
               l=abs(l)
            else
               l_neg=l
            endif
            ks=prop(index+5)
            form_fact=prop(index+6)          
            xflow_oil=prop(index+7)
            k_oil=nint(prop(index+8))
!
         elseif (lakon(nelem)(2:7).eq.'GAPFF2') then
            icase=0
            nodea=nint(prop(index+1))
            nodeb=nint(prop(index+2))
            nodec=nint(prop(index+3))
c            iaxial=nint(prop(index+4))
            radius=dsqrt((co(1,nodeb)+vold(1,nodeb)-
     &           co(1,nodea)-vold(1,nodea))**2)
            initial_radius=dsqrt((co(1,nodeb)-co(1,nodea))**2)
            d=2*radius
c           if(iaxial.ne.0) then
c              A=pi*radius**2/iaxial
c            else
               A=pi*radius**2
c            endif
            l_initial=dsqrt((co(2,nodec)-co(2,nodeb))**2)
            l=-dsqrt((co(2,nodec)+vold(2,nodec)-
     &           co(2,nodeb)-vold(2,nodeb))**2)
            if(l.lt.0d0) then
               l_neg=l
               l=abs(l)
            else
               l_neg=l
            endif
            ks=prop(index+5)
            form_fact=prop(index+6)         
            xflow_oil=prop(index+7)
            k_oil=nint(prop(index+8))
         endif
!
         pt1=v(2,node1)
         pt2=v(2,node2)
         xflow=v(1,nodem)*iaxial
!
         if((pt1.gt.pt2).or.(xflow.ge.0d0)) then
            inv=1
            Tt1=v(0,node1)+physcon(1)
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
            if(icase.eq.0) then
               Tt2=Tt1
               call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
            else
               t2=t1
               Tt2=v(0,node2)+physcon(1)
            endif
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
            if(icase.eq.0) then
               Tt2=Tt1
            else
               Tt2=v(0,node1)+physcon(1)
            endif
!
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
           if(dabs(dvi).lt.1E-30) then
              kgas=0
              call dynamic_viscosity(kgas,T1,dvi)
           endif        
!           
           reynolds=dabs(xflow)*d/(dvi*a)
!
           if(reynolds.lt.1) then
              reynolds = 1.d0
           endif
!
!     definition of the friction coefficient for 2 phase flows and pure air
!     
!     Friedel's Method
         if(lakon(nelem)(7:7).eq.'F') then
!
            if((k_oil.lt.0).or.(k_oil.gt.12)) then
               write(*,*) '*ERROR:in gaspipe.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the type of oil is not defined'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Current calculation stops here'
               call exit(201)
            elseif(xflow_oil.eq.0.d0) then
               write(*,*) '*WARNING:in gaspipe.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the oil mass flow rate is NULL'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Only pure air is considered'
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
            else
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!
               lambda=lambda*phi
!
            endif
!     
!     Alber's Method
!     
         elseif (lakon(nelem)(7:7).eq.'A') then
            if((k_oil.lt.0).or.(k_oil.gt.12)) then
               write(*,*) '*ERROR:in gaspipe_fanno.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the type of oil is not defined'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Current calculation stops here'
               call exit(201)
            elseif(xflow_oil.eq.0) then
               write(*,*) '*WARNING:in gaspipe_fanno.f'
               write(*,*) ' using two phase flow'
               write(*,*) ' the oil mass flow rate is NULL'
               write(*,*) ' check element ',nelem,' definition'
               write(*,*) ' Only pure air is considered'
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
            else
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!     
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
!
               lambda=lambda*phi
!
            endif
!     
!     for pure air
!     
         else
!     
            phi=1.d0
            call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &           lambda)
         endif
!     
         call pt2zpt1_crit(pt2,pt1,Tt1,Tt2,lambda,kappa,r,l,d,A,iflag,
     &     inv,pt2zpt1_c,qred_crit,crit,qred_max1,icase)
!
         Qred1=xflow*dsqrt(Tt1)/(A*Pt1)
!
         if(dabs(xflow)*dsqrt(Tt1)/(A*Pt1).gt.qred_max1) then
            crit=.true.
         endif
!
         Qred2=xflow*dsqrt(Tt2)/(A*Pt2)
         if(icase.eq.0) then
            qred_crit_out=dsqrt(kappa/R)*(2/(kappa+1))**(0.5d0*
     &           (kappa+1)/(kappa-1))
         else
            qred_crit_out=R**(-0.5d0)*(2/(kappa+1))**(0.5d0*
     &           (kappa+1)/(kappa-1))
         endif
!     
!     definition of the coefficients
!
         lld=lambda*l/d
!
         M2=dsqrt(2/km1*((Tt2/T2)-1))
         if(icase.eq.0) then
            if((M2.lt.1)) then
               crit=.false.
               if((M2.ge.1.d0).or.(dabs(M2-1).lt.1E-5)) then
                  pt2=pt1*pt2zpt1_c
               endif
            endif
         elseif (icase.eq.1) then
            if(M2.lt.1/dsqrt(kappa)) then
               crit=.false.
            else
               crit=.true.
            endif
         endif
!
!     adiabatic case
!     
         if(icase.eq.0) then
!
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
!            C1=2.d0*cp*A**2*X1_den*(-1.d0+2.d0*kdkm1*T1dTt1)
!     &           -2.d0*xflow**2*R**2*T1
            C1=2.d0*cp*A**2*X1_den*(-1.d0+2.d0*kdkm1*(T1dTt1-1))
     &           -2.d0*xflow**2*R**2*T1
!     
!            C2=2.d0*cp*A**2*X2_den*(-1.d0+2.d0*kdkm1*T2dTt2)
!     &           -2.d0*xflow**2*R**2*T2
            C2=2.d0*cp*A**2*X2_den*(-1.d0+2.d0*kdkm1*(T2dTt2-1))
     &           -2.d0*xflow**2*R**2*T2
!     
            expon1=(kappa+1)/km1
            expon2=2*kappa/(km1)
!     
            cte=0.5d0*(kappa+1)/kappa
!
            term1=pt1**2*T1**expon1*Tt1**(-expon2)*A**2
!
            if(.not.crit) then
               term1=pt1**2*T1**expon1*Tt1**(-expon2)*A**2
               term2=pt2**2*T2**expon1*Tt2**(-expon2)*A**2        
!     
!     simplified version
               term3=Tt2dT2
               term4=Tt1dT1
!     
               term5=T1**(expon1)*Tt1**(-expon2)*(pt1**2)               
               term6=T2**(expon1)*Tt2**(-expon2)*(pt2**2)               
!
               B1=1/(R*xflow**2)*term1*expon1/T1
     &              +cte*(-(2/km1)*1/T1)
!     
               B2=1/(R*xflow**2)*term2*(-expon1/T2)
     &              +cte*(2/km1*1/T2)
!     
!     residual
!    
!     Simplified version
!     
               f=1/(R*xflow**2)*(term1-term2)
     &              +cte*(log(term3)-log(term4)-log(term5)+log(term6))
     &              -lld 
     &              +b2/c2*(2*cp*A**2*(Tt2-T2)
     &              *X2_den-xflow**2*R**2*T2**2)
     &              +b1/c1*(2*cp*A**2*(Tt1-T1)
     &              *X1_den-xflow**2*R**2*T1**2)
!     
!     pressure node1
!     
               df(1)=1/(R*xflow**2)*(term1*2/pt1)
     &              +cte*(-2/pt1)
     &              +B1/C1*(4.d0*cp*A**2*(Tt1-T1)*pt1*X_T1dTt1)
!     
!     temperature node1
!     
               df(2)=1/(R*xflow**2)*term1*(-expon2)/Tt1
     &              +cte*(expon1*1/Tt1)
     &              +b1/c1*(2*cp*A**2*X1_den
     &              *(1.d0-2.d0*kdkm1*(Tt1-T1)/Tt1))
!     
!     mass flow
!     
               df(3)=-2.d0/(R*(inv*xflow)**3)*(term1-term2)
     &              +B2/C2*(-2.d0*inv*xflow*R*R*T2**2.d0)
     &              +B1/C1*(-2.d0*inv*xflow*R*R*T1**2.d0)
!     
!     pressure node2
!     
               df(4)=1/(R*xflow**2)*(-term2*2/pt2)
     &              +cte*(2/pt2)
     &              +B2/C2*(4.d0*cp*A**2*(Tt2-T2)*pt2*X_T2dTt2)
!     
!     temperature node2
!     
               df(5)=1/(R*xflow**2)*term2*(expon2/Tt2)
     &              +cte*(-expon1*1/Tt2)
     &              +b2/c2*(2*cp*A**2*X2_den
     &              *(1.d0-2.d0*kdkm1*(Tt2-T2)/Tt2))
!               
            else
!
               term=kappa*term1/(xflow**2*R)
               B1=expon1*1/T1*(1/kappa*term-1)+cte*1/T1
!     f=1/kappa*(term1-1)+cte*(log(T1dTt1)-log(2/kp1*term))
               f=1/kappa*(term-1)+cte*(log(T1dTt1)-log(2/kp1*term))
     &              -lld
     &              +b1/c1*(2*cp*A**2*(Tt1-T1)
     &              *X1_den-xflow**2*R**2*T1**2)
!     
!     pressure node1
!     
               df(1)=2/pt1*(1/kappa*term-cte)
     &              +B1/C1*(4.d0*cp*A**2*(Tt1-T1)*pt1*X_T1dTt1)
!     
!     temperature node1
!     
               df(2)=expon2*1/Tt1*(-1/kappa*term+1)-cte*1/Tt1
     &              +b1/c1*(2*cp*A**2*X1_den
     &              *(1.d0-2.d0*kdkm1*(Tt1-T1)/Tt1))
!     
!     mass flow
!     
               df(3)=2.d0/(inv*xflow)*(-term/kappa+cte)
     &              +B1/C1*(-2.d0*inv*xflow*R*R*T1**2.d0)
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
!     isothermal icase
!     
         elseif(icase.eq.1) then
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
               C1=2.d0*cp*A**2*X1_den*(1.d0-2.d0*kdkm1*(Tt1dT1-1.d0))
     &              +2.d0*xflow**2*R**2*T1
!     
               C2=2.d0*cp*A**2*X2_den*(1.d0-2.d0*kdkm1*(Tt2dT2-1.d0))
     &              +2.d0*xflow**2*R**2*T2
!     
               expon1=(kappa+1)/km1
               expon2=2*kappa/(kappa-1)
!     
               cte=0.5d0*(kappa+1)/kappa
!     
               term1=pt1**2*T1**expon1*Tt1**(-expon2)*A**2
               term2=pt2**2*T2**expon1*Tt2**(-expon2)*A**2
!     
               term5=T1**(expon1)*Tt1**(-expon2)*(pt1**2*A**2)
               term6=T2**(expon1)*Tt2**(-expon2)*(pt2**2*A**2)
!
            if(.not.crit) then
               B1=1/(R*xflow**2)*term1*expon1/T1
     &              -expon1/T1   
!     
               B2=1/(R*xflow**2)*term2*(-expon1/T2)
     &              +expon1/T2  
!     
!     Simplified version
!     
               f=1/(R*xflow**2)*(term1-term2)
     &              +(-log(term5)+log(term6))            
     &              -lld 
     &              +b2/c2*(2*cp*A**2*(Tt2-T2)
     &              *X2_den-xflow**2*R**2*T2**2)
     &              +b1/c1*(2*cp*A**2*(Tt1-T1)
     &              *X1_den-xflow**2*R**2*T1**2)
!     
!     pressure node1
!     
               df(1)=1/(R*xflow**2)*(term1*2/pt1)
     &              +(-(2/pt1))
     &              +B1/C1*(4.d0*cp*A**2*(Tt1-T1)*pt1*X_T1dTt1)
!     
!     temperature node1
!     
               df(2)=1/(R*xflow**2)*term1*(-expon2)/Tt1
     &              +(expon2/Tt1)            
     &              +b1/c1*(2*cp*A**2*X1_den
     &              *(1.d0-2.d0*kdkm1*(Tt1-T1)/Tt1))
!     
!     mass flow
!     
               df(3)=-2.d0/(R*xflow**3)*(term1-term2)
     &              +B2/C2*(-2.d0*inv*xflow*R*R*T2**2.d0)
     &              +B1/C1*(-2.d0*inv*xflow*R*R*T1**2.d0)
!     
!     pressure node2
!     
               df(4)=1/(R*xflow**2)*(-term2*2/pt2)
     &              +(2/pt2)
     &              +B2/C2*(4.d0*cp*A**2*(Tt2-T2)*pt2*X_T2dTt2)
!     
!     
!     temperature node2
!     
               df(5)=1/(R*xflow**2)*term2*(expon2/Tt2)
     &              +(-expon2/Tt2)            
     &              +b2/c2*(2*cp*A**2*X2_den
     &              *(1.d0-2.d0*kdkm1*(Tt2-T2)/Tt2))
               
            else
               term=term1/(xflow**2*R)
               B1=expon1/T1*(term-1)
!     alternate critical equation
!  
               f=term-1-log(term)            
     &              -lld 
     &              +b1/c1*(2*cp*A**2*(Tt1-T1)
     &              *X1_den-xflow**2*R**2*T1**2)
!     
!     pressure node1
!     
               df(1)=2/pt1*(term-1)
     &              +B1/C1*(4.d0*cp*A**2*(Tt1-T1)*pt1*X_T1dTt1)
!     
!     temperature node1
!     
               df(2)=expon2/Tt1*(-term+1)            
     &              +b1/c1*(2*cp*A**2*X1_den
     &              *(1.d0-2.d0*kdkm1*(Tt1-T1)/Tt1))
!     
!     mass flow
!     
               df(3)=2/xflow*(-term+1)
     &              +B1/C1*(-2.d0*inv*xflow*R*R*T1**2.d0)
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
         endif
!
!     output
!
      elseif(iflag.eq.3) then
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
!
         lambda=0.5
!
         if(l.lt.0d0) then
            l_neg=l
            l=abs(l)
         else
            l_neg=l
         endif
         ks=prop(index+4)
         if(lakon(nelem)(2:6).eq.'GAPFA') then
            icase=0
         elseif(lakon(nelem)(2:6).eq.'GAPFI') then
            icase=1
         endif
         form_fact=prop(index+5)
         xflow_oil=prop(index+6)
         k_oil=nint(prop(index+7))
!
         pt1=v(2,node1)
         pt2=v(2,node2)
!     
         if(xflow.ge.0d0) then
            inv=1
            xflow=v(1,nodem)*iaxial
            Tt1=v(0,node1)+physcon(1)
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
            if(icase.eq.0) then
               Tt2=Tt1
               call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
            else
               T2=T1
               Tt2=v(0,node2)+physcon(1)
               call tt_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase,iflag)
            endif
!      
!            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
!!
!            call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
!     
         else
            inv=-1
            pt1=v(2,node2)
            pt2=v(2,node1)
            xflow=v(1,nodem)*iaxial
!           
            Tt1=v(0,node2)+physcon(1)
            if(icase.eq.0) then
               Tt2=Tt1
            else
               Tt2=v(0,node1)+physcon(1)
            endif
!               
            call ts_calc(xflow,Tt1,Pt1,kappa,r,a,T1,icase)
            call ts_calc(xflow,Tt2,Pt2,kappa,r,a,T2,icase)
!     
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
! 
         reynolds=dabs(xflow)*d/(dvi*a)
!
        if(reynolds.lt.1.d0) then
            reynolds= 1.d0
         endif
!     
!     definition of the friction coefficient for 2 phase flows and pure air
!     
!     Friedel's Method
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
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
            else
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
!
            endif
!     
         elseif (lakon(nelem)(7:7).eq.'A') then
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
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
            else
               call two_phase_flow(Tt1,pt1,T1,Tt2,pt2,T2,xflow,
     &              xflow_oil,nelem,lakon,kon,ipkon,ielprop,prop,
     &              v,dvi,cp,r,k_oil,phi,lambda,nshcon,nrhcon,
     &              shcon,rhcon,ntmat_,mi)
!     
               call friction_coefficient(l_neg,d,ks,reynolds,form_fact,
     &              lambda)
!
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
!     
            write(1,*) ''
            write(1,55) ' from node ',node1,
     &           ' to node ', node2,':   air massflow rate = ',xflow,
     &           ' , oil massflow rate = ',xflow_oil
!     
            if(inv.eq.1) then
               write(1,53)'       Inlet node ',node1,' :    Tt1 = ',Tt1,
     &              ' , Ts1 = ',T1,'  , Pt1 = ',Pt1,
     &              ' ,M1 = ',M1
               write(1,*)'             Element ',nelem,lakon(nelem)
               write(1,57)'             eta = ',dvi,' kg/(m*s) , Re = '
     &              ,reynolds
               write(1,58)'             PHI = ',phi,' , LAMBDA = ',
     &              lambda,
     &              ', LAMBDA*l/d = ',lambda*l/d,' , ZETA_PHI = ',
     &              phi*lambda*l/d
               write(1,53)'       Outlet node ',node2,' :    Tt2 = ',
     &              Tt2,
     &              '  , Ts2 = ',T2,'  , Pt2 = ',Pt2,
     &              ' , M2 = ',M2
!    
            else if(inv.eq.-1) then
               write(1,53)'       Inlet node ',node2,':    Tt1= ',Tt1,
     &              ' , Ts1= ',T1,' , Pt1= ',Pt1,
     &              ' , M1= ',M1
               write(1,*)'             Element ',nelem,lakon(nelem)
               write(1,57)'             eta = ',dvi,' kg/(m*s) , Re = '
     &              ,reynolds
               write(1,58)'             PHI = ',phi,' , LAMBDA = ',
     &              lambda,
     &              ', LAMBDA*l/d = ',lambda*l/d,' , ZETA_PHI = ',
     &              phi*lambda*l/d
               write(1,53)'       Outlet node ',node1,' :    Tt2 = ',
     &              Tt2,
     &              '  , Ts2 = ',T2,'  , Pt2 =',Pt2,
     &              ' , M2 = ',M2
            endif
      endif
!     
 55   format(1X,a,i6,a,i6,a,e11.4,a,a,e11.4,a)
 53   format(1X,a,i6,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
 57   format(1X,a,e11.4,a,e11.4)
 58   format(1X,a,e11.4,a,e11.4,a,e11.4,a,e11.4)
 77   format(3x,a,e11.1,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4
     &     ,a,e11.4)	 
 78   format(a,i4,a,a,a)
 79   format(3X,a,i4,a,i4,a,i4)
 80   format(3x,a,e11.4,a,e11.4,a,e11.4)
!     
      xflow=xflow/iaxial
      df(3)=df(3)*iaxial
!     
      return
      end
