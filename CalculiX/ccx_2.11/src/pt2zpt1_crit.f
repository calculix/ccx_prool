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
!     calculate the maximal admissible pressure ratio pt2/pt1
!
!     1) assuming M2=1 for adiabatic respectively M2=1/dsqrt(kappa) for isotherm pipe choking 
!      M1 is calculated iteratively using a dichotomy scheme
!
!     2)the ratio of the critical pressure ratio  Qred_1/Qred_2crit=Pt2/Pt1
!     =D(M1)/D(M2_crit)is computed [D(M)=M*(1+0.5*(kappa-1)*M)**(-0.5*(kappa+1)/(kappa-1))]
!
!     author: Yannick Muller
!   
      subroutine pt2zpt1_crit(pt2,pt1,Tt1,Tt2,lambda,kappa,r,l,d,A,
     &     iflag,inv,pt2zpt1_c,qred_crit,crit,qred_max1,icase)
!     
      implicit none
!
      logical crit
!
      integer iflag,inv,icase,i
!     
      real*8 pt2,pt1,lambda,kappa,l,d,M1,pt2zpt1,pt2zpt1_c,
     &     km1,kp1,km1zk,kp1zk,Tt1,Tt2,r,A,
     &     xflow_crit,qred_crit,f1,f2,f3,m1_ac,m1_min,m1_max,
     &     expon1,qred_max1,lld
!
!     useful variables and constants
! 
      km1=kappa-1.d0
      kp1=kappa+1.d0
      km1zk=km1/kappa
      kp1zk=kp1/kappa
      lld=lambda*l/d
      expon1=-0.5d0*kp1/km1
!
!     adiabatic case
!
      if(icase.eq.0) then
!     
!     computing M1 using dichotomy method
!     
         i=1
          m1_max=1
          m1_min=0.001d0
          do
             i=i+1
             m1_ac=(m1_min+m1_max)*0.5d0
!
             f1=(1.d0-M1_min**2)*(kappa*M1_min**2)**(-1)
     &            +0.5d0*kp1zk*log((0.5d0*kp1)*M1_min**2
     &            *(1+0.5d0*km1*M1_min**2)**(-1))-lld
!     
             f2=(1.d0-M1_ac**2)*(kappa*M1_ac**2)**(-1)
     &            +0.5d0*kp1zk*log((0.5d0*kp1)*M1_ac**2
     &            *(1+0.5d0*km1*M1_ac**2)**(-1))-lld
!     
             f3=(1.d0-M1_max**2)*(kappa*M1_max**2)**(-1)
     &            +0.5d0*kp1zk*log((0.5d0*kp1)*M1_max**2
     &            *(1+0.5d0*km1*M1_max**2)**(-1))-lld
!     
             if(abs(f2).le.1E-6) then
                M1=m1_ac
                exit
             endif
             if(i.gt.50) then
                M1=M1_ac
                exit
             endif
!     
             if((f3.gt.f2).and.(f2.ge.f1)) then
                if((f1.lt.0d0).and.(f2.lt.0d0)) then
                   m1_min=m1_ac
                else
                   m1_max=m1_ac
                endif
             elseif((f3.lt.f2).and.(f2.le.f1)) then  
                if((f3.lt.0d0).and.(f2.lt.0d0) )then
                   m1_max=m1_ac
                else
                   m1_min=m1_ac
                endif
             endif
          enddo
!
          Pt2zpt1_c=M1*(0.5d0*kp1)**(0.5*kp1/km1)
     &         *(1+0.5d0*km1*M1**2)**(-0.5d0*kp1/km1)
!     
!     isotherm case
!
       elseif (icase.eq.1) then
!     
!     computing M1 using dichotomy method for choked conditions M2=1/dsqrt(kappa)
!     (1.d0-kappa*M1**2)/(kappa*M1**2)+log(kappa*M1**2)-lambda*l/d=0
!     
          m1_max=1/dsqrt(kappa)
          m1_min=0.1d0
          i=1
          do
             i=i+1
             m1_ac=(m1_min+m1_max)*0.5d0
!     
             f1=(1.d0-kappa*M1_min**2)/(kappa*M1_min**2)
     &            +log(kappa*M1_min**2)-lambda*l/d
!     
             f2=(1.d0-kappa*M1_ac**2)/(kappa*M1_ac**2)
     &            +log(kappa*M1_ac**2)-lambda*l/d
!     
             f3=(1.d0-kappa*M1_max**2)/(kappa*M1_max**2)
     &            +log(kappa*M1_max**2)-lambda*l/d
!     
             if((abs(f2).le.1E-5).or.(i.ge.50)) then
                M1=m1_ac
                exit
             endif
!     
             if((f3.gt.f2).and.(f2.ge.f1)) then
                if((f1.lt.0d0).and.(f2.lt.0d0)) then
                   m1_min=m1_ac
                else
                   m1_max=m1_ac
                endif
             elseif((f3.lt.f2).and.(f2.le.f1)) then  
                if((f3.lt.0d0).and.(f2.lt.0d0) )then
                   m1_max=m1_ac
                else
                   m1_min=m1_ac
                endif
             endif
          enddo
!     
!        computing the critical pressure ratio in the isothermal case
!     pt=A*dsqrt(kappa)/(xflow*dsqrt(kappa Tt))*M*(1+0.5d0*(kappa-1)M**2)**(-0.5d0*(kappa+1)/(kappa-1))
!     and forming the pressure ratio between inlet and outlet(choked)
!     
          Pt2zPt1_c=dsqrt(Tt2/Tt1)*M1*dsqrt(kappa)*((1+0.5d0*km1/kappa)
     &         /(1+0.5d0*km1*M1**2))**(0.5d0*(kappa+1)/km1)
!     
       endif
!     
       pt2zpt1=pt2/pt1
       if(Pt2zPt1.le.Pt2zPt1_c) then
          crit=.true.
       endif
!     
       if (iflag.eq.1) then
          xflow_crit=inv*M1*Pt1*A/dsqrt(Tt1)*dsqrt(kappa/r)
     &         *(1+0.5d0*km1*M1**2)**(-0.5d0*kp1/km1) 
       elseif(iflag.eq.2) then
             qred_max1=M1*dsqrt(kappa/r)
     &            *(1+0.5d0*km1*M1**2)**(-0.5d0*kp1/km1)
       endif
!     
       Qred_crit=M1*dsqrt(kappa/r)
     &      *(1+0.5d0*km1*M1**2)**(-0.5d0*kp1/km1)
!     
      return
      end      
      
      
