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
      subroutine tt_calc(xflow,Tt,Pt,kappa,r,a,Ts,icase,iflag)
!     
!     this subroutine solves the implicit equation
!     f=xflow*dsqrt(Tt)/(a*Pt)-C*(TtdT)**expon*(Ttdt-1)**0.5d0
!     in order to find Tt when ts , xflow, pt and a are given
!
!     author: Yannick Muller
!     
      implicit none
!
      integer inv,icase,i,iflag
!     
      real*8 xflow,Tt,Pt,Ts,kappa,r,f,df,a,expon,Tt_old,C,TtzTs,
     &     deltaTt,TtzTs_crit, Qred_crit,Qred,h1,h2,h3,Tt_crit
!     
      expon=-0.5d0*(kappa+1.d0)/(kappa-1.d0)
!     
      C=dsqrt(2.d0/r*kappa/(kappa-1.d0))
!     
!     f=xflow*dsqrt(Tt)/(a*Pt)-C*(TtdT)**expon*(Ttdt-1)**0.5d0
!
!     df=-C*Ttdt**expon*(expon/Ts*(TtdT-1)**0.5d0
!     &     -0.5d0*TtdT/Ts*(TtdT-1.d0)**(-0.5d0))
!     
      Tt_old=Tt
!    
!
      if(xflow.lt.0d0) then 
         inv=-1
      else
         inv=1
      endif
!    
      if(dabs(xflow).le.1e-9) then
         Tt=Ts
         return
      endif
!
      Qred=abs(xflow)*dsqrt(Tt)/(a*Pt)
!
c      write(*,*) 'epsilon',(Qred/C)**2
!
!     optimised estimate of T static
!
      Tt=Ts*(1+(Qred**2/C**2))
!     
!     adiabatic
!     
      if(icase.eq.0) then
!
         TtzTs_crit=(kappa+1.d0)/2.d0
!         
!     isothermal
!
      else
!     
!         if(iflag.ne.3) then    
         TtzTs_crit=(1d0+(kappa-1.d0)/(2.d0*kappa))
         if(iflag.ne.3) then
            Tt_crit=(A*pt/(dabs(xflow))*C*Ttzts_crit**expon*
     &           (TtzTs_crit-1)**0.5d0)**2
!         if(iflag.ne.3) then
            if(Tt.gt.Tt_crit) then
               Tt=Tt_crit
            endif
         endif
!         Tt=Tt_crit
!         Qred=abs(xflow)*dsqrt(Tt)/(a*Pt)
!
      endif
!
      Qred_crit=C*(TtzTs_crit)**expon*(Ttzts_crit-1.d0)**0.5d0
!     
!     
      if(Qred.ge.Qred_crit) then
!     
         Tt=Ts*TtzTs_crit
!     
         return
!     
      endif
      i=0
!     
      do 
         i=i+1
         Ttzts=Tt/Ts
         h1=Ttzts-1.d0
         h2= dsqrt(h1)
         h3=Ttzts**expon
!     
         f=C*h2*h3
!     
         df=0.5*inv*xflow/(A*Pt*dsqrt(Tt))
     &        -1/Tt*C*h2*h3*(expon+0.5d0*(h1)**(-1))
!
         f=Qred-f
         deltaTt=-f/df
c         write(*,*) 'deltaTs=',deltaTs
!     
         Tt=Tt+deltaTt
c         write(*,*) 'Ts',Ts
!     
         if( (((dabs(Tt-Tt_old)/tt_old).le.1.E-8))
     &        .or.((dabs(Tt-Tt_old)).le.1.E-10) 
     &    .or.((f.le.1E-5).and.(deltaTt.lt.1E-3))) then
c            write(*,*) 'f=',f
c            write(*,*) 'Ts=',Ts
c            write(*,*) 'i',i
c$$$c            write(*,*) ''
c$$$            write(*,*) 'f',f
c$$$            write(*,*) ''
c$$$            write(*,*) 'Ts',Ts
c$$$            write(*,*) 'Tt',Tt
c$$$            write(*,*) ''
            Qred_crit=C*(TtzTs_crit)**expon*(Ttzts_crit-1.d0)**0.5d0
            Qred=abs(xflow)*dsqrt(Tt)/(a*Pt)
!     
            if((Qred.ge.Qred_crit).and.(iflag.eq.3)) then
!     
               Tt=Ts*TtzTs_crit
!     
            endif
            exit
         else if((i.gt.40)) then
            Tt=0.99*Ts*TtzTs_crit
c$$$            Tt=1.2*Ts
c$$$            write(*,*) 'Break'
c$$$            write(*,*) 'f',f
c$$$            write(*,*) 'Ts',Ts
c$$$            write(*,*) 'Tt',Tt
c$$$            write(*,*) ''
            exit
         endif
         Tt_old=Tt
      enddo
!     
C      write(*,*) 'end of ts_clac.f'
c      write(*,*) ''
      return
      end
      
      
      
