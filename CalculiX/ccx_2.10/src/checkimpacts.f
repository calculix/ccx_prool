!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine checkimpacts(ne,neini,temax,sizemaxinc,
     & energyref,tmin,tper,idivergence,
     & iforceincsize,istab,dtheta,enres,energy,energyini,allwk,
     & allwkini,dampwk,dampwkini,emax,mortar,maxdecay,enetoll)
!     
!     # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
!     Routine that contains the implementation of the logic to
!       rule the increment size during contact conditions.
!     The values of tolerances have been tested for the ball
!       model, the two sliding blades (simplified model of 
!       blade from Mr. Wittig) and the real blade model.
!     Friction has not been tested deeply.
! 
!     Main variables and meaning
!
!     sizemaxinc    : maximum size of the current increment
!     iforceincsize : flag to set "dtheta=sizemaxinc" in 
!                       checkconvergence.
!     cstate        : vector containing contact data
!     temax            : max. natural period of oscillation 
!     icase         : flag to print debug informations
!     emax          : maximum energy of the system over the ti-
!                       me history
!     enres         : energy residual before contact (or re-
!                       adapted)
!     delta         : eneres normalized
!     fact          : factor to set sizemaxinc according to the  
!                       contact formulation
!     stab_th       : \hat{r}_{e}(t_n) -> mod belytschko before
!                       contact (or initial value). This is the 
!                       stability threshold
!     delta_r_hat   : \hat{r}_{e}(t) - \hat{r}_{e}(t-1) -> varia-
!                       tion of the modified belitschko criterion 
!                       used to control jumps
!     r_hat         : \hat{r}_{e}(t) actual value of the modified  
!                       belitschko criterion
! 
!     Proposed by Matteo Pacher
! 
!     # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
!     
      implicit none
!     
      integer idivergence,
     & iforceincsize,ne,neini,istab,mortar
!     
      real*8 temax,energyref,sizemaxinc,tmin,tper,dtheta,
     & delta_r_hat,r_hat,delta,allwk,allwkini,energy(4),
     & energyini(4),dampwk,dampwkini,emax,fact,
     & r_hat_bc,maxdecay,enres,enetoll
!     
      intent(in) ne,neini,tmin,tper,energyref,dtheta,
     & temax,mortar,allwk,allwkini,energy,energyini,dampwk,
     & dampwkini,emax
!     
      intent(out) sizemaxinc,idivergence,iforceincsize,
     & istab
!     
!     Initialization
!     
c      icase=0
      iforceincsize=0
!     
!     Adaption of the energy residual (absolute/relative check)
!     
      if((enres.ne.0.0).and.(dabs(enres).lt.(enetoll/4.0)))then
         delta=enres*emax
      else
         delta=enres
      endif
!     
      if(mortar.eq.0)then
         fact=10.0
      else
         fact=1.0
      endif
!     
!     Compute thresholds and energy values
!     
      delta_r_hat=energy(1)+energy(2)+energy(3)+energy(4)-allwk-
     &     dampwk-(energyini(1)+energyini(2)+energyini(3)+
     &     energyini(4)-allwkini-dampwkini)
!
      if(emax.le.0.0)then
!     
!     No kinetic energy in the structure: energyref is the internal energy
!     this happens at the beginning of the calculation
!     
         r_hat=(energy(1)+energy(2)+energy(3)+energy(4)-allwk-
     &        dampwk-energyref)/energyref
         delta_r_hat=delta_r_hat/energyref
         r_hat_bc=delta/energyref
      else
         r_hat=(energy(1)+energy(2)+energy(3)+energy(4)-allwk-
     &        dampwk-energyref)/emax
         delta_r_hat=delta_r_hat/emax
         r_hat_bc=delta/emax
      endif
!     
!     Logic to adapt the increment size
!     
      if(mortar.eq.0)then
!     
!     Energy conservation rules for NODE TO SURFACE penalty contact
!     
         if((delta_r_hat.lt.(-0.008)).and.(ne.ge.neini))then
!     
!     Impact (or too high variation during pers. contact)
!     delta_r_hat = r_hat-r_hat_ini
!     
            idivergence=1
            sizemaxinc=dtheta*0.25
            iforceincsize=1
c            icase=1
         elseif((r_hat-r_hat_bc.gt.0.0025).and.(ne.le.neini))then
!     
!     Rebound (or too high variation during pers. contact)
!     r_hat_bc is r_hat before contact
!     
            idivergence=1
            sizemaxinc=dtheta*0.5
            iforceincsize=1
c            icase=3
         else
!     
!     Persistent Contact
!     
c            icase=2
            if(r_hat.gt.(-0.9*maxdecay))then
               sizemaxinc=max(fact*temax/tper,1.01*dtheta)
               sizemaxinc=min(sizemaxinc,100.0*temax/tper)
            else
               sizemaxinc=max(temax/tper/10.0,0.5*dtheta)
               istab=1
            endif
            
         endif
!     
      elseif(mortar.eq.1)then
!     
!     Energy conservation rules for SURFACE TO SURFACE penalty contact
!     
         if((delta_r_hat.lt.(-0.008)).and.(ne.ge.neini))then
!     
!     Impact (or too high variation during pers. contact)
!     delta_r_hat = r_hat-r_hat_ini
!     
            idivergence=1
            sizemaxinc=dtheta*0.25
            iforceincsize=1
c            icase=1
!     
         elseif((r_hat-r_hat_bc.gt.0.0025).and.(ne.le.neini))then     
!     
!     Rebound (or too high variation during pers. contact)
!     r_hat_bc is r_hat before contact
!     
            idivergence=1
            sizemaxinc=dtheta*0.5
            iforceincsize=1
c            icase=3
!     
         else
!     
!     Persistent Contact
!     
c            icase=2
            if(r_hat.gt.(-0.9*maxdecay))then
               sizemaxinc=min(fact*temax/tper,1.1*dtheta)
               sizemaxinc=min(sizemaxinc,100.0*temax/tper)
            else
               sizemaxinc=max(temax/tper/10.0,0.5*dtheta)
               istab=1
            endif
         endif
      endif                     !(mortar)
!     
      if(sizemaxinc.lt.tmin)then
         sizemaxinc=tmin
      endif
!     
!     Debug prints
!     
c      if(icase.eq.1)then
c         write(*,*)"# # # # # # # # # # # # # # # # # # # # # # # # # #"
c         write(*,*)"icase=1"
c         write(*,*)"Impact detected, increment reattempted"
c      elseif(icase.eq.2)then
c         write(*,*)"# # # # # # # # # # # # # # # # # # # # # # # # # #"
c         write(*,*)"icase=2"
c         write(*,*)"Persistent Contact conditions"
c      elseif(icase.eq.3)then
c         write(*,*)"# # # # # # # # # # # # # # # # # # # # # # # # # #"
c         write(*,*)"icase=3"
c         write(*,*)"Rebound, the increment is reattempted (stability)"
c      endif
!     
      return
      end
