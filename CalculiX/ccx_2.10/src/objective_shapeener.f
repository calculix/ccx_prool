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
      subroutine objective_shapeener(co,kon,ipkon,lakon,ne,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,iperturb,fn,iout,qa,vold,nmethod,
     &  veold,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
     &  xstateini,xstiff,xstate,npmat_,matname,mi,ielas,icmd,
     &  ncmat_,nstate_,stiini,vini,ener,enerini,istep,iinc,
     &  springarea,reltime,calcul_qa,iener,ikin,nal,ne0,thicke,
     &  emeini,pslavsurf,pmastsurf,mortar,clearini,nea,neb,
     &  ielprop,prop,distmin,ndesi,nodedesi,ndirdesi,nobject,
     &  g0,dgdx,numobject,sti,dgdxtot,dfextminds,df,neq,nactdof,
     &  nk)
!
!     calculates the total differential of the objective function 
!     shape energy
!
      implicit none
!
      character*8 lakon(*)
      character*80 matname(*)
!
      integer kon(*),nea,neb,mi(*),nelcon(2,*),nrhcon(*),
     &  nalcon(2,*),ielmat(mi(3),*),ielorien(mi(3),*),
     &  ntmat_,ipkon(*),ne0,istep,iinc,ne,ithermal(2),
     &  iprestr,iener,norien,iperturb(*),iout,nal,icmd,nmethod,
     &  ielas,ncmat_,nstate_,ikin,ielprop(*),nplicon(0:ntmat_,*),
     &  nplkcon(0:ntmat_,*),npmat_,calcul_qa,mortar,ndesi,
     &  nodedesi(*),ndirdesi(*),nobject,numobject,neq,
     &  nactdof(0:mi(2),*),nk
!
      real*8 co(3,*),stiini(6,mi(1),*),stx(6,mi(1),*),
     &  prop(*),elcon(0:ncmat_,ntmat_,*),rhcon(0:1,ntmat_,*),
     &  alcon(0:6,ntmat_,*),vini(0:mi(2),*),alzero(*),orab(7,*),
     &  fn(0:mi(2),*),t0(*),t1(*),prestr(6,mi(1),*),
     &  vold(0:mi(2),*),veold(0:mi(2),*),springarea(2,*),
     &  ener(mi(1),*),enerini(mi(1),*),qa(3),dtime,
     &  time,ttime,plicon(0:2*npmat_,ntmat_,*),
     &  plkcon(0:2*npmat_,ntmat_,*),xstiff(27,mi(1),*),
     &  xstate(nstate_,mi(1),*),xstateini(nstate_,mi(1),*),reltime,
     &  thicke(mi(3),*),emeini(6,mi(1),*),clearini(3,9,*),
     &  pslavsurf(3,*),pmastsurf(6,*),distmin,g0(nobject),
     &  dgdx(ndesi,nobject),dgdxtot(ndesi,nobject),sti(6,mi(1),*), 
     &  dfextminds(ndesi,neq),df(ndesi,neq)
!
!     -------------------------------------------------------------
!     Calculation of the derivative w.r.t to the designvaribales
!     (For the objective function mass this is equal to the
!     total difference
!     -------------------------------------------------------------
!
      call objective_shapeener_dx(co,kon,ipkon,lakon,ne,
     &  stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
     &  ielmat,ielorien,norien,orab,ntmat_,t0,t1,ithermal,prestr,
     &  iprestr,iperturb,fn,iout,qa,vold,nmethod,veold,dtime,
     &  time,ttime,plicon,nplicon,plkcon,nplkcon,xstateini,xstiff,
     &  xstate,npmat_,matname,mi,ielas,icmd,ncmat_,nstate_,
     &  stiini,vini,ener,enerini,istep,iinc,springarea,reltime,
     &  calcul_qa,iener,ikin,nal,ne0,thicke,emeini,pslavsurf,
     &  pmastsurf,mortar,clearini,nea,neb,ielprop,prop,distmin,
     &  ndesi,nodedesi,ndirdesi,nobject,g0,dgdx,numobject,sti)
!
!     -------------------------------------------------------------
!     Calculation of the total derivative of the shape energy w.r.t.
!     design variables
!     -------------------------------------------------------------
!
      call objective_shapeener_tot(dgdxtot,dgdx,dfextminds,df,vold,
     &  neq,ndesi,nobject,numobject,mi,nactdof,nk)
      
      return
      end
