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
      subroutine resultsforc_se(nk,f,dfn,nactdof,ipompc,nodempc,
     &  coefmpc,labmpc,nmpc,mi,fmpc,calcul_fn,calcul_f,ndesi,df)
!
!     calculating the equation system internal force vector
!     (one entry for each active degree of freedom)
!
      implicit none
!
      character*20 labmpc(*)
!
      integer mi(*),nactdof(0:mi(2),*),ipompc(*),nodempc(3,*),nk,i,j,
     &  nmpc,ist,ndir,node,index,calcul_fn,calcul_f,ndesi,desvar
!
      real*8 f(*),dfn(ndesi,0:mi(2),*),coefmpc(*),fmpc(*),forcempc,
     &  df(ndesi,*)
!
!     subtracting the mpc force (for each linear mpc there is one
!     force; the actual force in a node belonging to the mpc is
!     obtained by multiplying this force with the nodal coefficient.
!     The force has to be subtracted from f, since it does not
!     appear on the rhs of the equation system
!
!     Loop over all designvariables
      do desvar=1,ndesi
!
      if(calcul_fn.eq.1)then
        do i=1,nmpc
            ist=ipompc(i)
            node=nodempc(1,ist)
            ndir=nodempc(2,ist)
            if(ndir.gt.3) cycle
            forcempc=dfn(desvar,ndir,node)/coefmpc(ist)
            fmpc(i)=forcempc
            dfn(desvar,ndir,node)=0.d0
            index=nodempc(3,ist)
            if(index.eq.0) cycle
            do
               node=nodempc(1,index)
               ndir=nodempc(2,index)
               dfn(desvar,ndir,node)=dfn(desvar,ndir,node)-
     &coefmpc(index)*forcempc
               index=nodempc(3,index)
               if(index.eq.0) exit
            enddo
         enddo
      endif
!
!     calculating the system force vector
!
      if(calcul_f.eq.1) then
         do i=1,nk
            do j=0,mi(2)
               if(nactdof(j,i).ne.0) then
                  df(desvar,nactdof(j,i))=dfn(desvar,j,i)
               endif
            enddo
         enddo
      endif
!
!     adding the mpc force again to fn
!
c      if(calcul_fn.eq.1)then
c         do i=1,nmpc
c            ist=ipompc(i)
c            node=nodempc(1,ist)
c            ndir=nodempc(2,ist)
c            if(ndir.gt.3) cycle
c            forcempc=fmpc(i)
c            dfn(desvar,ndir,node)=forcempc*coefmpc(ist)
c            index=nodempc(3,ist)
!
!           nodes not belonging to the structure have to be
!           taken out
!
c            if(labmpc(i)(1:7).eq.'MEANROT') then
c               if(nodempc(3,nodempc(3,index)).eq.0) cycle
c            elseif(labmpc(i)(1:10).eq.'PRETENSION') then
c               if(nodempc(3,index).eq.0) cycle
c            elseif(labmpc(i)(1:5).eq.'RIGID') then
c               if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,inde
c     &x))))).eq.0) cycle
c            else
c               if(index.eq.0) cycle
c            endif
c            do
c               node=nodempc(1,index)
c               ndir=nodempc(2,index)
c               dfn(desvar,ndir,node)=dfn(desvar,ndir,node)+
c     &coefmpc(index)*forcempc
c               index=nodempc(3,index)
c               if(labmpc(i)(1:7).eq.'MEANROT') then
c                  if(nodempc(3,nodempc(3,index)).eq.0) exit
c               elseif(labmpc(i)(1:10).eq.'PRETENSION') then
c                  if(nodempc(3,index).eq.0) exit
c               elseif(labmpc(i)(1:5).eq.'RIGID') then
c                  if(nodempc(3,nodempc(3,nodempc(3,nodempc(3,nodempc(3,i
c     &ndex))))).eq.0) exit
c               else
c                  if(index.eq.0) exit
c               endif
c            enddo
c         enddo
c      endif
!
!     end of loop over all designvariables
      enddo
!
      return
      end
