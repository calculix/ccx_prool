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
!     construction of the B matrix      
!
!     author: Yannick Muller
!     
      subroutine flowoutput(itg,ieg,ntg,nteq,bc,lakon,ntmat_,
     &     v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
     &     ielmat,prop,ielprop,nactdog,nacteq,iin,physcon,
     &     camt,camf,camp,uamt,uamf,uamp,rhcon,nrhcon,
     &     vold,jobnamef,set,istartset,iendset,ialset,nset,mi)
!     
      implicit none
!     
      character*8 lakon(*)
      character*132 jobnamef(*)
      character*81 set(*)
!     
      integer mi(*),itg(*),ieg(*),ntg,nflow,ielmat(mi(3),*),i,
     &     nrhcon(*),node,
     &     ntmat_,nteq,nshcon(*),nelem,index,ipkon(*),kon(*),iin,
     &     nactdog(0:3,*),nacteq(0:3,*),ielprop(*),
     &     istartset(*),iendset(*),ialset(*),nset
!     
      real*8 physcon(3),v(0:mi(2),*),shcon(0:3,ntmat_,*),co(3,*),
     &     prop(*),dtime,ttime,time,xflow,camt(*),camf(*),camp(*),
     &     rhcon(0:1,ntmat_,*),vold(0:mi(2),*),uamt,uamf,uamp,eta,
     &     bc(*)
!
      do i=1,nflow
         if(lakon(ieg(i))(2:5).eq.'LICH') then
            if((lakon(ieg(i))(6:7).eq.'SG').or.
     &           (lakon(ieg(i))(6:7).eq.'WE').or.
     &           (lakon(ieg(i))(6:7).eq.'DS')) then
               index=ipkon(ieg(i))
               node=kon(index+2)
               if(nactdog(3,node).eq.0) cycle
               index=ielprop(ieg(i))
               if(lakon(ieg(i))(6:7).eq.'SG') then
                  eta=prop(index+4)
                  nelem=int(prop(index+7))      
               elseif(lakon(ieg(i))(6:7).eq.'WE') then
                  eta=prop(index+4)
                  nelem=int(prop(index+7))      
               elseif(lakon(ieg(i))(6:7).eq.'DS') then
                  eta=prop(index+7)
                  nelem=int(prop(index+9))      
               endif
               if(nelem.ne.0) then
                  write(*,*) '     *INFO in flowoutput: hydraulic jump'
                  write(*,*) '           in element ',nelem,'.'
                  write(*,*) '           relative location:',eta
                  write(*,*)
               endif
            endif
         endif
      enddo
!
      return
      end
