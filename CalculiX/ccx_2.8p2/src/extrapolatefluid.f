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
      subroutine extrapolatefluid(nk,iponofa,inofa,inum,vfa,v,ielfa,
     &  ithermal)
!
!     extrapolates the field values at the center of the faces to
!     the nodes
!
      implicit none
!
      integer nk,iponofa(*),inofa(2,*),inum(*),ielfa(4,*),i,l,indexf,
     &  iface,ithermal
!
      real*8 vfa(0:5,*),v(0:4,*)
!
      do i=1,nk
         if(ithermal.eq.0) then
            do l=1,4
               v(l,i)=0.d0
            enddo
            inum(i)=0
            indexf=iponofa(i)
            do
               if(indexf.eq.0) exit
               iface=inofa(1,indexf)
               do l=1,4
                  v(l,i)=v(l,i)+vfa(l,iface)
c                  if(i.eq.23324) then
c                     if(l.eq.3) then
c                     write(*,*) 'extrapolatefluid ',i,iface,vfa(l,iface)
c                     endif
c                  endif
               enddo
               inum(i)=inum(i)+1
               indexf=inofa(2,indexf)
            enddo
            if(inum(i).gt.0) then
               do l=1,4
                  v(l,i)=v(l,i)/inum(i)
               enddo
            endif
         else
            do l=0,4
               v(l,i)=0.d0
            enddo
            inum(i)=0
            indexf=iponofa(i)
            do
               if(indexf.eq.0) exit
               iface=inofa(1,indexf)
               do l=0,4
                  v(l,i)=v(l,i)+vfa(l,iface)
               enddo
               inum(i)=inum(i)+1
               indexf=inofa(2,indexf)
            enddo
            if(inum(i).gt.0) then
               do l=0,4
                  v(l,i)=v(l,i)/inum(i)
               enddo
            endif
         endif
      enddo
!  
      return
      end
