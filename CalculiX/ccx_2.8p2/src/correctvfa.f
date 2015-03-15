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
      subroutine correctvfa(nface,ielfa,area,vfa,ap,bp,xxn,
     &  ifabou,ipnei,nef,neifa,hfa,vel,xboun)
!
!     correction of v due to the balance of mass
!     the correction is in normal direction to the face
!
      implicit none
!
      integer i,nface,ielfa(4,*),iel1,iel2,ifabou(*),j,indexf,k,
     &  ipnei(*),ifa,nef,neifa(*),indexb
!
      real*8 ap(*),bp(*),area(*),vfa(0:5,*),xxn(3,*),
     &  flux,hfa(3,*),vel(nef,0:5),dh,xboun(*)
!
      do i=1,nface
!
!        first neighboring element
!
         iel1=ielfa(1,i)
         j=ielfa(4,i)
         indexf=ipnei(iel1)+j
!
!        second neighboring element
!
         iel2=ielfa(2,i)
!
!        factor between mass flow and velocity
!
         ap(i)=ap(i)/(area(i)*vfa(5,i))
!
!        internal face
!
         if(iel2.gt.0) then
            dh=hfa(1,i)*xxn(1,indexf)+
     &         hfa(2,i)*xxn(2,indexf)+
     &         hfa(3,i)*xxn(3,indexf)
            do k=1,3
               vfa(k,i)=(dh-ap(i)*(vel(iel2,4)-vel(iel1,4)+bp(i)))
     &                  *xxn(k,indexf)
            enddo
!
!        external face with pressure boundary condition
!
         elseif(iel2.lt.0) then
            dh=hfa(1,i)*xxn(1,indexf)+
     &         hfa(2,i)*xxn(2,indexf)+
     &         hfa(3,i)*xxn(3,indexf)
            indexb=ifabou(-iel2+4)
            if(indexb.ne.0) then
               do k=1,3
                  vfa(k,i)=(dh-ap(i)*(xboun(indexb)-vel(iel1,4)+bp(i)))
     &                          *xxn(k,indexf)
               enddo
            endif
         endif
c         write(*,*) 'correctvfa ',iel1,j,vfa(1,i),vfa(2,i),vfa(3,i)
c         if(iel2.lt.0) write(*,*) 'correctvfa1 ',iel1,j,
c     &      vfa(1,i)*xxn(1,indexf)+vfa(2,i)*xxn(2,indexf)+
c     &      vfa(3,i)*xxn(3,indexf)
      enddo
!
!     check conservation of mass (not satisfied if alpha != 1)
!
      do i=1,nef
         flux=0.d0
         indexf=ipnei(i)
         do j=1,6
            indexf=indexf+1
            ifa=neifa(indexf)
            flux=flux+area(ifa)*(vfa(1,ifa)*xxn(1,indexf)+
     &                vfa(2,ifa)*xxn(2,indexf)+
     &                vfa(3,ifa)*xxn(3,indexf))
         enddo
c         write(*,*) 'correctvfa mass check ',i,flux
      enddo
!  
      return
      end
