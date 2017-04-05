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
      subroutine extrapolate_vel(nface,ielfa,xrlfa,vel,vfa,
     &  ifabou,xboun,ipnei,nef,icyclic,c,ifatie,xxn)
!
!     inter/extrapolation of v at the center of the elements
!     to the center of the faces
!
      implicit none
!
      integer nface,ielfa(4,*),ifabou(*),iel1,iel2,iel3,i,j,ipointer,
     &  indexf,ipnei(*),nef,icyclic,ifatie(*)
!
      real*8 xrlfa(3,*),vel(nef,0:7),vfa(0:7,*),xboun(*),xl1,xl2,
     &  c(3,3),xxn(3,*),dd
!
c$omp parallel default(none)
c$omp& shared(nface,ielfa,xrlfa,vfa,vel,ipnei,ifabou,xboun,
c$omp&        icyclic,c,ifatie,xxn)
c$omp& private(i,iel1,xl1,iel2,xl2,j,iel3,ipointer,indexf,dd)
c$omp do
      do i=1,nface
         iel1=ielfa(1,i)
         xl1=xrlfa(1,i)
         iel2=ielfa(2,i)
         if(iel2.gt.0) then
            xl2=xrlfa(2,i)
            if((icyclic.eq.0).or.(ifatie(i).eq.0)) then
               do j=1,3
                  vfa(j,i)=xl1*vel(iel1,j)+xl2*vel(iel2,j)
               enddo
            elseif(ifatie(i).gt.0) then
               do j=1,3
                  vfa(j,i)=xl1*vel(iel1,j)+xl2*
     &              (c(j,1)*vel(iel2,1)
     &              +c(j,2)*vel(iel2,2)
     &              +c(j,3)*vel(iel2,3))
               enddo
            else
               do j=1,3
                  vfa(j,i)=xl1*vel(iel1,j)+xl2*
     &              (c(1,j)*vel(iel2,1)
     &              +c(2,j)*vel(iel2,2)
     &              +c(3,j)*vel(iel2,3))
               enddo
            endif
         else
            indexf=ipnei(iel1)+ielfa(4,i)
            iel3=ielfa(3,i)
c            write(*,*) 'extrapolate_vel ',iel1,iel2,iel3
            if(iel2.lt.0) then
               ipointer=-iel2
!
!              global x-direction
!
               if(ifabou(ipointer+1).gt.0) then
!     
!     v_1 given
!     
                  vfa(1,i)=xboun(ifabou(ipointer+1))
               elseif((ifabou(ipointer+4).ne.0).and.(iel3.ne.0)) then
!     
!     p given; linear interpolation
!     
                  vfa(1,i)=xl1*vel(iel1,1)+xrlfa(3,i)*vel(iel3,1)
               else
!     
!     constant extrapolation
!     
                  vfa(1,i)=vel(iel1,1)
               endif
!     
!     global y-direction
!     
               if(ifabou(ipointer+2).gt.0) then
!     
!     v_2 given
!     
                  vfa(2,i)=xboun(ifabou(ipointer+2))
               elseif((ifabou(ipointer+4).ne.0).and.(iel3.ne.0)) then
!     
!     p given; linear interpolation
!     
                  vfa(2,i)=xl1*vel(iel1,2)+xrlfa(3,i)*vel(iel3,2)
               else
!     
!     constant extrapolation
!     
                  vfa(2,i)=vel(iel1,2)
               endif
!     
!     global z-direction
!     
               if(ifabou(ipointer+3).gt.0) then
!     
!     v_3 given
!     
                  vfa(3,i)=xboun(ifabou(ipointer+3))
               elseif((ifabou(ipointer+4).ne.0).and.(iel3.ne.0)) then
!     
!     p given; linear interpolation
!     
                  vfa(3,i)=xl1*vel(iel1,3)+xrlfa(3,i)*vel(iel3,3)
               else
!     
!     constant extrapolation
!     
                  vfa(3,i)=vel(iel1,3)
               endif
!
!     correction for sliding boundary conditions        
!
c               if(ifabou(ipointer+5).eq.2) then
               if(ifabou(ipointer+5).lt.0) then
                  dd=vfa(1,i)*xxn(1,indexf)+
     &               vfa(2,i)*xxn(2,indexf)+
     &               vfa(3,i)*xxn(3,indexf)
                  do j=1,3
                     vfa(j,i)=vfa(j,i)-dd*xxn(j,indexf)
                  enddo
               endif
!     
            else
!     
!     constant extrapolation
!     
               do j=1,3
                  vfa(j,i)=vel(iel1,j)
               enddo
            endif
         endif
      enddo
c$omp end do
c$omp end parallel
!
c      write(*,*) 'extrapolate_vel '
c      do i=1,nef
c         write(*,*) i,(vel(i,j),j=0,5)
c      enddo
c      do i=1,nface
c         write(*,*) i,(vfa(j,i),j=0,5)
c      enddo
      return
      end
