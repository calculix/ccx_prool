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
      subroutine objective_shapeener_tot(dgdxtot,dgdx,dfextminds,
     &  df,vold,neq,ndesi,nobject,numobject,mi,nactdof,nk)
!
      integer neq,ndesi,nobject,numobject,mi(*),i,j,jj,
     &  nactdof(0:mi(2),*)
!      
      real*8 dgdxtot(ndesi,nobject),dgdx(ndesi,nobject),
     &  dfextminds(ndesi,neq),df(ndesi,neq),vold(0:mi(2),*),
     &  dummy1(ndesi,neq),dummy2(ndesi)
!
!     ----------------------------------------------------------------
!     Calculation of the total differential:
!     dgdxtot = dgdx + vold^(T) * ( dfextminds - df )
!     ----------------------------------------------------------------
!     
      do i=1,ndesi
         do j=1,neq
            dummy1(i,j)=dfextminds(i,j)-df(i,j)
         enddo
      enddo
!     
      do jj=1,ndesi
         dummy2(jj)=0.d0
         do i=1,nk            
            do j=1,3
               if(nactdof(j,i).eq.0) cycle
          dummy2(jj)=dummy2(jj)+
     &             vold(j,i)*dummy1(jj,nactdof(j,i))
            enddo         
         enddo
      enddo
!          
      do i=1,ndesi
         dgdxtot(i,numobject)=dgdx(i,numobject)+dummy2(i)
      enddo
!      
      return
      end
