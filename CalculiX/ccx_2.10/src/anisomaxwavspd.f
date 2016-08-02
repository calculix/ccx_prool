!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2007 Guido Dhondt
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
      subroutine anisomaxwavspd(elas,rho,iorth,maxwavspd)
!
!     Calculates the propagation wave speed in a material, up to its 21 
!     constants. Subroutine for calcmatwavsps.f

!     Based on the procedure in:
!     C. Lane. The Development of a 2D Ultrasonic Array Inspection 
!     for Single Crystal Turbine Blades.
!     Switzerland: Springer International Publishing, 2014.
!
!      CARLO MONJARAZ TEC (CMT)
!
!       INPUT:
!       
!       elas: double(21) - The elasticity vector, containing 21 entries. 
!             Non used are zero. If material is iorthtropic, values are
!             rearranged in middle step to match indexes from anisotropic
!             material card.
!             
!        rho: double - Density of the material
!        
!       iorth: INTEGER - if the value is 1 : material is iorthtropic
!                        for other vaules:   material is anisotropic
!                        
!       maxwavspd: double(*) - contains the list of maximum wave speeds
!                  per each material defined in the model. Only materials
!                  that are used in the model are updated, else are left
!                  as -1.0                 
!        
!       OUTPUT:
!       
!       maxwavspd
!
      implicit none
!
      real*8 elas(21),c(3,3,3,3),rho,n(3),cm(3,3,3),
     &       cmm(3,3),dd,al(3),alz(3,3),fv1(3),fv2(3),
     &       theta(2*100+1),phi(100+1),pi,p1(3),p2(3),p3(3),v(3),
     &       maxwavspd,speed,nmax(3)
!     
      integer i,j,k,l,it,jt,kt,lt,imax,jmax,itheta,jphi,
     &        itensor(4,21),ier,matz,ndim,iorth
!
      data itensor /1,1,1,1,
     &              1,1,2,2,
     &              2,2,2,2,
     &              1,1,3,3,
     &              2,2,3,3,
     &              3,3,3,3,
     &              1,1,1,2,
     &              2,2,1,2,
     &              3,3,1,2,
     &              1,2,1,2,
     &              1,1,1,3,
     &              2,2,1,3,
     &              3,3,1,3,
     &              1,2,1,3,
     &              1,3,1,3,
     &              1,1,2,3,
     &              2,2,2,3,
     &              3,3,2,3,
     &              1,2,2,3,
     &              1,3,2,3,
     &              2,3,2,3/      
!     
      pi=4.d0*datan(1.d0)
!     
      speed=-1.d0
      maxwavspd=-1.d0
!     
      write(*,*)'++cMT: calculating max. speed in ANISOTROPIc...'
!     
!--------IF IORTHTROPIc-----------------------------     
      if(iorth.eq.1)then
!     
         elas(10)=elas(9)
         elas(15)=elas(8)
         elas(21)=elas(7)
         elas(9)=0.d0
         elas(8)=0.d0
         elas(7)=0.d0
!     
      endif
!     
      maxwavspd=-1.d0
      do i=1,101
         phi(i)=(i-1.d0)/100.d0*pi
      enddo
!     
      do i=1,201
         theta(i)=((i-1.d0)/100.d0-1.d0 )*pi
      enddo
!--------FIlling  c voigt Matrix-----------------------------            
      do i=1,21
         it=itensor(1,i)
         jt=itensor(2,i)
         kt=itensor(3,i)
         lt=itensor(4,i)          
!     
         c(lt,kt,jt,it)=elas(i)
         c(lt,kt,it,jt)=elas(i)
         c(kt,lt,jt,it)=elas(i)
         c(jt,it,lt,kt)=elas(i)
!     2/3, 1/2, 1/3 swap
         it=itensor(2,i)
         jt=itensor(1,i)
         kt=itensor(4,i)
         lt=itensor(3,i)
         c(lt,kt,jt,it)=elas(i)
         c(lt,kt,it,jt)=elas(i)
         c(kt,lt,jt,it)=elas(i)
         c(jt,it,lt,kt)=elas(i)
      enddo
!     
!-----for each direction in the unit sphere------------------
!     
      itheta=170
      jphi=30
      do itheta=1,201
         do jphi=1,101
!     
            do l=1,3
               n(l)=0.
               v(l)=0.
               do k=1,3
                  cmm(k,l)=0.
                  do j=1,3
                     cm(j,l,k)=0.
                  enddo
               enddo
            enddo            
!     
            n(1)=cos(theta(itheta))*Sin(phi(jphi))
            n(2)=Sin(theta(itheta))*Sin(phi(jphi))
            n(3)=cos(phi(jphi))
!     
!     c ------------ PER EAcH DIREcTION find wave speed-----------------------
!     
            dd=dsqrt(n(1)*n(1)+n(2)*n(2) +n(3)*n(3))
            
            n(1)=n(1) / dd
            n(2)=n(2) / dd
            n(3)=n(3) / dd
!     
            do l=1,3
               do k=1,3
                  do i=1,3
                     do j=1,3
                        cm(l,k,i)=cm(l,k,i)+c(l,k,j,i)*n(j)
                     enddo
                  enddo        
               enddo
            enddo
!     
            do k=1,3
               do i=1,3
                  do l=1,3
                     cmm(k,i)=cmm(k,i)+cm(l,k,i)*n(l)
                  enddo
               enddo        
            enddo
!     
            ndim=3
            matz=1
            ier=0
!
!     ---------reset vars for EIGvALUES
!
            do j=1,3
               al(j)=0.
               fv1(j)=0.
               fv2(j)=0.
               do i=1,3
                  alz(j,i)=0.
               enddo
            enddo
!     
            call rs(ndim,ndim,cmm,al,matz,alz,fv1,fv2,ier)
!     
!           ------normalizing eigenvectors to P vectors----------
!     
            dd=dsqrt(alz(1,1)**2+alz(2,1)**2+alz(3,1)**2)
            p1(1)=alz(1,1) / dd
            p1(2)=alz(2,1) / dd
            p1(3)=alz(3,1) / dd
            dd=dsqrt(alz(1,2)**2+alz(2,2)**2+alz(3,2)**2)
            p2(1)=alz(1,2) / dd
            p2(2)=alz(2,2) / dd
            p2(3)=alz(3,2) / dd
            dd=dsqrt(alz(1,3)**2+alz(2,3)**2+alz(3,3)**2)
            p3(1)=alz(1,3) / dd
            p3(2)=alz(2,3) / dd
            p3(3)=alz(3,3) / dd
!     
            do l=1,3
               do k=1,3
                  cmm(k,l)=0.
                  do j=1,3
                     cm(j,l,k)=0.
                  enddo
               enddo
            enddo
!     
            do l=1,3
               do j=1,3
                  do i=1,3
                     do k=1,3
                        cm(l,j,i)=cm(l,j,i)+c(l,k,j,i)*n(k);
                     enddo
                  enddo        
               enddo
            enddo
!     
            do  l=1,3
               do  j=1,3        
                  do  i=1,3
                     cmm(l,j)=cmm(l,j)+cm(l,j,i)*p3(i);        
                  enddo
               enddo        
            enddo
!     
            do j=1,3
               do i=1,3
                  v(j)=v(j)+cmm(j,i)*p3(i)
               enddo
            enddo
!     
            dd=dsqrt(v(1)**2+v(2)**2+v(3)**2) 
            speed=dsqrt(dd/rho) 
!     
            if(speed.gt.maxwavspd)then
               maxwavspd=speed 
               imax=itheta
               jmax=jphi
               nmax=n
            endif
            
         enddo
      enddo        
!-------END for each direction in the unit sphere-------------
!     
      return
      end
      
