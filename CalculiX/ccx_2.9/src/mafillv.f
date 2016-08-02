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
      subroutine mafillv(nef,ipnei,neifa,neiel,vfa,xxn,area,
     &  auv,adv,jq,irow,nzs,bv,vel,cosa,umfa,xlet,xle,gradvfa,xxi,
     &  body,volume,ielfa,lakon,ifabou,nbody,neq,
     &  dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,a1,
     &  a2,a3,flux)
!
      implicit none
!
      character*2 one,two,three
      character*8 lakon(*)
!
      integer i,nef,jdof1,indexf,ipnei(*),j,ifa,iel,neifa(*),
     &  neiel(*),jdof2,jq(*),irow(*),nzs,iwall,ielfa(4,*),
     &  ipointer,ifabou(*),nbody,neq,k,indexb,numfaces,nactdohinv(*)
!
      real*8 xflux,vfa(0:5,*),xxn(3,*),area(*),auv(*),adv(*),bv(neq,3),
     &  vel(nef,0:5),cosa(*),umfa(*),xlet(*),xle(*),coef,gradvfa(3,3,*),
     &  xxi(3,*),body(0:3,*),volume(*),coef2,dtimef,velo(nef,0:5),
     &  veloo(nef,0:5),rhovel,constant,sel(3,*),xrlfa(3,*),gamma(*),
     &  xxj(3,*),a1,a2,a3,flux(*)
!
      one='b1'
      two='b2'
      three='b3'
!
      do i=1,nef
         jdof1=i
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakon(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
!
!              convection
!
            indexf=indexf+1
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) jdof2=iel
            xflux=flux(indexf)
c            xflux=(vfa(1,ifa)*xxn(1,indexf)+
c     &           vfa(2,ifa)*xxn(2,indexf)+
c     &           vfa(3,ifa)*xxn(3,indexf))
c     &           *vfa(5,ifa)*area(ifa)
!
!              vfa sometimes undefined?
!
            if(xflux.ge.0.d0) then
!
!                 outflowing xflux
!
               call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,
     &              xflux,nzs)
               bv(jdof1,1)=bv(jdof1,1)-(vfa(1,ifa)-vel(i,1))*xflux
               bv(jdof1,2)=bv(jdof1,2)-(vfa(2,ifa)-vel(i,2))*xflux
               bv(jdof1,3)=bv(jdof1,3)-(vfa(3,ifa)-vel(i,3))*xflux
c               write(*,*) 'mafillv1 ',bv(1,1),bv(1,2),bv(1,3)
            else
               if(iel.gt.0) then
!
!                    incoming flux from neighboring element
!
                  call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof2,xflux,
     &                 nzs)
                  bv(jdof1,1)=bv(jdof1,1)-(vfa(1,ifa)-vel(iel,1))*xflux
                  bv(jdof1,2)=bv(jdof1,2)-(vfa(2,ifa)-vel(iel,2))*xflux
                  bv(jdof1,3)=bv(jdof1,3)-(vfa(3,ifa)-vel(iel,3))*xflux
c               write(*,*) 'mafillv2 ',bv(1,1),bv(1,2),bv(1,3)
               else
!
!                    incoming flux through boundary
!
                  if(ielfa(2,ifa).lt.0) then
                     indexb=-ielfa(2,ifa)
                     if(((ifabou(indexb+1).ne.0).and.
     &                    (ifabou(indexb+2).ne.0).and.
     &                    (ifabou(indexb+3).ne.0)).or.
c     &                    (dabs(xflux).lt.1.d-3)) then
     &                    (dabs(xflux).lt.1.d-10)) then
                        bv(jdof1,1)=bv(jdof1,1)-vfa(1,ifa)*xflux
                        bv(jdof1,2)=bv(jdof1,2)-vfa(2,ifa)*xflux
                        bv(jdof1,3)=bv(jdof1,3)-vfa(3,ifa)*xflux
c               write(*,*) 'mafillv3 ',bv(1,1),bv(1,2),bv(1,3)
                     else
                        write(*,*) '*ERROR in mafillv: not all'
                        write(*,*) '       components of an incoming'
                        write(*,*) '       flux through face ',j
                        write(*,*)'       of element ',nactdohinv(i),
     &                        ' are given'
                     endif
                  else
                     write(*,*) '*ERROR in mafillv: not all'
                     write(*,*) '       components of an incoming'
                     write(*,*) '       flux through face ',j
                     write(*,*)'       of element ',nactdohinv(i),
     &                     ' are given'
                  endif
               endif
            endif
!
!              diffusion
!
            if(iel.ne.0) then
!
!                 neighboring element
!
               coef=umfa(ifa)*area(ifa)/xlet(indexf)
               call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,
     &              coef,nzs)
               call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof2,
     &              -coef,nzs)
!
!                 correction for non-orthogonal grid
!
               bv(jdof1,1)=bv(jdof1,1)+umfa(ifa)*area(ifa)*
     &              (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradvfa(1,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradvfa(1,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
               bv(jdof1,2)=bv(jdof1,2)+umfa(ifa)*area(ifa)*
     &              (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradvfa(2,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradvfa(2,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
               bv(jdof1,3)=bv(jdof1,3)+umfa(ifa)*area(ifa)*
     &              (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradvfa(3,2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradvfa(3,3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
c               write(*,*) 'mafillv4 ',bv(1,1),bv(1,2),bv(1,3)
            else
!
!                 boundary; check whether wall (specified by user),
!                           outlet (no velocity boundary conditions or
!                           none of those
!
               iwall=0
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  iwall=ifabou(ipointer+5)
               endif
               if(iwall.eq.0) then
!
!                    external face, but no wall
!
                  if((ifabou(ipointer+1).ne.0).or.
     &                 (ifabou(ipointer+2).ne.0).or.
     &                 (ifabou(ipointer+3).ne.0)) then
!
!                       no outlet: face velocity fixed
!
                     coef=umfa(ifa)*area(ifa)/xle(indexf)
                     call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,
     &                    coef,nzs)
                     bv(jdof1,1)=bv(jdof1,1)+coef*vfa(1,ifa)
                     bv(jdof1,2)=bv(jdof1,2)+coef*vfa(2,ifa)
                     bv(jdof1,3)=bv(jdof1,3)+coef*vfa(3,ifa)
c               write(*,*) 'mafillv5 ',bv(1,1),bv(1,2),bv(1,3)
                  else
!
!                       outlet: no diffusion
!
                  endif
!
!                    correction for non-orthogonal grid
!
                  bv(jdof1,1)=bv(jdof1,1)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(1,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradvfa(1,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradvfa(1,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  bv(jdof1,2)=bv(jdof1,2)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(2,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradvfa(2,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradvfa(2,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
                  bv(jdof1,3)=bv(jdof1,3)+umfa(ifa)*area(ifa)*
     &                 (gradvfa(3,1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradvfa(3,2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradvfa(3,3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
c               write(*,*) 'mafillv6 ',bv(1,1),bv(1,2),bv(1,3)
               else
!     
!                    wall
!     
                  coef=umfa(ifa)*area(ifa)/(xle(indexf)*cosa(indexf))
                  call add_sm_fl_as(auv,adv,jq,irow,jdof1,jdof1,coef,
     &                 nzs)
!
!                    correction for non-orthogonal grid and nonzero
!                    wall velocity
!
                  coef2=((vel(i,1)-vfa(1,ifa))*xxn(1,indexf)+
     &                 (vel(i,2)-vfa(2,ifa))*xxn(2,indexf)+
     &                 (vel(i,3)-vfa(3,ifa))*xxn(3,indexf))*coef
                  bv(jdof1,1)=bv(jdof1,1)+coef*vfa(1,ifa)+
     &                 coef2*xxn(1,indexf)
                  bv(jdof1,2)=bv(jdof1,2)+coef*vfa(2,ifa)+
     &                 coef2*xxn(2,indexf)
                  bv(jdof1,3)=bv(jdof1,3)+coef*vfa(3,ifa)+
     &                 coef2*xxn(3,indexf)
c               write(*,*) 'mafillv7 ',bv(1,1),bv(1,2),bv(1,3)
               endif
            endif
!     
!     pressure
!     
         enddo
!     
!           body force
!     
         rhovel=vel(i,5)*volume(i)
!
         if(nbody.gt.0) then
            bv(jdof1,1)=bv(jdof1,1)+rhovel*body(1,i)
            bv(jdof1,2)=bv(jdof1,2)+rhovel*body(2,i)
            bv(jdof1,3)=bv(jdof1,3)+rhovel*body(3,i)
c               write(*,*) 'mafillv9 ',bv(1,1),bv(1,2),bv(1,3)
         endif
!
!           transient term
!
         constant=rhovel
         bv(jdof1,1)=bv(jdof1,1)-(a2*velo(i,1)+a3*veloo(i,1))*constant
         bv(jdof1,2)=bv(jdof1,2)-(a2*velo(i,2)+a3*veloo(i,2))*constant
         bv(jdof1,3)=bv(jdof1,3)-(a2*velo(i,3)+a3*veloo(i,3))*constant
c               write(*,*) 'mafillv10 ',bv(1,1),bv(1,2),bv(1,3)
c               write(*,*) 'mafillv10 ',velo(1,1),veloo(1,1),constant
c               write(*,*) 'mafillv10 ',rhovel,dtimef,vel(i,5)
         constant=a1*constant
         call add_sm_fl(auv,adv,jq,irow,jdof1,jdof1,constant,nzs)
!
!           copying b into sel (rhs without pressure)
!
         do j=1,3
            sel(j,jdof1)=bv(jdof1,j)
         enddo
c         write(*,*) 'mafillv sel ',sel(1,1),sel(2,1),sel(3,1)
!
!           pressure contribution to b
!
         indexf=ipnei(i)
         do j=1,numfaces
            indexf=indexf+1
            ifa=neifa(indexf)
            bv(jdof1,1)=bv(jdof1,1)
     &           -vfa(4,ifa)*xxn(1,indexf)*area(ifa)
            bv(jdof1,2)=bv(jdof1,2)
     &           -vfa(4,ifa)*xxn(2,indexf)*area(ifa)
            bv(jdof1,3)=bv(jdof1,3)
     &           -vfa(4,ifa)*xxn(3,indexf)*area(ifa)
         enddo
!            
      enddo
!
c      do i=1,nzs
c         write(*,*) 'mafillv auv,irow',i,auv(i),auv(i+nzs)
c      enddo
c      do i=1,nef
c         write(*,*) 'mafillv b adv',i,bv(i,1),adv(i)
c      enddo
!     
      return
      end
