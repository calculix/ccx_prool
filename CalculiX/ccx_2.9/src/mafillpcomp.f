c!
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
      subroutine mafillpcomp(nef,lakon,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,nzs,hfa,gradpel,bp,xxi,neij,
     &  xlen,cosb,ielmat,mi,a1,a2,a3,velo,veloo,dtimef,shcon,
     &  ntmat_,vel,nactdohinv,xrlfa)
!
!     filling the lhs and rhs to calculate p
!
      implicit none
!
      character*8 lakon(*)
!
      integer i,nef,jdof1,indexf,ipnei(*),j,neifa(*),iel3,
     &  neiel(*),iel,ifa,jdof2,irow(*),ielfa(4,*),compressible,
     &  ifabou(*),neq,jq(*),iel2,indexb,knownflux,indexf2,
     &  iatleastonepressurebc,j2,neij(*),nzs,numfaces,imat,
     &  mi(*),ielmat(mi(3),*),ntmat_,nactdohinv(*),knownpressure
!
      real*8 coef,vfa(0:5,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(*),au(*),xle(*),xxn(3,*),ap(*),b(*),cosb(*),
     &  hfa(3,*),gradpel(3,*),bp(*),xxi(3,*),xlen(*),r,a1,a2,a3,
     &  flux,constant,velo(nef,0:5),veloo(nef,0:5),dtimef,
     &  shcon(0:3,ntmat_,*),vel(nef,0:5),dd,convec,fluxisobar,
     &  coef1,coef3,xrlfa(3,*),coef2,gamma,coefp,coefn,xmach
!
      iatleastonepressurebc=0
!
      do i=1,nef
c         if(i.eq.2840) write(*,*) 'mafillpcomp ',vel(i,4)
!
c         xmach=dsqrt((vel(i,1)**2+vel(i,2)**2+vel(i,3)**2)/
c     &            (1.4*r*vel(i,0)))
c         write(*,*) dsqrt(vel(i,1)**2+vel(i,2)**2+vel(i,3)**2),
c     &        r,vel(i,0)
         jdof1=i
         imat=ielmat(1,i)
         r=shcon(3,1,imat)
         indexf=ipnei(i)
         if(lakon(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakon(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
            knownflux=0
            knownpressure=0
!     
!     diffusion
!     
            indexf=indexf+1
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) then
               jdof2=iel
               coef=vfa(5,ifa)*(volume(i)+volume(iel))*area(ifa)/
     &              (advfa(ifa)*2.d0*xlet(indexf)*cosb(indexf))
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &              coef,nzs)
c               write(*,*) 'mafillpcomp ',jdof1,jdof1,coef
c               coef1=coef
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,
     &              -coef,nzs)
c               write(*,*) 'mafillpcomp ',jdof1,jdof2,-coef
c               coef2=-coef
               b(jdof1)=b(jdof1)+coef*(velo(jdof2,4)-velo(jdof1,4))
               convec=coef*(velo(jdof2,4)-velo(jdof1,4))
!     
!     correction for non-orthogonal meshes
!     
               j2=neij(indexf)
               indexf2=ipnei(iel)+j2
               bp(ifa)=((gradpel(1,iel)*(xxi(1,indexf2)
     &              -cosa(indexf2)*xxn(1,indexf2))+
     &              gradpel(2,iel)*(xxi(2,indexf2)
     &              -cosa(indexf2)*xxn(2,indexf2))+
     &              gradpel(3,iel)*(xxi(3,indexf2)
     &              -cosa(indexf2)*xxn(3,indexf2)))
     &              *xle(indexf2)
     &              -(gradpel(1,i)*(xxi(1,indexf)
     &              -cosa(indexf)*xxn(1,indexf))+
     &              gradpel(2,i)*(xxi(2,indexf)
     &              -cosa(indexf)*xxn(2,indexf))+
     &              gradpel(3,i)*(xxi(3,indexf)
     &              -cosa(indexf)*xxn(3,indexf)))
     &              *xle(indexf))
               b(jdof1)=b(jdof1)+coef*bp(ifa)
               convec=convec+coef*bp(ifa)
c               convec=-convec/vfa(4,ifa)
!
!              following line is correct if the temperature
!              changes from one pressure correction iteration to
!              the next 
!
               convec=-convec/(vfa(5,ifa)*r*vfa(0,ifa))
c               write(*,*) 'mafillpcomp ',vfa(4,ifa),
c     &             (vfa(5,ifa)*r*vfa(0,ifa))
               if(i.gt.iel) bp(ifa)=-bp(ifa)
            else
               iel2=ielfa(2,ifa)
               if(iel2.lt.0) then
                  if((ifabou(-iel2+1).ne.0).and.
     &                 (ifabou(-iel2+2).ne.0).and.
     &                 (ifabou(-iel2+3).ne.0)) then
!     
!     all velocity components given
!     
                     knownflux=1
                  elseif(ifabou(-iel2+5).eq.2) then
!
!                    sliding conditions
!
                     knownflux=2
                  elseif(ifabou(-iel2+4).ne.0) then
                     knownpressure=1
                     iatleastonepressurebc=1
!     
!     pressure given (only if not all velocity
!     components are given)
!     
                     coef=vfa(5,ifa)*volume(i)*area(ifa)/
     &                    (advfa(ifa)*xle(indexf)*cosa(indexf))
                     call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &                    coef,nzs)
c               write(*,*) 'mafillpcomp ',jdof1,jdof1,coef
c                     coef1=coef
                     b(jdof1)=b(jdof1)+coef*vfa(4,ifa)
     &                                -coef*velo(jdof1,4)
                     convec=coef*vfa(4,ifa)-coef*velo(jdof1,4)
!     
!     correction for non-orthogonal meshes
!     
                     bp(ifa)=(-(gradpel(1,i)*(xxi(1,indexf)
     &                    -cosa(indexf)*xxn(1,indexf))+
     &                    gradpel(2,i)*(xxi(2,indexf)
     &                    -cosa(indexf)*xxn(2,indexf))+
     &                    gradpel(3,i)*(xxi(3,indexf)
     &                    -cosa(indexf)*xxn(3,indexf)))
     &                    *xle(indexf))
                     b(jdof1)=b(jdof1)+coef*bp(ifa)
                     convec=convec+coef*bp(ifa)
c                     convec=-convec/vfa(4,ifa)
!
!              following line is correct if the temperature
!              changes from one pressure correction iteration to
!              the next 
!
                     convec=-convec/(vfa(5,ifa)*r*vfa(0,ifa))
                  endif
               endif
            endif
!     
!     save the coefficient for correctvfa.f
!     
            ap(ifa)=coef
!
!           convection
!
!           flux
!
c            if(knownflux.eq.1) then
c               flux=area(ifa)*
c     &             (vfa(1,ifa)*xxn(1,indexf)+
c     &              vfa(2,ifa)*xxn(2,indexf)+
c     &              vfa(3,ifa)*xxn(3,indexf))
c            endif
            if(knownflux.ne.2) then
               flux=area(ifa)*
     &             (vfa(1,ifa)*xxn(1,indexf)+
     &              vfa(2,ifa)*xxn(2,indexf)+
     &              vfa(3,ifa)*xxn(3,indexf))
            else
               flux=0.d0
            endif
!
!           flux based on constant pressure
!
            if(knownflux.eq.0) then
               fluxisobar=area(ifa)*
     &             (hfa(1,ifa)*xxn(1,indexf)+
     &              hfa(2,ifa)*xxn(2,indexf)+
     &              hfa(3,ifa)*xxn(3,indexf))
            endif
!     
!           rhs
!
            if(knownflux.eq.0) then
               b(jdof1)=b(jdof1)-vfa(5,ifa)*fluxisobar
            elseif(knownflux.eq.1) then
               b(jdof1)=b(jdof1)-vfa(5,ifa)*flux
            endif
!
c            write(*,*) 'xmach ',xmach
c            if(xmach.lt.1.d0) then
c               coef=0.d0
c            else
               if(knownflux.eq.0) then
c                  coef=convec
!     
!              following line leads to oscillations in the solution
!              (only for subsonic and transonic solutions)
!
c                  coef=fluxisobar/(r*vfa(0,ifa))+convec
c                 coef=fluxisobar/(r*vfa(0,ifa))
                  coef=0.d0
               elseif(knownflux.eq.1) then
c                  coef=flux/(r*vfa(0,ifa))
                  coef=0.d0
               else
                  coef=0.d0
               endif
c            endif
!
ccccc
C            if(iel.ne.0) then
C               if(ielfa(1,ifa).eq.i) then
C                  coefp=coef*xrlfa(1,ifa)
C                  coefn=coef*xrlfa(2,ifa)
C               else
C                  coefp=coef*xrlfa(2,ifa)
C                  coefn=coef*xrlfa(1,ifa)
C               endif
C            endif
cccc
c            if(flux.ge.0.d0) then
C            gamma=0.5d0
            if(coef.ge.0.d0) then
!     
!     outflowing flux
!     
C               if(iel.gt.0) then
C                  coefp=gamma*coefp+(1.d0-gamma)*coef
C                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
C     &                 coefp,nzs)
C                  coefn=gamma*coefn
C                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,
C     &                 coefn,nzs)
C               else
C                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
C     &              coef,nzs)
C               endif
c
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &              coef,nzs)
c
c               retarded central difference     
c               b(jdof1)=b(jdof1)-(vfa(4,ifa)-vel(i,4))*coef
c
            else
               if(iel.gt.0) then
!     
!                    incoming flux from neighboring element
!
C                  coefp=gamma*coefp
C                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,coefp,
C     &                 nzs)
C                  coefn=gamma*coefn+(1.d0-gamma)*coef
C                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,coefn,
C     &                 nzs)
c
                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,coef,
     &                 nzs)
!
c
c               retarded central difference     
c                  b(jdof1)=b(jdof1)-(vfa(4,ifa)-vel(iel,4))*coef
c
               elseif(knownpressure.eq.0) then
c???                  iel3=ielfa(3,ifa)
c                  if(iel3.gt.0) then
c                     coef1=coef*xrlfa(1,ifa)
c                     call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,coef1,
c     &                    nzs)
c                     coef3=coef*xrlfa(3,ifa)
c                     call add_sm_fl_as(au,ad,jq,irow,jdof1,iel3,coef3,
c     &                    nzs)
c                  else
                     call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,coef,
     &                    nzs)
c                  endif
               endif
            endif
!
         enddo
!
!        transient term
!
         a1=1.d0/dtimef
         a2=-1.d0/dtimef
         a3=0.d0/dtimef
C         a1=1.5d0/dtimef
C         a2=-2.d0/dtimef
C         a3=0.5d0/dtimef
         constant=volume(i)/(r*vel(i,0))
         b(jdof1)=b(jdof1)-
     &        (a1*velo(i,4)+a2*velo(i,4)+a3*veloo(i,4))*constant
         constant=a1*constant
         call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,constant,nzs)
!
      enddo
!     
      return
      end
