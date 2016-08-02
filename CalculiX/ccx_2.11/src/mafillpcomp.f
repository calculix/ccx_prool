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
      subroutine mafillpcomp(nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,neq,nzs,hfa,gradpel,bp,xxi,neij,
     &  xlen,cosb,ielmatf,mi,a1,a2,a3,velo,veloo,dtimef,shcon,
     &  ntmat_,vel,nactdohinv,xrlfa,flux,nefa,nefb)
!
!     filling the lhs and rhs to calculate p
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,nef,jdof1,indexf,ipnei(*),j,neifa(*),iel3,
     &  neiel(*),iel,ifa,jdof2,irow(*),ielfa(4,*),compressible,
     &  ifabou(*),neq,jq(*),iel2,indexb,knownflux,indexf2,
     &  j2,neij(*),nzs,numfaces,imat,nefa,nefb,
     &  mi(*),ielmatf(mi(3),*),ntmat_,nactdohinv(*),knownpressure
!
      real*8 coef,vfa(0:5,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(neq),au(nzs),xle(*),xxn(3,*),ap(*),b(neq),cosb(*),
     &  hfa(3,*),gradpel(3,*),bp(*),xxi(3,*),xlen(*),r,a1,a2,a3,
     &  xflux,constant,velo(nef,0:5),veloo(nef,0:5),dtimef,
     &  shcon(0:3,ntmat_,*),vel(nef,0:5),dd,convec,fluxisobar,
     &  coef1,coef3,xrlfa(3,*),coef2,gamma,coefp,coefn,xmach,
     &  flux(*),bp_ifa,aa(8,8)
!
      intent(in) nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,jq,irow,ielfa,ifabou,xle,
     &  xxn,nzs,hfa,gradpel,xxi,neij,
     &  xlen,cosb,ielmatf,mi,a1,a2,a3,velo,veloo,dtimef,shcon,
     &  ntmat_,vel,nactdohinv,xrlfa,flux
!
      intent(inout) ad,au,b,ap,bp
!
      do i=nefa,nefb
         jdof1=i
         imat=ielmatf(1,i)
         r=shcon(3,1,imat)
         indexf=ipnei(i)
         if(lakonf(i)(4:4).eq.'8') then
            numfaces=6
         elseif(lakonf(i)(4:4).eq.'6') then
            numfaces=5
         else
            numfaces=4
         endif
         do j=1,numfaces
            knownflux=0
            knownpressure=0
            convec=0
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
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,
     &              -coef,nzs)
c               write(*,*) 'mafillpcomp1',i,j,jdof1,jdof2,-coef
               b(jdof1)=b(jdof1)+coef*(vel(jdof2,4)-vel(jdof1,4))
               convec=coef*(vel(jdof2,4)-vel(jdof1,4))
!     
!     correction for non-orthogonal meshes
!     
               j2=neij(indexf)
               indexf2=ipnei(iel)+j2
               bp_ifa=((gradpel(1,iel)*(xxi(1,indexf2)
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
               b(jdof1)=b(jdof1)+coef*bp_ifa
               convec=convec+coef*bp_ifa
!
!              following line is correct if the temperature
!              changes from one pressure correction iteration to
!              the next 
!
               convec=-convec/(vfa(5,ifa)*r*vfa(0,ifa))
c               if(i.gt.iel) bp_ifa=-bp_ifa
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
!     
!     pressure given (only if not all velocity
!     components are given)
!     
                     coef=vfa(5,ifa)*volume(i)*area(ifa)/
     &                    (advfa(ifa)*xle(indexf)*cosa(indexf))
                     call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &                    coef,nzs)
                     b(jdof1)=b(jdof1)+coef*vfa(4,ifa)
     &                                -coef*vel(jdof1,4)
                     convec=coef*vfa(4,ifa)-coef*vel(jdof1,4)
!     
!     correction for non-orthogonal meshes
!     
                     bp_ifa=(-(gradpel(1,i)*(xxi(1,indexf)
     &                    -cosa(indexf)*xxn(1,indexf))+
     &                    gradpel(2,i)*(xxi(2,indexf)
     &                    -cosa(indexf)*xxn(2,indexf))+
     &                    gradpel(3,i)*(xxi(3,indexf)
     &                    -cosa(indexf)*xxn(3,indexf)))
     &                    *xle(indexf))
                     b(jdof1)=b(jdof1)+coef*bp_ifa
                     convec=convec+coef*bp_ifa
!
!              following line is correct if the temperature
!              changes from one pressure correction iteration to
!              the next 
!
                     convec=-convec/(vfa(5,ifa)*r*vfa(0,ifa))
c                  else
c
c                    check!
c
c                     convec=0.d0
                  endif
               endif
            endif
!     
!     save coefficients for correctvfa.f
!     
            if((iel.eq.0).or.(i.lt.iel)) then
               ap(ifa)=coef
               bp(ifa)=bp_ifa
            endif
!     
!     save the coefficient for correctvfa.f
!     
c            ap(ifa)=coef
!
!           convection
!
!           flux
!
            if(knownflux.eq.1) then
               xflux=flux(indexf)
c               xflux=area(ifa)*vfa(5,ifa)*
c               xflux=area(ifa)*
c     &             (vfa(1,ifa)*xxn(1,indexf)+
c     &              vfa(2,ifa)*xxn(2,indexf)+
c     &              vfa(3,ifa)*xxn(3,indexf))
c               write(*,*) 'mafillpcomp ',i,j,ifa
c               write(*,*) vfa(5,ifa)
c               write(*,*) vfa(1,ifa)
c               write(*,*) vfa(2,ifa)
c               write(*,*) vfa(3,ifa)
c               write(*,*) 'mafillpcomp ',xflux,flux(indexf)
            elseif(knownflux.eq.2) then
               xflux=0.d0
            endif
c            write(*,*) 'mafillpcomp8 ',knownflux,xflux
!
!           flux based on constant pressure
!
            if(knownflux.eq.0) then
               fluxisobar=area(ifa)*
     &             (hfa(1,ifa)*xxn(1,indexf)+
     &              hfa(2,ifa)*xxn(2,indexf)+
     &              hfa(3,ifa)*xxn(3,indexf))
c               write(*,*) 'mafillpcomp5 ',
c     &             i,j,jdof1,jdof2,hfa(1,ifa),hfa(2,ifa),hfa(3,ifa)
c               write(*,*) 'mafillpcomp6 ',
c     &             i,j,jdof1,jdof2,xxn(1,indexf),xxn(2,indexf),
c     &                 xxn(3,indexf)
            endif
c            write(*,*) 'mafillpcomp9 ',knownflux,fluxisobar
!     
!           rhs
!
            if(knownflux.eq.0) then
               b(jdof1)=b(jdof1)-vfa(5,ifa)*fluxisobar
            elseif(knownflux.eq.1) then
               b(jdof1)=b(jdof1)-xflux
            endif
!
            if(knownflux.eq.0) then
!     
!              following line leads to oscillations in the solution
!              (only for subsonic and transonic solutions)
!
               coef=fluxisobar/(r*vfa(0,ifa))+convec
c               write(*,*) 'mafillpcomp3 ',
c     &             i,j,jdof1,jdof2,fluxisobar,r,vfa(0,ifa)
c               write(*,*) 'mafillpcomp4 ',
c     &             i,j,jdof1,jdof2,convec
            elseif(knownflux.eq.1) then
               coef=xflux/(r*vfa(0,ifa)*vfa(5,ifa))
            else
               coef=0.d0
            endif
c            write(*,*) 'mafillpcomp10 ',knownflux,coef
!
            if(coef.ge.0.d0) then
!     
!     outflowing flux
!     
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
                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,coef,
     &                 nzs)
c               write(*,*) 'mafillpcomp2',i,j,jdof1,jdof2,coef
!
c
c               retarded central difference     
c                  b(jdof1)=b(jdof1)-(vfa(4,ifa)-vel(iel,4))*coef
c
               elseif(knownpressure.eq.0) then
                     call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,coef,
     &                    nzs)
               endif
            endif
!
         enddo
!
!        transient term
!
c         a1=1.d0/dtimef
c         a2=-1.d0/dtimef
c         a3=0.d0/dtimef
c         constant=volume(i)/(r*vel(i,0))
c         b(jdof1)=b(jdof1)-
c     &        (a1*vel(i,5)+a2*velo(i,5)+a3*veloo(i,5))*volume(i)
         b(jdof1)=b(jdof1)-
     &        (vel(i,5)-velo(i,5))*volume(i)/dtimef
c         constant=a1*constant
         constant=volume(i)/(r*vel(i,0)*dtimef)
         call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,constant,nzs)
!
      enddo
!   
c      write(*,*) 'mafillpcomp '
c      do i=1,8
c         do j=1,8
c            aa(i,j)=0.d0
c         enddo
c      enddo
c      do i=1,8
c         do j=jq(i),jq(i+1)-1
c            write(*,*) i,irow(j)
c            aa(irow(j),i)=au(j)
c            aa(i,irow(j))=au(j+nzs)
c         enddo
c         aa(i,i)=ad(i)
c      enddo
c      do i=1,8
c         write(*,'(9(1x,e11.4))') (aa(i,j),j=1,8),b(i)
c      enddo
      return
      end
