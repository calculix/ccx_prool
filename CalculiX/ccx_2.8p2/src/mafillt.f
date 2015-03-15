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
      subroutine mafillt(nef,ipnei,neifa,neiel,vfa,xxn,area,
     &  au,ad,jq,irow,nzs,b,vel,umel,xlet,xle,gradtfa,xxi,
     &  body,volume,compressible,ielfa,lakon,ifabou,nbody,neq,
     &  dtimef,velo,veloo,cpfa,hcfa,cpel,gradvel,xload,gammat,xrlfa,
     &  xxj,nactdohinv,a1,a2,a3)
!
!     filling the matrix for the conservation of energy
!
      implicit none
!
      logical knownflux
!
      character*8 lakon(*)
!
      integer i,nef,jdof1,indexf,ipnei(*),j,ifa,iel,neifa(*),
     &  neiel(*),jdof2,jq(*),irow(*),nzs,compressible,ielfa(4,*),
     &  ipointer,ifabou(*),nbody,neq,indexb,numfaces,nactdohinv(*)
!
      real*8 flux,vfa(0:5,*),xxn(3,*),area(*),au(*),ad(*),b(neq),
     &  vel(nef,0:5),umel(*),xlet(*),xle(*),coef,gradtfa(3,*),
     &  xxi(3,*),body(0:3,*),volume(*),dtimef,velo(nef,0:5),
     &  veloo(nef,0:5),rhovel,constant,cpel(*),gradvel(3,3,*),
     &  cpfa(*),hcfa(*),div,xload(2,*),gammat(*),xrlfa(3,*),
     &  xxj(3,*),a1,a2,a3
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
!     convection
!
            indexf=indexf+1
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) jdof2=iel
            flux=(vfa(1,ifa)*xxn(1,indexf)+
     &           vfa(2,ifa)*xxn(2,indexf)+
     &           vfa(3,ifa)*xxn(3,indexf))
     &           *vfa(5,ifa)*cpfa(ifa)*area(ifa)
!     
            if(flux.ge.0.d0) then
!     
!     outflowing flux
!     
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &              flux,nzs)
!     centdiff
ccccccc                   b(jdof1)=b(jdof1)-(vfa(0,ifa)-vel(i,0))*flux
!end centdiff
            else
               if(iel.gt.0) then
!
!                    incoming flux from neighboring element
!
                  call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,flux,
     &                 nzs)
!centdiff
cccccc                    b(jdof1)=b(jdof1)-(vfa(0,ifa)-vel(iel,0))*flux
!end centdiff
               else
!
!                    incoming flux through boundary
!
                  if(ielfa(2,ifa).lt.0) then
                     indexb=-ielfa(2,ifa)
                     if((ifabou(indexb).ne.0).or.
     &                    (dabs(flux).lt.1.d-10)) then
                        b(jdof1)=b(jdof1)-vfa(0,ifa)*flux
                     else
                        write(*,*) '*ERROR in mafillt: the tempera-'
                        write(*,*) '       ture of an incoming flux'
                        write(*,*) '       through face ',j,'of'
                        write(*,*)'       element ',nactdohinv(i),
     &                          ' is not given'
c                        call exit(201)
                     endif
                  else
                     write(*,*) '*ERROR in mafillt: the tempera-'
                     write(*,*) '       ture of an incoming flux'
                     write(*,*) '       through face ',j,'of'
                     write(*,*)'       element ',nactdohinv(i),
     &                   ' is not given'
c                     call exit(201)
                  endif
               endif
            endif
!
!           diffusion
!
            if(iel.ne.0) then
!     
!              neighboring element
!     
               coef=hcfa(ifa)*area(ifa)/xlet(indexf)
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &              coef,nzs)
               call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof2,
     &              -coef,nzs)
!     
!              correction for non-orthogonal grid
!     
               b(jdof1)=b(jdof1)+hcfa(ifa)*area(ifa)*
     &              (gradtfa(1,ifa)*(xxn(1,indexf)-xxj(1,indexf))+
     &              gradtfa(2,ifa)*(xxn(2,indexf)-xxj(2,indexf))+
     &              gradtfa(3,ifa)*(xxn(3,indexf)-xxj(3,indexf)))
            else
!     
!              boundary; either temperature given or adiabatic
!              or outlet
!     
               knownflux=.false.
               ipointer=abs(ielfa(2,ifa))
               if(ipointer.gt.0) then
                  if(ifabou(ipointer+6).gt.0) then
!     
!     heat flux is known
!     
                     b(jdof1)=b(jdof1)-xload(1,ifabou(ipointer+6))
                     knownflux=.true.
!     
                  elseif((ifabou(ipointer).ne.0).or.
     &                    (ifabou(ipointer+1).ne.0).or.
     &                    (ifabou(ipointer+2).ne.0).or.
     &                    (ifabou(ipointer+3).ne.0)) then
!     
!                    temperature given or no outlet:
!                    temperature is assumed fixed
!     
                     coef=hcfa(ifa)*area(ifa)/xle(indexf)
                     call add_sm_fl_as(au,ad,jq,irow,jdof1,jdof1,
     &                    coef,nzs)
                     b(jdof1)=b(jdof1)+coef*vfa(0,ifa)
                  else
!     
!                     outlet: no diffusion
!     
                  endif
               endif
!     
!              correction for non-orthogonal grid
!     
               if(.not.knownflux) then
                  b(jdof1)=b(jdof1)+hcfa(ifa)*area(ifa)*
     &                 (gradtfa(1,ifa)*(xxn(1,indexf)-xxi(1,indexf))+
     &                 gradtfa(2,ifa)*(xxn(2,indexf)-xxi(2,indexf))+
     &                 gradtfa(3,ifa)*(xxn(3,indexf)-xxi(3,indexf)))
               endif
            endif
         enddo
!     
!        -p*div(v) term
!     
         if(compressible.eq.1) then
            div=gradvel(1,1,i)+gradvel(2,2,i)+gradvel(3,3,i)
            b(jdof1)=b(jdof1)-vel(i,4)*div*volume(i)
         endif
!
!           viscous dissipation
!
         b(jdof1)=b(jdof1)+umel(i)*volume(i)*
     &        (2.d0*(gradvel(1,1,i)**2+gradvel(2,2,i)**2+
     &        gradvel(3,3,i)**2)+
     &        (gradvel(1,2,i)+gradvel(2,1,i))**2+
     &        (gradvel(1,3,i)+gradvel(3,1,i))**2+
     &        (gradvel(2,3,i)+gradvel(3,2,i))**2)
         if(compressible.eq.1) then
            b(jdof1)=b(jdof1)-2.d0*umel(i)*volume(i)/3.d0*div**2
         endif
!     
!           body heat source and body sources
!     
         rhovel=vel(i,5)*volume(i)
!
         if(nbody.gt.0) then
c               b(jdof1)=b(jdof1)+rhovel*(body(0,i)+body(1,i)*vel(i,1)+
c     &                body(2,i)*vel(i,2)+body(3,i)*vel(i,3))
            b(jdof1)=b(jdof1)+rhovel*body(0,i)
c                  write(*,*) 'body source ',body(0,i),body(1,i),
c     &             body(2,i),body(3,i)
         endif
!
!           transient term
!
         constant=rhovel*cpel(i)/dtimef
         b(jdof1)=b(jdof1)-(a2*velo(i,0)+a3*veloo(i,0))*constant
         constant=a1*constant
         call add_sm_fl(au,ad,jq,irow,jdof1,jdof1,constant,nzs)
!     
      enddo
!     
      return
      end
