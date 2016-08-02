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
      subroutine rhsp(nef,lakonf,ipnei,neifa,neiel,vfa,area,
     &  advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
     &  b,xxn,neq,nzs,hfa,gradpel,bp,xxi,neij,xlen,nefa,nefb)
!
!     filling the lhs and rhs to calculate the first correction to the
!     pressure p'
!
      implicit none
!
      character*8 lakonf(*)
!
      integer i,nef,jdof1,indexf,ipnei(*),j,neifa(*),
     &  neiel(*),iel,ifa,jdof2,irow(*),ielfa(4,*),nefa,nefb,
     &  ifabou(*),neq,nzs,jq(*),iel2,indexb,knownflux,numfaces,
     &  iatleastonepressurebc,j2,indexf2,neij(*)
!
      real*8 coef,vfa(0:5,*),volume(*),area(*),advfa(*),xlet(*),
     &  cosa(*),ad(*),au(*),xle(*),xxn(3,*),ap(*),b(*),bp(*),
     &  hfa(3,*),xxi(3,*),gradpel(3,*),xlen(*),bp_ifa
!
c      iatleastonepressurebc=0
!
      do i=nefa,nefb
         jdof1=i
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
!
!              diffusion
!
            indexf=indexf+1
            ifa=neifa(indexf)
            iel=neiel(indexf)
            if(iel.ne.0) then
!
!           element has a neighbor across the face
!
!           correction for non-orthogonal meshes
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
               b(jdof1)=b(jdof1)-ap(ifa)*bp_ifa
c               if(i.gt.iel) bp_ifa=-bp_ifa
            else
!
!                 external face
!
               iel2=ielfa(2,ifa)
               if(iel2.lt.0) then
                  if((ifabou(-iel2+1).ne.0).and.
     &                 (ifabou(-iel2+2).ne.0).and.
     &                 (ifabou(-iel2+3).ne.0)) then
!
!                    all velocity components given
!
                     knownflux=1
                  elseif(ifabou(-iel2+5).eq.2) then
!
!                    sliding conditions
!
                     knownflux=2
                  elseif(ifabou(-iel2+4).gt.0) then
c                     iatleastonepressurebc=1
!     
!                    pressure given
!                        
!     
!                    correction for non-orthogonal meshes
!     
                     bp_ifa=(-(gradpel(1,i)*(xxi(1,indexf)
     &                    -cosa(indexf)*xxn(1,indexf))+
     &                    gradpel(2,i)*(xxi(2,indexf)
     &                    -cosa(indexf)*xxn(2,indexf))+
     &                    gradpel(3,i)*(xxi(3,indexf)
     &                    -cosa(indexf)*xxn(3,indexf)))
     &                    *xle(indexf))
!
                     b(jdof1)=b(jdof1)-ap(ifa)*
     &                    (vfa(4,ifa)+bp_ifa)
                  endif
               endif
            endif
!     
!     save coefficients for correctvfa.f
!     
            if((iel.eq.0).or.(i.lt.iel)) then
               bp(ifa)=bp_ifa
            endif
!
            if(knownflux.eq.1) then
               b(jdof1)=b(jdof1)+vfa(5,ifa)*area(ifa)*
     &              (vfa(1,ifa)*xxn(1,indexf)+
     &              vfa(2,ifa)*xxn(2,indexf)+
     &              vfa(3,ifa)*xxn(3,indexf))
            elseif(knownflux.ne.2) then
               b(jdof1)=b(jdof1)+vfa(5,ifa)*area(ifa)*
     &              (hfa(1,ifa)*xxn(1,indexf)+
     &              hfa(2,ifa)*xxn(2,indexf)+
     &              hfa(3,ifa)*xxn(3,indexf))
            endif
         enddo
      enddo
!
!     at least one pressure bc is needed. If none is applied,
!     the last dof is set to 0
!
!     a pressure bc is only recognized if not all velocity degrees of
!     freedom are prescribed on the same face
!
c      if(iatleastonepressurebc.eq.0) b(nef)=0.d0
!
      return
      end
