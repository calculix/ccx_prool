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
      subroutine errorestimator(yi,yn,ipkon,kon,lakon,nk,
     &  ne,mi,ielmat,nterms)
!
!     the error in the node is calculated based on the maximum difference
!     between the max principal stress (mechanical calculations)
!     or heat flux size (thermal calculations) at the integration
!     points in the elements belonging to the node
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),kon(*),mi(*),ne,indexe,null,
     &  nonei20(3,12),nonei10(3,6),nk,i,j,k,
     &  nonei15(3,9),nopev,nterms,ishift,mint3d,
     &  m,jj,ielmat(mi(3),*),nlayer,nopeexp
!
      real*8 yi(nterms,mi(1),*),yn(nterms,*),size,wpsmin,wpsmax,
     &  absdiff,reldiff,sizemax,al(3),sizemin,firstprin,c(3,3)
!
      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
      null=0
!
      do i=1,nk
         do j=1,nterms
            yn(j,i)=0.d0
         enddo
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(1:1).eq.'F') then
            cycle
         elseif(lakonl(4:4).eq.'2') then
            nopev=8
         elseif(lakonl(4:4).eq.'8') then
            nopev=8
         elseif(lakonl(4:5).eq.'10') then
            nopev=4
         elseif(lakonl(4:4).eq.'4') then
            nopev=4
         elseif(lakonl(4:5).eq.'15') then
            nopev=6
         elseif(lakonl(4:4).eq.'6') then
            nopev=6
         else
            cycle
         endif
!
         if(lakonl(4:5).eq.'8R') then
            mint3d=1
c         elseif(lakonl(4:7).eq.'20RB') then
c            if((lakonl(8:8).eq.'R').or.(lakonl(8:8).eq.'C')) then
c               mint3d=50
c            else
c               call beamintscheme(lakonl,mint3d,ielprop(i),prop,
c     &              null,xi,et,ze,weight)
c            endif
         elseif((lakonl(4:4).eq.'8').or.
     &           (lakonl(4:6).eq.'20R')) then
            if(lakonl(7:8).eq.'LC') then
               cycle
            else
               mint3d=8
            endif
         elseif(lakonl(4:4).eq.'2') then
            mint3d=27
         elseif(lakonl(4:5).eq.'10') then
            mint3d=4
         elseif(lakonl(4:4).eq.'4') then
            mint3d=1
         elseif(lakonl(4:5).eq.'15') then
            mint3d=9
         elseif(lakonl(4:5).eq.'6') then
            mint3d=2
         elseif(lakonl(1:2).eq.'ES') then
            cycle
         endif
!
!        calculating the maximal differences of first principal
!        stress or heat flux size across the integration points
!
         absdiff=0.d0
         reldiff=0.d0
!
         if(nterms.eq.6) then
!
!           mechanical calculation: max principal stress
!
            wpsmin=1.d30
            wpsmax=-1.e30
            do j=1,mint3d
               c(1,1)=yi(1,j,i)
               c(2,2)=yi(2,j,i)
               c(3,3)=yi(3,j,i)
               c(1,2)=yi(4,j,i)
               c(1,3)=yi(5,j,i)
               c(2,3)=yi(6,j,i)
!     
!              calculate the eigenvalues
!     
               call calceigenvalues(c,al)
!     
c               if(dabs(al(3)).gt.dabs(al(1))) then
c                  firstprin=al(3)
c               else
c                  firstprin=al(1)
c               endif
               firstprin=al(3)
!     
               wpsmin=min(wpsmin,firstprin)
               wpsmax=max(wpsmax,firstprin)
!
c            write(*,*) 'errorestimator',i,j,firstprin
            enddo
            absdiff=wpsmax-wpsmin
            if(max(dabs(wpsmax),dabs(wpsmin)).lt.1.d-30) then
               reldiff=0.d0
            else
               reldiff=absdiff/(max(dabs(wpsmax),dabs(wpsmin)))
            endif
c            write(*,*) 'errorestimator',i,absdiff,reldiff
         else
!
!           thermal calculation: size of heat flux
!            
            sizemin=1.e30
            sizemax=0.d0
!
            do j=1,mint3d
               c(1,1)=yi(1,j,i)
               c(2,2)=yi(2,j,i)
               c(3,3)=yi(3,j,i)
!
               size=dsqrt(c(1,1)**2+c(2,2)**2+c(3,3)**2)
               sizemin=min(sizemin,size)
               sizemax=max(sizemax,size)
!
            enddo
            absdiff=sizemax-sizemin
            if(max(sizemax,sizemin).lt.1.d-30) then
               reldiff=0.d0
            else
               reldiff=absdiff/(max(sizemax,sizemin))
            endif
         endif
!
!        transferring the maximum to the nodes belonging to the
!        element
!
         do j=1,nopev
            yn(3,kon(indexe+j))=max(yn(3,kon(indexe+j)),absdiff)
            yn(5,kon(indexe+j))=max(yn(5,kon(indexe+j)),reldiff)
         enddo
!
      enddo
!
!        determining the field values in the midside nodes
!
      if(nterms.eq.6) then
         ishift=2
      else
         ishift=4
      endif
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') then
            nlayer=0
            do j=1,mi(3)
               if(ielmat(j,i).gt.0) then
                  nlayer=nlayer+1
               else
                  exit
               endif
            enddo
!
            if(lakonl(4:4).eq.'2') then
               nopeexp=28
            elseif(lakonl(4:5).eq.'15') then
               nopeexp=21
            endif
         endif
!
         if(lakonl(4:5).eq.'20') then
            if(lakonl(7:8).ne.'LC') then
               do j=9,20
                  do k=3,5,ishift
                     yn(k,kon(indexe+j))=(
     &                    yn(k,kon(indexe+nonei20(2,j-8)))+
     &                    yn(k,kon(indexe+nonei20(3,j-8))))/2.d0
                  enddo
               enddo
            else
               do m=1,nlayer
                  jj=20*(m-1)
                  do j=9,20
                     do k=3,5,ishift
                        yn(k,kon(indexe+nopeexp+jj+j))=(
     &                      yn(k,kon(indexe+nopeexp+jj+nonei20(2,j-8)))+
     &                      yn(k,kon(indexe+nopeexp+jj+nonei20(3,j-8))))
     &                      /2.d0
                     enddo
                  enddo
               enddo
            endif
         elseif(lakonl(4:5).eq.'10') then
            do j=5,10
               do k=3,5,ishift
                  yn(k,kon(indexe+j))=(yn(k,kon(indexe+nonei10(2,j-4)))+
     &                 yn(k,kon(indexe+nonei10(3,j-4))))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            do j=7,15
               do k=3,5,ishift
                  yn(k,kon(indexe+j))=(yn(k,kon(indexe+nonei15(2,j-6)))+
     &                 yn(k,kon(indexe+nonei15(3,j-6))))/2.d0
               enddo
            enddo
         endif
      enddo
!     
      return
      end
