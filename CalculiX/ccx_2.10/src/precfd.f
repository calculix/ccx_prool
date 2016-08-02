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
      subroutine precfd(ne,ipkon,kon,lakon,ipnei,neifa,neiel,ipoface,
     &  nodface,ielfa,nkonnei,nface,ifaext,nfaext,
     &  isolidsurf,nsolidsurf,set,nset,istartset,iendset,ialset,
     &  vel,vold,mi,neij,nef,nactdoh,ipkonf,lakonf,ielmatf,ielmat,
     &  ielorienf,ielorien,norien)
!
!     gathering topological information for computational fluid 
!     dynamics calculations
!
      implicit none
!
      character*8 lakon(*),lakonf(*)
      character*81 set(*),noset
!
      integer ne,ipkon(*),ipnei(*),ipoface(*),nodface(5,*),neifa(*),
     &  ielfa(4,*),nkonnei,nface,i,j,k,index,indexe,neiel(*),ithree,
     &  nfaext,ifaext(*),isolidsurf(*),nsolidsurf,indexf,
     &  nset,istartset(*),iendset(*),ialset(*),iaux,kflag,ifour,iel2,
     &  ifaceq(8,6),ifacet(7,4),ifacew(8,5),kon(*),nodes(4),iel1,j2,
     &  indexold,ifree,ifreenew,ifreenei,mi(*),neij(*),ifreenei2,nef,
     &  nactdoh(*),ipkonf(*),ielmatf(mi(3),*),ielmat(mi(3),*),nf(5),
     &  nope,ielorien(mi(3),*),ielorienf(mi(3),*),norien,jopposite8(6),
     &  jopposite6(5)
!
      real*8 vel(nef,0:5),vold(0:mi(2),*)
!
!     nodes belonging to the element faces
!
      data ifaceq /4,3,2,1,11,10,9,12,
     &            5,6,7,8,13,14,15,16,
     &            1,2,6,5,9,18,13,17,
     &            2,3,7,6,10,19,14,18,
     &            3,4,8,7,11,20,15,19,
     &            4,1,5,8,12,17,16,20/
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
      data ifacew /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             4,6,3,1,12,15,9,13/
      data jopposite6 /2,1,0,0,0/
      data jopposite8 /2,1,5,6,3,4/
      data nf /3,3,4,4,4/
!
!     transfer element information from ipkon, lakon, ielmat
!     and ielorien to ipkonf, lakonf, ielmatf and ielorienf. 
!     The latter fields were created
!     for CFD applications in which the elements are numerated in
!     ascending order (without gaps)
!
      do i=1,ne
c         write(*,*) 'precfd ',i,nactdoh(i)
         if(nactdoh(i).ne.0) then
            ipkonf(nactdoh(i))=ipkon(i)
            lakonf(nactdoh(i))=lakon(i)
            do j=1,mi(3)
               ielmatf(j,nactdoh(i))=ielmat(j,i)
               if(norien.gt.0) ielorienf(j,nactdoh(i))=ielorien(j,i)
            enddo
         endif
      enddo
!
      kflag=1
      ithree=3
      ifour=4
!
!     determining the external element faces of the fluid mesh 
!     the faces are catalogued by the three lowes nodes numbers
!     in ascending order. ipoface(i) points to a face for which
!     node i is the lowest node and nodface(1,ipoface(i)) and
!     nodface(2,ipoface(i)) are the next lower ones. 
!     nodface(3,ipoface(i)) contains the element number,
!     nodface(4,ipoface(i)) the face number and nodface(5,ipoface(i))
!     is a pointer to the next surface for which node i is the
!     lowest node; if there are no more such surfaces the pointer
!     has the value zero
!     An external element face is one which belongs to one element
!     only
!
      ifree=1
      ifreenei=0
!
      do i=1,6*nef
         nodface(5,i)=i+1
      enddo
      do i=1,nef
         indexe=ipkonf(i)
         if(lakonf(i)(4:4).eq.'8') then
!
!           hex element
!
            ipnei(i)=ifreenei
            do j=1,6
               do k=1,4
                  nodes(k)=kon(indexe+ifaceq(k,j))
               enddo
               call isortii(nodes,iaux,ifour,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifreenei=ifreenei+1
                     neifa(ifreenei)=ifree
                     ielfa(1,ifree)=i
                     ielfa(4,ifree)=j
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     ifreenei=ifreenei+1
!
!                    completing the facial info in neifa
!
                     neifa(ifreenei)=index
                     ielfa(2,index)=i
!
!                    the neighboring elements to the face are i and iel2
!                    with local face number for the face of j and j2
!
                     iel2=ielfa(1,index)
                     j2=ielfa(4,index)
!
!                    completing the neighboring info for (element i,side j)
!
                     neiel(ifreenei)=iel2
                     neij(ifreenei)=j2
!
!                    completing the neighboring info for (element iel2,side j2)
!
                     ifreenei2=ipnei(iel2)+j2
                     neiel(ifreenei2)=i
                     neij(ifreenei2)=j
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         else if(lakonf(i)(4:4).eq.'6') then
!
!           wedge element
!
            ipnei(i)=ifreenei
            do j=1,5
               do k=1,nf(j)
                  nodes(k)=kon(indexe+ifacew(k,j))
               enddo
               call isortii(nodes,iaux,nf(j),kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifreenei=ifreenei+1
                     neifa(ifreenei)=ifree
                     ielfa(1,ifree)=i
                     ielfa(4,ifree)=j
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     ifreenei=ifreenei+1
!
!                    completing the facial info in neifa
!
                     neifa(ifreenei)=index
                     ielfa(2,index)=i
!
!                    the neighboring elements to the face are i and iel2
!                    with local face number for the face of j and j2
!
                     iel2=ielfa(1,index)
                     j2=ielfa(4,index)
!
!                    completing the neighboring info for (element i,side j)
!
                     neiel(ifreenei)=iel2
                     neij(ifreenei)=j2
!
!                    completing the neighboring info for (element iel2,side j2)
!
                     ifreenei2=ipnei(iel2)+j2
                     neiel(ifreenei2)=i
                     neij(ifreenei2)=j
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         else
!
!           tet element
!
            ipnei(i)=ifreenei
            do j=1,4
               do k=1,3
                  nodes(k)=kon(indexe+ifacet(k,j))
               enddo
               call isortii(nodes,iaux,ithree,kflag)
               indexold=0
               index=ipoface(nodes(1))
               do
!
!                 adding a surface which has not been 
!                 catalogued so far
!
                  if(index.eq.0) then
                     ifreenew=nodface(5,ifree)
                     nodface(1,ifree)=nodes(2)
                     nodface(2,ifree)=nodes(3)
                     nodface(3,ifree)=i
                     nodface(4,ifree)=j
                     nodface(5,ifree)=ipoface(nodes(1))
                     ipoface(nodes(1))=ifree
                     ifreenei=ifreenei+1
                     neifa(ifreenei)=ifree
                     ielfa(1,ifree)=i
                     ielfa(4,ifree)=j
                     ifree=ifreenew
                     exit
                  endif
!
!                 removing a surface which has already
!                 been catalogued
!
                  if((nodface(1,index).eq.nodes(2)).and.
     &               (nodface(2,index).eq.nodes(3))) then
                     ifreenei=ifreenei+1
!
!                    completing the facial info in neifa
!
                     neifa(ifreenei)=index
                     ielfa(2,index)=i
!
!                    the neighboring elements to the face are i and iel2
!                    with local face number for the face of j and j2
!
                     iel2=ielfa(1,index)
                     j2=ielfa(4,index)
!
!                    completing the neighboring info for (element i,side j)
!
                     neiel(ifreenei)=iel2
                     neij(ifreenei)=j2
!
!                    completing the neighboring info for (element iel2,side j2)
!
                     ifreenei2=ipnei(iel2)+j2
                     neiel(ifreenei2)=i
                     neij(ifreenei2)=j
                     exit
                  endif
                  indexold=index
                  index=nodface(5,index)
               enddo
            enddo
         endif
      enddo
!
      nkonnei=ifreenei
      nface=ifree-1
!
!     catalogueing external faces
!
      nfaext=0
!
      do i=1,nface
         if(ielfa(2,i).ne.0) cycle
         nfaext=nfaext+1
         ifaext(nfaext)=i
c         j=ielfa(4,i)
c         iel1=ielfa(1,i)
c         indexf=ipnei(iel1)
c         j=jopposite(j)
         iel1=ielfa(1,i)
         if(lakonf(iel1)(4:4).eq.'8') then
            j=ielfa(4,i)
            j=jopposite8(j)
         elseif(lakonf(iel1)(4:4).eq.'6') then
            j=ielfa(4,i)
            j=jopposite6(j)
            if(j.eq.0) cycle
         else
            cycle
         endif
         indexf=ipnei(iel1)
         ielfa(3,i)=neiel(indexf+j)
      enddo
!
c      do i=1,ne
c         k=nactdoh(i)
c         if(k.eq.0) cycle
c         write(*,*)'precfd neiel',i,k,ipnei(k),(neiel(ipnei(k)+j),j=1,6)
c         write(*,*)'precfd neij',i,k,ipnei(k),(neij(ipnei(k)+j),j=1,6)
c         write(*,*)'precfd neifa',i,k,ipnei(k),(neifa(ipnei(k)+j),j=1,6)
c      enddo
c      do k=1,nface
c         write(*,*) 'precfd ielfa',k,(ielfa(j,k),j=1,4)
c      enddo
c      do k=1,nfaext
c         write(*,*) 'precfd ifaext',k,ifaext(k)
c      enddo
c      write(*,*)
!
!     faces belonging to solid surfaces
!
      nsolidsurf=0
      noset(1:13)='SOLIDSURFACET'
      do i=1,nset
         if(set(i)(1:13).eq.noset(1:13)) exit
      enddo
      if(i.gt.nset) then
         write(*,*) '*WARNING in precfd: facial surface SOLID SURFACE '
         write(*,*) '         has not been defined.'
         write(*,*)
      else
c         nsolidsurf=0
         do j=istartset(i),iendset(i)
            nsolidsurf=nsolidsurf+1
            isolidsurf(nsolidsurf)=ialset(j)
         enddo
         call isortii(isolidsurf,iaux,nsolidsurf,kflag)
      endif
c      do i=1,nsolidsurf
c         write(*,*) 'precfd solidsurf ',i,isolidsurf(i)
c      enddo
!
!     initial conditions: element values
!     vel was initialized at allocation
!
      do i=1,nef
         indexe=ipkonf(i)
         if(lakonf(i)(4:4).eq.'8') then
            nope=8
         elseif(lakonf(i)(4:4).eq.'6') then
            nope=6
         else
            nope=4
         endif
         do j=0,4
            do k=1,nope
               vel(i,j)=vel(i,j)+vold(j,kon(indexe+k))
            enddo
            vel(i,j)=vel(i,j)/nope
         enddo
      enddo
c      do i=1,nef
c         write(*,*) 'precfd vel ',i,(vel(j,i),j=1,3)
c      enddo
!
      return
      end
