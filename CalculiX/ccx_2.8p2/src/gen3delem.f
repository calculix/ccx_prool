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
      subroutine gen3delem(kon,ipkon,lakon,ne,ipompc,nodempc,coefmpc,
     &  nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,ikboun,ilboun,nboun,
     &  nboun_,nodeboun,ndirboun,xboun,iamboun,nam,
     &  inotr,trab,nk,nk_,iponoel,inoel,iponor,xnor,thicke,thickn,
     &  knor,istep,offset,t0,t1,ikforc,ilforc,rig,nforc,
     &  nforc_,nodeforc,ndirforc,xforc,iamforc,nelemload,sideload,
     &  nload,ithermal,ntrans,co,ixfree,ikfree,inoelfree,iponoelmax,
     &  iperturb,tinc,tper,tmin,tmax,ctrl,typeboun,nmethod,nset,set,
     &  istartset,iendset,ialset,prop,ielprop,vold,mi,nkon,ielmat,
     &  icomposite,t0g,t1g,idefforc,iamt1)
!
!     generates three-dimensional elements:
!         for isochoric elements
!         for plane stress
!         for plane strain
!         for plate and shell elements
!         for beam elements
!
      implicit none
!
      logical isochoric
!
      character*1 typeboun(*)
      character*8 lakon(*)
      character*20 labmpc(*),sideload(*),label
      character*81 set(*)
!
      integer ipompc(*),nodempc(3,*),nmpc,nmpc_,mpcfree,ikmpc(*),
     &  ilmpc(*),kon(*),ipkon(*),ne,mpc,indexe,i,j,k,node,idof,
     &  id,mpcfreeold,ikboun(*),ilboun(*),nboun,nboun_,kflag,idummy,
     &  iterm(500),nterm,neigh(7,8),l,m,nodeboun(*),ndirboun(*),nk,
     &  nk_,index,iponoel(*),inoel(3,*),inoelfree,istep,nmpcold,
     &  ikforc(*),ilforc(*),nodeforc(2,*),ndirforc(*),iamforc(*),
     &  nelemload(*),nforc,nforc_,ithermal(2),nload,iamboun(*),
     &  ntrans,inotr(2,*),nam,iponoelmax,iperturb,numnod,itransaxial,
     &  rig(*),nmethod,nset,istartset(*),iendset(*),ialset(*),nkon,
     &  ielprop(*),idir,indexref,indexold,idofold,idold,indexnew,
     &  idofnew,idnew,ksol,lsol,nmpc0,nmpc01,nmpcdif,mi(*),nope,
     &  ielmat(mi(3),*),iponor(2,*),knor(*),ixfree,ikfree,icomposite,
     &  idefforc(*),idim,iamt1(*)
!
      real*8 coefmpc(*),thicke(mi(3),*),xnor(*),thickn(2,*),tinc,
     &  tper,tmin,t0g(2,*),t1g(2,*),e1(3),e2(3),xt1(3),
     &  tmax,offset(2,*),t0(*),t1(*),xforc(*),trab(7,*),co(3,*),b(3,3),
     &  xboun(*),pi,ctrl(*),prop(*),vold(0:mi(2),*),xlag(3,20),
     &  xeul(3,20),a(3,3),xi,et,ze,coloc(3,8),xj
!
      data neigh /1,9,2,12,4,17,5,2,9,1,10,3,18,6,
     &            3,11,4,10,2,19,7,4,11,3,12,1,20,8,
     &            5,13,6,16,8,17,1,6,13,5,14,7,18,2,
     &            7,15,8,14,6,19,3,8,15,7,16,5,20,4/
!
      data coloc /-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,
     &            -1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1./
!
      isochoric=.false.
      pi=4.d0*datan(1.d0)
!
!     catalogueing the element per node relationship for shell/beam
!     elements and transferring the nodal thickness to the elements
!
!     inoelfree=1 means that there is at least one 1D or 2D element
!     in the structure. Otherwise inoelfree=0.
!
      if((istep.eq.1).and.(inoelfree.eq.1)) then
!
!        shift of the connectivity for composite elements
!
         if(icomposite.eq.1) then
            call changekon(ne,ipkon,lakon,mi,nkon,thicke,ielmat,kon)
         endif
!
         itransaxial=0
!
         do i=1,ne
            if(ipkon(i).lt.0) cycle
            if((lakon(i)(1:2).ne.'C3').and.(lakon(i)(1:1).ne.'D').and.
     &         (lakon(i)(1:1).ne.'G').and.(lakon(i)(1:1).ne.'E')) then
!
!              number of nodes belonging to the element
!
               if(lakon(i)(1:3).eq.'B31') then
                  numnod=2
               elseif(lakon(i)(1:1).eq.'B') then
                  numnod=3
               elseif((lakon(i)(2:2).eq.'3').or.
     &                (lakon(i)(4:4).eq.'3')) then
                  numnod=3
               elseif((lakon(i)(2:2).eq.'4').or.
     &                (lakon(i)(4:4).eq.'4')) then
                  numnod=4
               elseif((lakon(i)(2:2).eq.'6').or.
     &                (lakon(i)(4:4).eq.'6')) then
                  numnod=6
               elseif((lakon(i)(2:2).eq.'8').or.
     &                (lakon(i)(4:4).eq.'8')) then
                  numnod=8
               endif
!
c               if(lakon(i)(1:1).eq.'B') then
c                  numnod=3
c               elseif((lakon(i)(2:2).eq.'6').or.
c     &                (lakon(i)(4:4).eq.'6')) then
c                  numnod=6
c               else
c                  numnod=8
c               endif
               indexe=ipkon(i)
               do j=1,numnod
                  node=kon(indexe+j)
                  iponoelmax=max(iponoelmax,node)
                  inoel(1,inoelfree)=i
                  inoel(2,inoelfree)=j
                  inoel(3,inoelfree)=iponoel(node)
                  iponoel(node)=inoelfree
                  inoelfree=inoelfree+1
                  if(lakon(i)(1:2).ne.'CA') then
                     if(thickn(1,node).gt.0.d0)
     &                    thicke(1,indexe+j)=thickn(1,node)
                     if(thickn(2,node).gt.0.d0)
     &                    thicke(2,indexe+j)=thickn(2,node)
                  endif
                  if(thicke(1,indexe+j).le.0.d0) then
                     if(lakon(i)(1:1).eq.'C') then
                        thicke(1,indexe+j)=1.d0
                     else
                        write(*,*)'*ERROR in gen3delem: first thickness'
                        write(*,*)'       in node ',j,' of element ',i
                        write(*,*)'       is zero'
                        call exit(201)
                     endif
                  endif
                  if((lakon(i)(1:1).eq.'B').and.
     &                 (thicke(2,indexe+j).le.0.d0)) then
                     write(*,*) '*ERROR in gen3delem: second thickness'
                     write(*,*)'       in node ',j,' of beam element ',i
                     write(*,*)'       is zero'
                     call exit(201)
                  endif
               enddo
            endif
         enddo
c!
c!        checking whether any rotational degrees of freedom are fixed
c!        by SPC's, MPC's or loaded by bending moments or torques
c!        in the end, rig(i)=0 if no rigid knot is defined in node i,
c!        else rig(i)=the rotational node of the knot. The value -1 is
c!        a dummy. 
c!
c         do i=1,nboun
c            if(ndirboun(i).gt.4) rig(nodeboun(i))=-1
c         enddo
c         do i=1,nforc
c            if(ndirforc(i).gt.4) rig(nodeforc(1,i))=-1
c         enddo
c         do i=1,nmpc
c            index=ipompc(i)
c            do
c               if(index.eq.0) exit
c               if(nodempc(2,index).gt.4) then
c                  rig(nodempc(1,index))=-1
c               endif
c               index=nodempc(3,index)
c            enddo
c         enddo
!
!     calculating the normals in nodes belonging to shells/beams
!
         nmpcold=nmpc
!
         call gen3dnor(nk,nk_,co,iponoel,inoel,iponoelmax,kon,ipkon,
     &     lakon,ne,thicke,offset,iponor,xnor,knor,rig,iperturb,tinc,
     &     tper,tmin,tmax,ctrl,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikmpc,ilmpc,labmpc,ikboun,ilboun,nboun,nboun_,
     &     nodeboun,ndirboun,xboun,iamboun,typeboun,nam,ntrans,inotr,
     &     trab,ikfree,ixfree,nmethod,ithermal,istep,mi,icomposite,
     &     ielmat,vold)
!
      endif
!
      if(istep.eq.1) then
c!
c!        incompressible elements
c!
c         nmpc0=nmpc
c         nmpc01=nmpc0+1
c         do i=1,ne
c            if(ipkon(i).lt.0) cycle
c            if(lakon(i)(1:7).eq.'C3D20RI') then
c               isochoric=.true.
c               indexe=ipkon(i)
c!
c               do j=1,20
c                  node=kon(indexe+j)
c                  do k=1,3
c                     xlag(k,j)=co(k,node)
c                     xeul(k,j)=xlag(k,j)+vold(k,node)
c                  enddo
c               enddo
c!
c               do j=1,8
c                  node=kon(indexe+j)
c                  mpc=0
c                  label(1:9)='ISOCHORIC'
c                  write(label(10:20),'(i11)') node
c                  nmpcdif=nmpc-nmpc0
c                  call cident20(labmpc(nmpc01),label,nmpcdif,id)
c                  id=id+nmpc0
c                  if(id.gt.0) then
c                     if(labmpc(id).eq.label) then
c                        mpc=id
c                     endif
c                  endif
c!
c!                 new MPC: look for suitable dependent dof
c!
c                  if(mpc.eq.0) then
c                     mpc=id+1
c                     ksol=0
c                     loop: do k=1,7
c                        do l=1,3
c                           idof=8*(kon(indexe+neigh(k,j))-1)+l
c!
c!                          check for SPC's using the same DOF
c!
c                           call nident(ikboun,idof,nboun,id)
c                           if(id.gt.0) then
c                              if(ikboun(id).eq.idof) cycle
c                           endif
c!
c!                          check for MPC's using the same DOF
c!     
c                           call nident(ikmpc,idof,nmpc,id)
c                           if(id.gt.0) then
c                              if(ikmpc(id).eq.idof) cycle
c                           endif
c!
c                           ksol=k
c                           lsol=l
c                           exit loop
c                        enddo
c                     enddo loop
c!
c!                    no mpc available
c!
c                     if(ksol.eq.0) then
c                        write(*,*) 
c     &                       '*WARNING in gen3delem: no free DOF in'
c                        write(*,*) 
c     &                       '         node ',node,' for isochoric'
c                        write(*,*) '         MPC application'
c                        cycle
c                     endif
c!
c!                    new mpc
c!
c                     nmpc=nmpc+1
c                     if(nmpc.gt.nmpc_) then
c                        write(*,*) '*ERROR in gen3delem: increase nmpc_'
c                        call exit(201)
c                     endif
c                     do l=1,nmpc
c                        if(ilmpc(l).ge.mpc) ilmpc(l)=ilmpc(l)+1
c                     enddo
c                     do l=nmpc,id+2,-1
c                        ikmpc(l)=ikmpc(l-1)
c                        ilmpc(l)=ilmpc(l-1)
c                     enddo
c                     ikmpc(id+1)=idof
c                     ilmpc(id+1)=mpc
c                     do l=nmpc,mpc+1,-1
c                        ipompc(l)=ipompc(l-1)
c                        labmpc(l)=labmpc(l-1)
c                     enddo
c!
c                     labmpc(mpc)(1:9)='ISOCHORIC           '
c                     write(labmpc(mpc)(10:20),'(i11)') node
c!
c!                    terms of the node itself and its neighbors
c!
c                     ipompc(mpc)=mpcfree
c                     do l=lsol,3
c                        nodempc(1,mpcfree)=kon(indexe+neigh(ksol,j))
c                        nodempc(2,mpcfree)=l
c                        mpcfree=nodempc(3,mpcfree)
c                     enddo
c!
c                     do k=ksol+1,7
c                        do l=1,3
c                           nodempc(1,mpcfree)=kon(indexe+neigh(k,j))
c                           nodempc(2,mpcfree)=l
c                           mpcfree=nodempc(3,mpcfree)
c                        enddo
c                     enddo
c!
c                     do k=1,ksol-1
c                        do l=1,3
c                           nodempc(1,mpcfree)=kon(indexe+neigh(k,j))
c                           nodempc(2,mpcfree)=l
c                           mpcfree=nodempc(3,mpcfree)
c                        enddo
c                     enddo
c!
c                     do l=1,lsol-1
c                        nodempc(1,mpcfree)=kon(indexe+neigh(ksol,j))
c                        nodempc(2,mpcfree)=l
c                        mpcfree=nodempc(3,mpcfree)
c                     enddo
c!
c!                 add nonhomogeneous term
c!                  
c                     nk=nk+1
c                     if(nk.gt.nk_) then
c                        write(*,*) '*ERROR in gen3delem: increase nk_'
c                        call exit(201)
c                     endif
c                     mpcfreeold=mpcfree
c                     mpcfree=nodempc(3,mpcfree)
c                     nodempc(1,mpcfreeold)=nk
c                     nodempc(2,mpcfreeold)=1
c                     nodempc(3,mpcfreeold)=0
c                     idof=8*(nk-1)+1
c                     call nident(ikboun,idof,nboun,id)
c                     nboun=nboun+1
c                     if(nboun.gt.nboun_) then
c                        write(*,*)'*ERROR in gen3delem: increase nboun_'
c                        call exit(201)
c                     endif
c                     nodeboun(nboun)=nk
c                     ndirboun(nboun)=1
c                     typeboun(nboun)='I'
c                     do l=nboun,id+2,-1
c                        ikboun(l)=ikboun(l-1)
c                        ilboun(l)=ilboun(l-1)
c                     enddo
c                     ikboun(id+1)=idof
c                     ilboun(id+1)=nboun
c!
c                  else
c!
c                     indexref=nodempc(3,nodempc(3,ipompc(mpc)))
c                     index=nodempc(3,indexref)
c                     nterm=0
c                     do
c                        if(index.eq.0) exit
c                        nterm=nterm+1
c                        if(nterm.gt.500) then
c                           write(*,*) '*ERROR in gen3delem:'
c                           write(*,*) '       increase nterm_'
c                           call exit(201)
c                        endif
c                        iterm(nterm)=
c     &                       8*(nodempc(1,index)-1)+nodempc(2,index)
c                        index=nodempc(3,index)
c                     enddo
c                     kflag=1
c                     call isortii(iterm,idummy,nterm,kflag)
c!
c                     do k=2,7
c                        do l=1,3
c                           m=8*(kon(indexe+neigh(k,j))-1)+l
c                           call nident(iterm,m,nterm,id)
c                           if(id.ne.0) then
c                              if(iterm(id).eq.m) then
c                                 cycle
c                              endif
c                           endif
c                           mpcfreeold=mpcfree
c                           mpcfree=nodempc(3,mpcfree)
c                           nodempc(3,mpcfreeold)=nodempc(3,indexref)
c                           nodempc(3,indexref)=mpcfreeold
c                           nodempc(1,mpcfreeold)=kon(indexe+neigh(k,j))
c                           nodempc(2,mpcfreeold)=l
c                        enddo
c                     enddo
c!
c                  endif
c!
c                  xi=coloc(1,j)
c                  et=coloc(2,j)
c                  ze=coloc(3,j)
c!     
c                  call deuldlag(xi,et,ze,xlag,xeul,xj,a)
c!     
c                  b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
c                  b(1,2)=a(3,1)*a(2,3)-a(2,1)*a(3,3)
c                  b(1,3)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
c                  b(2,1)=a(3,2)*a(1,3)-a(1,2)*a(3,3)
c                  b(2,2)=a(1,1)*a(3,3)-a(3,1)*a(1,3)
c                  b(2,3)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
c                  b(3,1)=a(1,2)*a(2,3)-a(2,2)*a(1,3)
c                  b(3,2)=a(2,1)*a(1,3)-a(1,1)*a(2,3)
c                  b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
c!     
c                  index=ipompc(mpc)
c                  do
c                     if(nodempc(3,index).eq.0) then
c                        coefmpc(index)=1.d0
c                        idof=8*(nodempc(1,index)-1)
c     &                       +nodempc(2,index)
c                        call nident(ikboun,idof,nboun,id)
c                        xboun(ilboun(id))=xboun(ilboun(id))+
c     &                       a(1,1)*b(1,1)+a(1,2)*b(1,2)+a(1,3)*b(1,3)
c     &                       -1.d0/xj
c                        exit
c                     else
c                        node=nodempc(1,index)
c                        idir=nodempc(2,index)
c                        do k=1,7
c                           if(kon(indexe+neigh(k,j)).eq.node) then
c                              if(k.eq.1) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)+1.5d0*
c     &                                  (xi*b(1,1)+et*b(1,2)+ze*b(1,3))
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)+1.5d0*
c     &                                  (xi*b(2,1)+et*b(2,2)+ze*b(2,3))
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)+1.5d0*
c     &                                  (xi*b(3,1)+et*b(3,2)+ze*b(3,3))
c                                 endif
c                              elseif(k.eq.2) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*xi*b(1,1)
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*xi*b(2,1)
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*xi*b(3,1)
c                                 endif
c                              elseif(k.eq.3) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*xi*b(1,1)
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*xi*b(2,1)
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*xi*b(3,1)
c                                 endif
c                              elseif(k.eq.4) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*et*b(1,2)
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*et*b(2,2)
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*et*b(3,2)
c                                 endif
c                              elseif(k.eq.5) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*et*b(1,2)
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*et*b(2,2)
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*et*b(3,2)
c                                 endif
c                              elseif(k.eq.6) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*ze*b(1,3)
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*ze*b(2,3)
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    -2.d0*ze*b(3,3)
c                                 endif
c                              elseif(k.eq.7) then
c                                 if(idir.eq.1) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*ze*b(1,3)
c                                 elseif(idir.eq.2) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*ze*b(2,3)
c                                 elseif(idir.eq.3) then
c                                    coefmpc(index)=coefmpc(index)
c     &                                    +0.5d0*ze*b(3,3)
c                                 endif
c                              endif
c                              exit
c                           endif
c                        enddo
c                     endif
c                     index=nodempc(3,index)
c                  enddo
c!
c               enddo
c            endif
c         enddo
!
!        if there is any plane stress, plane strain or axisymmetric
!        element the structure should lie in the z=0 plane
!
         if(inoelfree.ne.0) then
            do i=1,ne
               if(ipkon(i).lt.0) cycle
               if((lakon(i)(1:2).eq.'CP').or.
     &              (lakon(i)(1:2).eq.'CA')) then
                  indexe=ipkon(i)
                  read(lakon(i)(4:4),'(i1)') nope
c                  if(lakon(i)(4:4).eq.'6') then
c                     nope=6
c                  else
c                     nope=8
c                  endif
                  do j=1,nope
                     node=kon(indexe+j)
                     if(dabs(co(3,node)).gt.0.d0) then
                        write(*,*) '*ERROR in gen3delem. The structure'
                        write(*,*) '       contains plane stress, plane'
                        write(*,*) '       strain or axisymmetric'
                        write(*,*) '       elements and should lie in '
                        write(*,*) '       the z=0 plane. This is at'
                        write(*,*) '       least not the case for node',
     &                       node
                        call exit(201)
                     endif
                  enddo
               endif
            enddo
         endif
!
!        1D and 2D elements
!
         if(inoelfree.ne.0) then
            do i=1,ne
               if(ipkon(i).lt.0) cycle
               if((lakon(i)(1:2).eq.'CP').or.
     &              (lakon(i)(1:1).eq.'S').or.
     &              (lakon(i)(1:2).eq.'CA')) then
!
                 call gen3dfrom2d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &           thicke,offset,ntrans,inotr,trab,ikboun,ilboun,nboun,
     &           nboun_,nodeboun,ndirboun,xboun,iamboun,typeboun,ipompc,
     &           nodempc,coefmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,
     &           nk,nk_,co,rig,nmethod,iperturb,ithermal,mi,nam,
     &           icomposite,ielmat,vold)
!            
              elseif(lakon(i)(1:1).eq.'B') then
                 call gen3dfrom1d(i,kon,ipkon,lakon,ne,iponor,xnor,knor,
     &                thicke,ntrans,inotr,trab,nk,nk_,co,offset,mi)
              endif
!     
              if(lakon(i)(1:4).eq.'CPE3') then
                 lakon(i)(1:7)='C3D6  E'
              elseif(lakon(i)(1:5).eq.'CPE4R') then
                 lakon(i)(1:7)='C3D8R E'
              elseif(lakon(i)(1:4).eq.'CPE4') then
                 lakon(i)(1:7)='C3D8  E'
              elseif(lakon(i)(1:4).eq.'CPE6') then
                 lakon(i)(1:7)='C3D15 E'
              elseif(lakon(i)(1:5).eq.'CPE8R') then
                 lakon(i)(1:7)='C3D20RE'
              elseif(lakon(i)(1:4).eq.'CPE8') then
                 lakon(i)(1:7)='C3D20 E'
              elseif(lakon(i)(1:4).eq.'CPS3') then
                 lakon(i)(1:7)='C3D6  S'
              elseif(lakon(i)(1:5).eq.'CPS4R') then
                 lakon(i)(1:7)='C3D8R S'
              elseif(lakon(i)(1:4).eq.'CPS4') then
                 lakon(i)(1:7)='C3D8  S'
              elseif(lakon(i)(1:4).eq.'CPS6') then
                 lakon(i)(1:7)='C3D15 S'
              elseif(lakon(i)(1:5).eq.'CPS8R') then
                 lakon(i)(1:7)='C3D20RS'
              elseif(lakon(i)(1:4).eq.'CPS8') then
                 lakon(i)(1:7)='C3D20 S'
              elseif(lakon(i)(1:4).eq.'CAX3') then
                 lakon(i)(1:7)='C3D6  A'
              elseif(lakon(i)(1:5).eq.'CAX4R') then
                 lakon(i)(1:7)='C3D8R A'
              elseif(lakon(i)(1:4).eq.'CAX4') then
                 lakon(i)(1:7)='C3D8  A'
              elseif(lakon(i)(1:4).eq.'CAX6') then
                 lakon(i)(1:7)='C3D15 A'
              elseif(lakon(i)(1:5).eq.'CAX8R') then
                 lakon(i)(1:7)='C3D20RA'
              elseif(lakon(i)(1:4).eq.'CAX8') then
                 lakon(i)(1:7)='C3D20 A'
              elseif(lakon(i)(1:2).eq.'S3') then
                 lakon(i)(1:7)='C3D6  L'
              elseif(lakon(i)(1:3).eq.'S4R') then
                 lakon(i)(1:7)='C3D8R L'
              elseif(lakon(i)(1:2).eq.'S4') then
                 lakon(i)(1:7)='C3D8I L'
              elseif(lakon(i)(1:2).eq.'S6') then
                 lakon(i)(1:7)='C3D15 L'
              elseif(lakon(i)(1:3).eq.'S8R') then
                 lakon(i)(1:7)='C3D20RL'
              elseif(lakon(i)(1:2).eq.'S8') then
                 lakon(i)(1:7)='C3D20 L'
              elseif(lakon(i)(1:4).eq.'B31R') then
                 lakon(i)(1:7)='C3D8R B'
              elseif(lakon(i)(1:3).eq.'B31') then
                 lakon(i)(1:7)='C3D8I B'
              elseif(lakon(i)(1:4).eq.'B32R') then
                 lakon(i)(1:7)='C3D20RB'
              elseif(lakon(i)(1:3).eq.'B32') then
                 lakon(i)(1:7)='C3D20 B'
              endif
           enddo
c     Bernhardi start
        endif
        do i=1,ne
           if(lakon(i)(1:5).eq.'C3D8I') then
              call genmodes(i,kon,ipkon,lakon,ne,nk,nk_,co)
           endif
        enddo
c     Bernhardi end
c        do i=1,68
c           write(*,*) kon(i),",",co(1,kon(i)),",",co(2,kon(i)),
c     &          ",",co(3,kon(i))
c        enddo
!     
!     check whether the coefficient of the dependent
!     terms in ISOCHORIC MPC's is not zero
!
         if(isochoric) then
            do i=1,nmpc
               if(labmpc(i)(1:9).ne.'ISOCHORIC') cycle
               index=ipompc(i)
               if(dabs(coefmpc(index)).gt.1.d-10) cycle
!
!              coefficient of dependent term is zero: rearranging
!              the MPC
!
               indexold=index
               idofold=8*(nodempc(1,index)-1)+nodempc(2,index)
               call nident(ikmpc,idofold,nmpc,idold)
               do j=idold,nmpc-1
                  ikmpc(j)=ikmpc(j+1)
                  ilmpc(j)=ilmpc(j+1)
               enddo
               indexref=index
               index=nodempc(3,index)
!
               do
                  if(index.eq.0) then
                     write(*,*) '*ERROR in gen3delem: coefficient'
                     write(*,*) '       of dependent term is zero'
                     write(*,*) '       and no other DOF is available'
                     call exit(201)
                  endif
                  if(dabs(coefmpc(index)).gt.1.d-10) then
                     idofnew=8*(nodempc(1,index)-1)+nodempc(2,index)
!
!                    check whether DOF is not used in SPC
!
                     call nident(ikboun,idofnew,nboun,idnew)
                     if(idnew.gt.0) then
                        if(ikboun(idnew).eq.idofnew) then
                           indexref=index
                           index=nodempc(3,index)
                           cycle
                        endif
                     endif
!
!                    check whether DOF is not used in MPC
!
                     call nident(ikmpc,idofnew,nmpc,idnew)
                     if(idnew.gt.0) then
                        if(ikmpc(idnew).eq.idofnew) then
                           indexref=index
                           index=nodempc(3,index)
                           cycle
                        endif
                     endif
!
!                    DOF is OK: take it as dependent term
!
                     do j=nmpc,idnew+2,-1
                        ikmpc(j)=ikmpc(j-1)
                        ilmpc(j)=ilmpc(j-1)
                     enddo
                     ikmpc(idnew+1)=idofnew
                     ilmpc(idnew+1)=i
!
                     indexnew=index
                     index=nodempc(3,index)
                     ipompc(i)=indexnew
                     nodempc(3,indexnew)=indexold
                     nodempc(3,indexref)=index
                     exit
                  endif
                  indexref=index
                  index=nodempc(3,index)
               enddo
            enddo
         endif
!
!        filling the new KNOT MPC's (needs the coordinates
!        of the expanded nodes)
!
         if(inoelfree.ne.0) then
            call fillknotmpc(co,ipompc,nodempc,coefmpc,labmpc,
     &           nmpc,nmpcold,mpcfree,idim,e1,e2,xt1)
            call gen3dprop(prop,ielprop,iponoel,inoel,iponoelmax,kon,
     &           ipkon,lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,
     &           nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,
     &           trab,nam,nk,nk_,co,nmethod,iperturb)
         endif
!     
      endif
!
!        generating MPC's to connect shells and beams with solid
!        elements
!
      if((inoelfree.ne.0).and.(istep.eq.1)) then
         call gen3dconnect(kon,ipkon,lakon,ne,iponoel,inoel,
     &     iponoelmax,rig,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,
     &     nmpc_,mpcfree,ikmpc,ilmpc,labmpc)
      endif
!
      if(inoelfree.ne.0) then
!
!           multiplying existing boundary conditions
!
         call gen3dboun(ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,
     &     xboun,iamboun,typeboun,iponoel,inoel,iponoelmax,kon,ipkon,
     &     lakon,ne,iponor,xnor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &     mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,nk_,
     &     co,nmethod,iperturb,istep,vold,mi)
!
!        updating the nodal surfaces: establishing links between the user
!        defined nodes and the newly generated nodes (mid-nodes
!        for 2d elements, mean of corner nodes for 1d elements)
!
         if(istep.eq.1) then
            call gen3dsurf(iponoel,inoel,iponoelmax,kon,ipkon,
     &        lakon,ne,iponor,knor,ipompc,nodempc,coefmpc,nmpc,nmpc_,
     &        mpcfree,ikmpc,ilmpc,labmpc,rig,ntrans,inotr,trab,nam,nk,
     &        nk_,co,nmethod,iperturb,nset,set,istartset,iendset,ialset,
     &        ikboun,ilboun,nboun,nboun_,nodeboun,ndirboun,xboun,
     &        iamboun,typeboun,mi,vold)
         endif
!
!        updating the MPCs: establishing links between the user
!        defined nodes and the newly generated nodes (mid-nodes
!        for 2d elements, mean of corner nodes for 1d elements)
!
         if(istep.eq.1) then
            call gen3dmpc(ipompc,nodempc,coefmpc,nmpc,nmpc_,mpcfree,
     &        ikmpc,ilmpc,labmpc,iponoel,inoel,iponoelmax,kon,ipkon,
     &        lakon,ne,iponor,xnor,knor,rig)
         endif
!
!        updating the temperatures
!
         if(ithermal(1).gt.0) then
            call gen3dtemp(iponoel,inoel,iponoelmax,kon,ipkon,lakon,ne,
     &           iponor,xnor,knor,t0,t1,thicke,offset,rig,nk,nk_,co,
     &           istep,ithermal,vold,mi,t0g,t1g,nam,iamt1)
         endif
!
!        updating the concentrated loading
!
         call gen3dforc(ikforc,ilforc,nforc,nforc_,nodeforc,
     &     ndirforc,xforc,iamforc,ntrans,inotr,trab,rig,ipompc,nodempc,
     &     coefmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,labmpc,iponoel,inoel,
     &     iponoelmax,kon,ipkon,lakon,ne,iponor,xnor,knor,nam,nk,nk_,
     &     co,thicke,nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,
     &     iamboun,typeboun,xboun,nmethod,iperturb,istep,vold,mi,
     &     idefforc)
      endif
!
      return
      end


