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
      subroutine gaps(inpc,textpart,set,istartset,iendset,
     &  ialset,nset,nset_,nalset,nalset_,ipompc,nodempc,coefmpc,
     &  labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,lakon,ipkon,kon,nk,nk_,
     &  nodeboun,ndirboun,ikboun,ilboun,nboun,nboun_,iperturb,ne_,
     &  co,xboun,ctrl,typeboun,istep,istat,n,iline,ipol,inl,ipoinp,
     &  inp,iamboun,nam,inotr,trab,ntrans,nmethod,ipoinpc,mi,vold)
!
!     reading the input deck: *GAP
!
!     a gap between nodes a and b is formulated by a nonlinear MPC
!     linking node a and b. To simulate the gap feature an extra node
!     c is introduced. The first DOF of this node is fixed to zero by
!     a boundary SPC, the second DOF is left free. If the gap is closed
!     the first DOF of node c is used in the MPC leading to a linear,
!     tied MPC. If the gap is open, the second DOF of node c is used,
!     leading to no constraint at all.
!
      implicit none
!
      logical fixed,calcnormal
!
      character*1 typeboun(*),type,inpc(*)
      character*8 lakon(*)
      character*20 labmpc(*),label,labelb
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer istartset(*),iendset(*),ialset(*),ipompc(*),
     &  nodempc(3,*),noderef,idirref,
     &  nset,nset_,nalset,nalset_,nmpc,nmpc_,mpcfree,nk,nk_,ikmpc(*),
     &  ilmpc(*),ipkon(*),kon(*),i,node,ipos,istep,istat,n,ne_,
     &  j,k,nodeboun(*),ndirboun(*),ikboun(*),ilboun(*),iamboun(*),
     &  nboun,nboun_,key,iperturb(2),inode,iline,ipol,inl,ipoinpc(0:*),
     &  ipoinp(2,*),inp(3,*),l,index1,ibounstart,ibounend,iamplitude,
     &  nam,inotr(2,*),ntrans,nmethod,idummy,mi(*),node1,node2
!
      real*8 coefmpc(3,*),co(3,*),xboun(*),ctrl(*),xn(3),clearance,
     &  bounval,trab(7,*),dd,vold(0:mi(2),*)
!
      labelb='                    '
      fixed=.false.
      type='B'
      iamplitude=0
!
      if(istep.gt.0) then
         write(*,*) 
     &     '*ERROR reading *GAP: *GAP should be placed'
         write(*,*) '  before all step definitions'
         call exit(201)
      endif
!
!     reading the element set
!
      elset='
     &                      '
      ipos=0
!
      do i=2,n
         if(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &        '*WARNING reading *GAP: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*GAP%")
         endif
      enddo
!
!     checking whether the element set exists
!
      if(ipos.eq.0) then
         write(*,*) '*ERROR reading *GAP: no element set ',elset
         write(*,*) '       was been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*GAP%")
         call exit(201)
      endif
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *GAP: element set ',elset
         write(*,*) '  has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*GAP%")
         call exit(201)
      endif
!
!     the *GAP option implies a nonlinear geometric 
!     calculation
!
      iperturb(2)=1
      write(*,*) '*INFO reading *GAP: nonlinear geometric'
      write(*,*) '      effects are turned on'
      write(*,*)
      if(iperturb(1).eq.0) then
         iperturb(1)=2
      elseif(iperturb(1).eq.1) then
         write(*,*) '*ERROR reading *GAP: the *MPC option'
         write(*,*) '       cannot be used in a perturbation step'
         call exit(201)
      endif
!
      label='GAP                 '
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      read(textpart(1)(1:20),'(f20.0)',iostat=istat) clearance
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*GAP%")
      read(textpart(2)(1:20),'(f20.0)',iostat=istat) xn(1)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*GAP%")
      read(textpart(3)(1:20),'(f20.0)',iostat=istat) xn(2)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*GAP%")
      read(textpart(4)(1:20),'(f20.0)',iostat=istat) xn(3)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*GAP%")
!
!     check whether size of gap normal is zero; if so, the 
!     gap normal is calculated from the coordinates
!
      calcnormal=.false.
      dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
      if(dabs(dd).eq.0.d0) calcnormal=.true.
!
!     generating the gap MPC's
!
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            if(lakon(ialset(j))(1:1).ne.'G') then
               write(*,*) '*ERROR gaps: *GAP can only be used for'
               write(*,*) '       GAPUNI elements'
               write(*,*) '       Faulty element: ',ialset(j)
               call exit(201)
            endif
            index1=ipkon(ialset(j))
!
!           three terms for node 1
!
            node1=kon(index1+1)
            inode=0
            do l=1,3
               inode=inode+1
               call usermpc(ipompc,nodempc,coefmpc,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &              nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &              nboun,nboun_,inode,node1,co,label,
     &              typeboun,iperturb,noderef,idirref,xboun)
            enddo
!
!           three terms for node 2
!
            node2=kon(index1+2)
            do l=1,3
               inode=inode+1
               call usermpc(ipompc,nodempc,coefmpc,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &              nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &              nboun,nboun_,inode,node2,co,label,
     &              typeboun,iperturb,noderef,idirref,xboun)
            enddo
!
!           extra node for the gap DOF
!
            nk=nk+1
            if(nk.gt.nk_) then
               write(*,*) '*ERROR reading *GAP: increase nk_'
               call exit(201)
            endif
            node=nk
            call usermpc(ipompc,nodempc,coefmpc,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &           nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &           nboun,nboun_,inode,node,co,label,typeboun,
     &           iperturb,noderef,idirref,xboun)
!
!           calculating the gap normal
!
            if(calcnormal) then
               do l=1,3
                  xn(l)=co(l,node2)-co(l,node1)
               enddo
               dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
               if(dabs(dd).eq.0.d0) then
                  write(*,*) '*ERROR reading *GAP: gap normal cannot '
                  write(*,*) '       determined'
                  call exit(201)
               endif
               do l=1,3
                  xn(l)=xn(l)/dd
               enddo
            endif
!               
            do l=1,3
               co(l,nk)=xn(l)
            enddo
!
!           restraining the first DOF of the extra node
!
            ibounstart=1
            ibounend=1
            bounval=0.d0
            call bounadd(node,ibounstart,ibounend,bounval,
     &        nodeboun,ndirboun,xboun,nboun,nboun_,
     &        iamboun,iamplitude,nam,ipompc,nodempc,
     &        coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &        ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &        type,typeboun,nmethod,iperturb,fixed,vold,
     &        idummy,mi,labelb)
!
!           nonhomogeneous term for user MPC
!
            node=0
            call usermpc(ipompc,nodempc,coefmpc,
     &           labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &           nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &           nboun,nboun_,inode,node,co,label,typeboun,
     &           iperturb,noderef,idirref,xboun)
            co(1,nk)=clearance
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               if(lakon(k)(1:1).ne.'G') then
                  write(*,*)'*ERROR reading *GAP: *GAP can only be used'
                  write(*,*) '       for GAPUNI elements'
                  write(*,*) '       Faulty element: ',k
                  call exit(201)
               endif
               index1=ipkon(k)
!
!              three terms for node 1
!
               node1=kon(index1+1)
               inode=0
               do l=1,3
                  inode=inode+1
                  call usermpc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,inode,node1,co,label,
     &                 typeboun,iperturb,noderef,idirref,xboun)
               enddo
!
!              three terms for node 2
!
               node2=kon(index1+2)
               do l=1,3
                  inode=inode+1
                  call usermpc(ipompc,nodempc,coefmpc,
     &                 labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &                 nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &                 nboun,nboun_,inode,node2,co,label,
     &                 typeboun,iperturb,noderef,idirref,xboun)
               enddo
!
!              extra node for the gap DOF
!
               nk=nk+1
               if(nk.gt.nk_) then
                  write(*,*) '*ERROR reading *GAP: increase nk_'
                  call exit(201)
               endif
               node=nk
               call usermpc(ipompc,nodempc,coefmpc,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &              nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &              nboun,nboun_,inode,node,co,label,typeboun,
     &              iperturb,noderef,idirref,xboun)
!
!              calculating the gap normal
!
               if(calcnormal) then
                  do l=1,3
                     xn(l)=co(l,node2)-co(l,node1)
                  enddo
                  dd=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
                  if(dabs(dd).eq.0.d0) then
                     write(*,*) '*ERROR reading *GAP: gap normal cannot'
                     write(*,*) '       determined'
                     call exit(201)
                  endif
                  do l=1,3
                     xn(l)=xn(l)/dd
                  enddo
               endif
!     
               do l=1,3
                  co(l,nk)=xn(l)
               enddo
!
!              restraining the first DOF of the extra node
!
               ibounstart=1
               ibounend=1
               bounval=0.d0
               call bounadd(node,ibounstart,ibounend,bounval,
     &              nodeboun,ndirboun,xboun,nboun,nboun_,
     &              iamboun,iamplitude,nam,ipompc,nodempc,
     &              coefmpc,nmpc,nmpc_,mpcfree,inotr,trab,
     &              ntrans,ikboun,ilboun,ikmpc,ilmpc,co,nk,nk_,labmpc,
     &              type,typeboun,nmethod,iperturb,fixed,vold,idummy,
     &              mi,labelb)
!
!              nonhomogeneous term for user MPC
!
               node=0
               call usermpc(ipompc,nodempc,coefmpc,
     &              labmpc,nmpc,nmpc_,mpcfree,ikmpc,ilmpc,
     &              nk,nk_,nodeboun,ndirboun,ikboun,ilboun,
     &              nboun,nboun_,inode,node,co,label,typeboun,
     &              iperturb,noderef,idirref,xboun)
               co(1,nk)=clearance
            enddo
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end
