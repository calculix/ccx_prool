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
      subroutine getdesiinfo(set,istartset,iendset,ialset,nset,
     &  mi,nactdof,ndesi,nodedesi,ntie,tieset,itmp,nmpc,nodempc,
     &  ipompc)    
!
!     reading the input deck: *DESIGNVARIABLES
!
      implicit none
!
      character*81 setname
      character*81 set(*)
      character*81 tieset(3,*)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ndesi,
     &  node,nodedesi(*),nset,ntie,i,j,k,l,nmpc,nodempc(3,*),
     &  nactdof(0:mi(2),*),itmp(*),ntmp,index,id,ipompc(*) 
!
      setname(1:1)=' '
      ndesi=0
!
!     Search for the set name
!      
      do i=1,ntie
         if(tieset(1,i)(81:81).eq.'D') then
            setname=tieset(2,i)
         endif
      enddo 
!
!     Check for the existence of the name
!
      if(setname(1:1).eq.' ') then
        write(*,*) '*ERROR in getdesiinfo: name of node set '
        write(*,*) '  has not yet been defined. '
        call exit(201)
      endif
!
!     catalogue all nodes (dependent and independent) which
!     belong to MPC's
!
      ntmp=0
      do i=1,nmpc
         index=ipompc(i)
         do
            if(index.eq.0) exit
            node=nodempc(1,index)
            call nident(itmp,node,ntmp,id)
            if(id.gt.0) then
               if(itmp(id).eq.node) then
                  index=nodempc(3,index)
                  cycle
               endif
            endif
            ntmp=ntmp+1
            do j=ntmp,id+2,-1
               itmp(j)=itmp(j-1)
            enddo
            itmp(id+1)=node
            index=nodempc(3,index)
         enddo
      enddo
!
!     Search the name of the node set in "set(i)" and
!     assign the nodes of the set to the appropriate variables
!
      do i=1,nset
         if(setname.eq.set(i)) then  
            loop1: do j=istartset(i),iendset(i)
               if(ialset(j).gt.0) then
                  node=ialset(j)
!
!                 check for SPC-constraints: if a node is constrained in all dofs
!                 it is removed from the design node set
!
                  do l=1,3
c                     if(nactdof(l,node).le.0) cycle loop1
                     if(nactdof(l,node).gt.0) exit
                     if(l.eq.3) then
                        write(*,*) '*WARNING in getdesiinfo:'
                        write(*,*) '         node ',node,' has no'
                        write(*,*) '         active dofs and'
                        write(*,*) '         is removed from the set'
                        write(*,*) '         of design variables'
                        cycle loop1
                     endif
                  enddo
!
!                 check for MPC-constraints
!
                  call nident(itmp,node,ntmp,id)
                  if(id.gt.0) then
                     if(itmp(id).eq.node) then
                        write(*,*) '*WARNING in getdesiinfo:'
                        write(*,*) '         node ',node,' is subject'
                        write(*,*) '         to MPC-constraints and'
                        write(*,*) '         is removed from the set'
                        write(*,*) '         of design variables'
                        cycle loop1
                     endif
                  endif
!
                  ndesi=ndesi+1
                  nodedesi(ndesi)=node
               else
                  k=ialset(j-2)
                  loop2: do
                     k=k-ialset(j)
                     if(k.ge.ialset(j-1)) exit
!
!                 check for SPC-constraints: if a node is constrained in all dofs
!                 it is removed from the design node set
!
                     do l=1,3
                        if(nactdof(l,k).gt.0) exit
                        if(l.eq.3) then
                           write(*,*) '*WARNING in getdesiinfo:'
                           write(*,*) '         node ',k,' has no'
                           write(*,*) '         active dofs and'
                           write(*,*) '         is removed from the set'
                           write(*,*) '         of design variables'
                           cycle loop2
                        endif
                     enddo
!
!                 check for MPC-constraints
!
                     call nident(itmp,k,ntmp,id)
                     if(id.gt.0) then
                        if(itmp(id).eq.k) then
                           write(*,*) '*WARNING in getdesiinfo:'
                           write(*,*) '         node ',k,' is subject'
                           write(*,*) '         to MPC-constraints and'
                           write(*,*) '         is removed from the set'
                           write(*,*) '         of design variables'
                           cycle loop2
                        endif
                     endif
!
                     ndesi=ndesi+1
                     nodedesi(ndesi)=k
                  enddo loop2
               endif
            enddo loop1
         endif
      enddo 
!
      return
      end

