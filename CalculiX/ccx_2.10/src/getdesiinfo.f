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
     &  mi,nactdof,ndesi,ndirdesi,nodedesi,ntie,tieset)    
!
!     reading the input deck: *DESIGNVARIABLES
!
      implicit none
!
      character*81 setname,nodeset
      character*81 set(*)
      character*81 tieset(3,*)
!
      integer mi(*),istartset(*),iendset(*),ialset(*),ndesi,
     &  node,nodedesi(*),ndirdesi(*),nset,ntie,i,j,k,l,
     &  ipos,nactdof(0:mi(2),*) 
!
      nodeset='
     &                      '
      ipos=0
      ndesi=0
!
!     Search for the set name
!      
      do i=1,ntie
         if(tieset(1,i)(81:81).eq.'D') then
            setname=tieset(2,i)
            ipos=index(setname,' ')
            setname(ipos:ipos)='N'
         endif
      enddo 
!
!     Check for the existence of the name
!
      if(setname.eq.' ') then
        write(*,*) '*ERROR in getdesiinfo: name of node set '
        write(*,*) '  has not yet been defined. '
        call exit(201)
      endif
!
!     Search the name of the node set in "set(i)" and
!     assign the nodes of the set to the appropriate variables
!
      do j=1,nset
         if(setname.eq.set(j)) then  
            do k=istartset(j),iendset(j)
               node=ialset(k)
               do l=1,3
                   if(nactdof(l,node).ne.0) then
                     ndesi=ndesi+1
                     nodedesi(ndesi)=node
                     ndirdesi(ndesi)=l
                   endif   
               enddo
            enddo              
         endif
      enddo 
!
      return
      end

