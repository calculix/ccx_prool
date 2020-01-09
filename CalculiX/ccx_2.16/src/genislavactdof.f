!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2019 Guido Dhondt
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
c>     \brief islavactdof is the inverse of nactdof for active slave nodes:
c>     it links an active slave degree of freedom to the 
c>     corresponding slave node position in field islavnode and the
c>     global (x-y-z) degree of freedom
c>
c> @param [in] ntie       	number of contraints
c> @param [in] tieset     	name and dependent surface of tie set
c> @param [in] neq		number of active degrees of freedom
c> @param [in] nactdof 		(i,j) actual degree of freedom for direction i of node j
c> @param [in] nslavnode  	(i) for contraint i pointer into field islavnode
c> @param [in] islavact		(i) indicates, if slave node i is active (=-3 no-slave-node, =-2 no-LM-node, =-1 no-gap-node, =0 inactive node, =1 sticky node, =2 slipping/active node) 
c> @param [in] nmastnode	(i)pointer into field imastnode for contact tie i 
c> @param [in] imastnode	field storing the nodes of the master surfaces
c> @param [out] islavactdof     (i)=10*slavenodenumber+direction for active dof i
c> @param [in] islavnode	field storing the nodes of the slave surface
c> @param [in] mi		(1) max # of integration points per element (2) max degree of freedom per element
c>
      subroutine genislavactdof(ntie,tieset,neq,nactdof,nslavnode,
     &     islavact,nmastnode,imastnode,islavactdof,islavnode,mi,
     &     ithermal)
!     
!     Author : Samoela Rakotonanahary, Saskia Sitzmann
!
!     genislavactdof get the field islavactdof in order to 
!     help calculating the tangential matrices.
!     
!     islavactdof is the inverse of nactdof for active slave nodes:
!     it links an active slave degree of freedom to the 
!     corresponding slave node position in field islavnode and the
!     global (x-y-z) degree of freedom
!
      implicit none     
!     
      integer i,j,k,ntie,neq(*),node,nslavnode(*),
     &     mi(*),nactdof(0:mi(2),*),nmastnode(*),imastnode(*),
     &     islavact(*),islavactdof(*),islavnode(*),ithermal(*)
!
      character*81 tieset(3,*)
!     
!     close the contact.fbd file
!     
      close(20)
      close(30)
      close(40)
! do not change the order here: slave nodes have to be treated last since two contact definitions can share an edge     
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         do j = nmastnode(i)+1,nmastnode(i+1)
            node=imastnode(j)
            do k=1,3
               if (nactdof(k,node).le.0) cycle
               islavactdof(nactdof(k,node))=-(10*j+k)
            enddo
            if(ithermal(1).gt.1)then
               if (nactdof(0,node).le.0) cycle
               islavactdof(nactdof(0,node))=-(10*j+4)               
            endif        
         enddo
      enddo
!
      do i=1,ntie
         if(tieset(1,i)(81:81).ne.'C') cycle
         do j = nslavnode(i)+1,nslavnode(i+1)
            node=islavnode(j)
               do k=1,3
                  if (nactdof(k,node).le.0) cycle
                  islavactdof(nactdof(k,node))=10*j+k
               enddo 
            if(ithermal(1).gt.1)then
               if (nactdof(0,node).le.0) cycle
               islavactdof(nactdof(0,node))=(10*j+4)               
            endif       
         enddo
      enddo

!     
      return
      end
      
