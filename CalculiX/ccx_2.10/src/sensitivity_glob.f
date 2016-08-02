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
      subroutine sensitivity_glob(dgdxtot,dgdxtotglob,nobject,ndesi,
     &  nodedesi,ndirdesi,nk)
               
!
!    prepares the sensitivities for the output in the frd-file      
!
      implicit none
!
      integer nobject,ndesi,nodedesi(*),ndirdesi(*),nk,numobject,
     &  numnode,numdesvar,i
!
      real*8 dgdxtot(ndesi,nobject),dgdxtotglob(4,nk,nobject)
!     
!     Loop over all nodes
      do numnode=1,nk
!     Loop over designvariables
         do numdesvar=1,ndesi
            if(nodedesi(numdesvar).eq.numnode) then
!     Loop over all objective functions
               do numobject=1,nobject    
!     loop over all coordinates     
                  do i=1,3
                     if(nodedesi(numdesvar).eq.nodedesi(numdesvar+
     &                    i-1)) then
                        dgdxtotglob(i,numnode,numobject)=
     &                       dgdxtot(numdesvar+i-1,numobject)
                     endif
                  enddo
               enddo
               exit
            endif
         enddo
      enddo
!      
      return        
      end




