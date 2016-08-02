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
      subroutine objectives(inpc,textpart,istep,istat,n,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc,nener,nobject,objectset)        
!
!     reading the input deck: *OBJECTIVE
!
!     criteria: MASS
!               STRESS
!               EIGENFREQUENCY
!               SHAPE ENERGY
!            
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 objectset(3,*)
!
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*),j,nener,nobject
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *OBJECTIVE: *OBJECTIVE can
     &only be used within a SENSITIVITY STEP'     
         call exit(201)
      endif
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
!
      
      do j=1,nobject
         i=1
            if(textpart(i)(1:11).eq.'SHAPEENERGY') then
               read(textpart(i)(1:85),'(a80)',iostat=istat) 
     &              objectset(1,j)(1:80)
                 if(textpart(2)(1:6).eq.'ELSET=') then
                     read(textpart(2)(7:85),'(a80)',iostat=istat) 
     &                    objectset(3,j)(1:80) 
                 endif
	       nener=1
            elseif(textpart(i)(1:4).eq.'MASS') then
               read(textpart(i)(1:4),'(a80)',iostat=istat) 
     &              objectset(1,j)(1:80)
                 if(textpart(2)(1:6).eq.'ELSET=') then
                     read(textpart(2)(7:85),'(a80)',iostat=istat) 
     &                    objectset(3,j)(1:80) 
                 endif     
            elseif(textpart(i)(1:6).eq.'STRESS') then
               read(textpart(i)(1:6),'(a80)',iostat=istat) 
     &              objectset(1,j)(1:80)
               do i=2,n
                 if(textpart(i)(1:7).eq.'ENTITY=') then
                     read(textpart(i)(8:85),'(a80)',iostat=istat) 
     &                    objectset(2,j)(1:80) 
                 elseif(textpart(i)(1:6).eq.'ELSET=') then
                     read(textpart(i)(7:85),'(a80)',iostat=istat) 
     &                    objectset(3,j)(1:80) 
                 endif
               enddo
            elseif(textpart(i)(1:14).eq.'EIGENFREQUENCY') then
               read(textpart(i)(1:14),'(a80)',iostat=istat) 
     &              objectset(1,j)(1:80)
            endif
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
!
      enddo
!
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
!
      return
      end

