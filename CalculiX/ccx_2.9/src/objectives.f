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
     &        ipoinp,inp,ipoinpc,objective)        
!
!     reading the input deck: *OBJECTIVE
!
!     criteria: MASS
!               STRESS
!               EIGENFREQUENCY
!               DISPLACEMENT
!            
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
      character*81 objective(3)
!
      integer istep,istat,n,key,i,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*)
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *OBJECTIVE: *OBJECTIVE can
     &only be used within a SENSITIVITY STEP'     
         call exit(201)
      endif
!
      do i=2,n
         if(textpart(i)(1:9).eq.'CRITERIA=') then
            read(textpart(i)(10:85),'(a80)',iostat=istat) 
     &           objective(1)(1:80)
         elseif(textpart(i)(1:7).eq.'ENTITY=') then
            read(textpart(i)(8:85),'(a80)',iostat=istat) 
     &           objective(2)(1:80) 
         elseif(textpart(i)(1:3).eq.'MIN') then
            read(textpart(i)(1:85),'(a80)',iostat=istat) 
     &           objective(3)(1:80) 
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
!
      return
      end

