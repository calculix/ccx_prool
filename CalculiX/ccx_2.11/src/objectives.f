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
     &  inp(3,*),ipoinpc(0:*),nener,nobject,k,ipos
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *OBJECTIVE: *OBJECTIVE can
     &only be used within a SENSITIVITY STEP'     
         call exit(201)
      endif
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      do
         if(textpart(1)(1:12).eq.'DISPLACEMENT') then
            nobject=nobject+1
            objectset(1,nobject)(1:12)='DISPLACEMENT'
            do k=13,81
               objectset(1,nobject)(k:k)=' '
            enddo
            if(textpart(2)(1:5).eq.'NSET=') then
               read(textpart(2)(6:85),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               objectset(3,nobject)(ipos:ipos)='N'
            endif     
         elseif(textpart(1)(1:14).eq.'EIGENFREQUENCY') then
            nobject=nobject+1
            objectset(1,nobject)(1:14)='EIGENFREQUENCY'
            do k=15,81
               objectset(1,nobject)(k:k)=' '
            enddo
         elseif(textpart(1)(1:4).eq.'MASS') then
            nobject=nobject+1
            objectset(1,nobject)(1:4)='MASS'
            do k=5,81
               objectset(1,nobject)(k:k)=' '
            enddo
            if(textpart(2)(1:6).eq.'ELSET=') then
               read(textpart(2)(7:86),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
            endif     
            objectset(3,nobject)(81:81)=' '
            ipos=index(objectset(3,nobject),' ')
            objectset(3,nobject)(ipos:ipos)='E'
         elseif(textpart(1)(1:11).eq.'SHAPEENERGY') then
            nobject=nobject+1
            objectset(1,nobject)(1:11)='SHAPEENERGY'
            do k=12,81
               objectset(1,nobject)(k:k)=' '
            enddo
            if(textpart(2)(1:6).eq.'ELSET=') then
               read(textpart(2)(7:86),'(a80)',iostat=istat) 
     &              objectset(3,nobject)(1:80) 
               objectset(3,nobject)(81:81)=' '
               ipos=index(objectset(3,nobject),' ')
               objectset(3,nobject)(ipos:ipos)='E'
            endif
            nener=1
         elseif(textpart(1)(1:6).eq.'STRESS') then
            nobject=nobject+1
            objectset(1,nobject)(1:6)='STRESS'
            do k=7,81
               objectset(1,nobject)(k:k)=' '
            enddo
            do i=2,n
               if(textpart(i)(1:7).eq.'ENTITY=') then
                  read(textpart(i)(8:87),'(a80)',iostat=istat) 
     &                 objectset(2,nobject)(1:80) 
               elseif(textpart(i)(1:6).eq.'ELSET=') then
                  read(textpart(i)(7:86),'(a80)',iostat=istat) 
     &                 objectset(3,nobject)(1:80) 
                  objectset(3,nobject)(81:81)=' '
                  ipos=index(objectset(3,nobject),' ')
                  objectset(3,nobject)(ipos:ipos)='E'
               endif
            enddo
         else
            write(*,*) '*ERROR reading *OBJECTIVE'
            write(*,*) '       objective function not known'
            call inputerror(inpc,ipoinpc,iline,
     &"*OBJECTIVE%")
         endif
!     
         call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
         if((istat.lt.0).or.(key.eq.1)) exit
!     
      enddo
!     
      return
      end
      
      
