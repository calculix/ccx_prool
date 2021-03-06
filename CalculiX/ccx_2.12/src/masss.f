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
      subroutine masss(inpc,textpart,nrhcon,nmat,ntmat_,
     &        rhcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,nmat_,set,istartset,iendset,ialset,
     &        nset,ielmat,ielorien,ipoinpc,mi)
!
!     reading the input deck: *MASS
!
      implicit none
!
      character*1 inpc(*)
      character*80 matname(*)
      character*81 set(*),elset
      character*132 textpart(16)
!
      integer mi(*),nrhcon(*),nmat,ntmat_,ntmat,istep,
     &  n,key,i,istat,istartset(*),
     &  iendset(*),irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*),nmat_,
     &  ialset(*),ipos,nset,j,k,ielmat(mi(3),*),ielorien(mi(3),*),
     &  ipoinpc(0:*)  
!
      real*8 rhcon(0:1,ntmat_,*)
!
      ntmat=0
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR reading *MASS: *MASS should be placed'
         write(*,*) '  before all step definitions'
         call exit(201)
      endif
!
      nmat=nmat+1
      if(nmat.gt.nmat_) then
         write(*,*) '*ERROR reading *MASS: increase nmat_'
         call exit(201)
      endif
      matname(nmat)(1:4)='MASS'
      do i=5,80
         matname(nmat)(i:i)=' '
      enddo
!
      do i=2,n
         if(textpart(i)(1:6).eq.'ELSET=') then
            elset=textpart(i)(7:86)
            elset(81:81)=' '
            ipos=index(elset,' ')
            elset(ipos:ipos)='E'
         else
            write(*,*) 
     &           '*WARNING reading *MASS: parameter not recognized:'
            write(*,*) '         ',
     &           textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*MASS%")
         endif
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*) '*ERROR reading *MASS: definition of the'
         write(*,*) '       mass is not complete'
         call inputerror(inpc,ipoinpc,iline,
     &"*MASS%")
         call exit(201)
      endif
!
      ntmat=ntmat+1
      nrhcon(nmat)=ntmat
      if(ntmat.gt.ntmat_) then
         write(*,*) '*ERROR reading *MASS: increase ntmat_'
         call exit(201)
      endif
      read(textpart(1)(1:20),'(f20.0)',iostat=istat)
     &     rhcon(1,ntmat,nmat)
      if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*MASS%")
      rhcon(0,ntmat,nmat)=0.d0
!     
      if(ntmat.eq.0) then
         write(*,*) '*ERROR reading *MASS: *MASS card without data'
         call exit(201)
      endif
      do i=1,nset
         if(set(i).eq.elset) exit
      enddo
      if(i.gt.nset) then
         elset(ipos:ipos)=' '
         write(*,*) '*ERROR reading *MASS: element set ',elset
         write(*,*) '       has not yet been defined. '
         call inputerror(inpc,ipoinpc,iline,
     &"*MASS%")
         call exit(201)
      endif
!     
!     assigning the elements of the set the appropriate material
!     
      do j=istartset(i),iendset(i)
         if(ialset(j).gt.0) then
            ielmat(1,ialset(j))=nmat
            ielorien(1,ialset(j))=0
         else
            k=ialset(j-2)
            do
               k=k-ialset(j)
               if(k.ge.ialset(j-1)) exit
               ielmat(1,k)=nmat
               ielorien(1,k)=0
            enddo
         endif
      enddo
!     
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!     
      return
      end

