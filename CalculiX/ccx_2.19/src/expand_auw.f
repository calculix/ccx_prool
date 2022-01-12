!     
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2021 Guido Dhondt
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

C     >    *MASSLESS DYNAMIC CONTACT*: Expand contact directions matrix W_b
C     >    from nodal size into DOF size

C     > @param auw       the param
C     > @param jqw       the param
C     > @param iroww     the param
C     > @param nslavs    the param
C     > @param auwnew    the param
C     > @param jqwnew    the param
C     > @param irowwnew  the param
C     > @param neqslav   the param
C     > @param nactdof   the param
C     > @param mi        the param
C     > @param ktot      the param
C     > @param neqtot    the param
C     > @param islavnode the param
C     > @param imastnode the param

C     > @see    massless.c
C     > @author Guido Dhondt, Carlo Monjaraz
C     
      subroutine expand_auw(auw,jqw,iroww,nslavs,auwnew,jqwnew,
     &     irowwnew,neqslav,nactdof,mi,ktot,neqtot,islavnode,imastnode)
!     
      implicit none
!     
      integer nzswnew,jqwnew(*),jqw(*),iroww(*),nslavs,
     &     irowwnew(*),i,neqslav,idof,nodes,idir,j,
     &     length,id,mi(*),nactdof(0:mi(2),*),ktot(*),neqtot,
     &     nodem,islavnode(*),imastnode(*),kflag,
     &     jdofm,idirs,idirm,nodemrel,filedebug
!     
      real*8 auw(*),auwnew(*)
!     
      filedebug=0
!     
      nzswnew=0
!     
!     loop over all columns of auwnew 
!     
      do i=1,neqslav
        jqwnew(i)=nzswnew+1
        nodes=islavnode((i-1)/3 +1 )
!     
!     first node in the MPC is the slave node; identify its row
!     position
!     
        do idirs=1,3
          idof=nactdof(idirs,nodes)
          call nident(ktot,idof,neqtot,id)
          if(id.gt.0) then
            if(ktot(id).eq.idof) then
              nzswnew=nzswnew+1
              auwnew(nzswnew)=auw(jqw(i)-1+idirs)
              irowwnew(nzswnew)=id
            endif
          endif
        enddo
!     
!     independent (master) dofs
!     
        do j=jqw(i)+3,jqw(i+1)-1
          jdofm=iroww(j)-3*nslavs
          nodemrel=(jdofm-1)/3+1 ! position of master node in field imastnode
          nodem=imastnode(nodemrel) ! master node
          idirm=jdofm-3*(nodemrel-1) ! coordinate direction
          idof=nactdof(idirm,nodem)
          call nident(ktot,idof,neqtot,id)
          if(id.gt.0) then
            if(ktot(id).eq.idof) then
              nzswnew=nzswnew+1
              auwnew(nzswnew)=auw(j)
              irowwnew(nzswnew)=id
              cycle
            endif
          endif
        enddo
      enddo
!     
      jqwnew(neqslav+1)=nzswnew+1
!     
!     sorting the column in auwnew
!     
      kflag=2
      do i=1,neqslav
        length=jqwnew(i+1)-jqwnew(i)
        call isortid(irowwnew(jqwnew(i)),auwnew(jqwnew(i)),length,
     &       kflag)
      enddo
!     
!     DEBUG PRINT WB MATRIX
      if(filedebug>0)then
        write(*,*) 'DEBUG expand_auw(): Printing Wb expanded to disk'
        open(unit=31416,file='WbNew_colsort.csv')
        do i=1,neqslav
          idir=i-((i-1)/3)*3
          if(idir==1) then
            write(31416,*) islavnode((i-1)/3 +1 ) , ',N'
          endif
          if(idir==2) then
            write(31416,*) islavnode((i-1)/3 +1 ) , ',T1'
          endif
          if(idir==3) then
            write(31416,*) islavnode((i-1)/3 +1 ) , ',T2'
          endif
        enddo
        close(31416)
!
        open(unit=31416,file='WbNew_size.csv')
        write(31416,*) neqtot
        write(31416,*) neqslav
        close(31416)

        open(unit=31416,file='WbNew_auw.csv')
        do i=1,nzswnew
          write(31416,*) auwnew(i)
        enddo
        close(31416)

        open(unit=31416,file='WbNew_iroww.csv')
        do i=1,nzswnew
          write(31416,*) irowwnew(i)
        enddo
        close(31416)

        open(unit=31416,file='WbNew_jqw.csv')
        do i=1,neqslav+1
          write(31416,*) jqwnew(i)
        enddo
        close(31416)
      endif
!
      return
      end

