C $Id: rd_js.f,v 1.1 2009/09/16 14:06:14 bmy Exp $
      subroutine rd_js(nj1,namfil)
C-----------------------------------------------------------------------
c  Reread the ratj.d file to map photolysis rate to reaction
c  Read in quantum yield 'jfacta' and fastj label 'jlabel'
C-----------------------------------------------------------------------
c
c     jfacta    Quantum yield (or multiplication factor) for photolysis
c     jlabel    Reference label identifying appropriate J-value to use
c     ipr       Photolysis reaction counter - should total 'jppj'
c
C-----------------------------------------------------------------------
      IMPLICIT NONE

#     include "cmn_fj.h"
#     include "jv_cmn.h"
c
      integer nj1, ipr, i
      character*6  namfil
      character*120 cline
c
c Reread the ratj.d file to map photolysis rate to reaction
c                     Read in quantum yield jfacta and fastj label jlabel
      ipr=0
      open(nj1,file=namfil,status='old',form='formatted')
 10   read(nj1,'(a)',err=20) cline
      if(cline(2:5).eq.'9999') then
         go to 20
      elseif(cline(1:1).eq.'#') then
         go to 10
      elseif(cline(5:5).eq.'$') then
         go to 10
      else
         ipr=ipr+1
         read(cline(79:83),'(f5.1)') jfacta(ipr)
         read(cline(86:92),'(a7)')   jlabel(ipr)
         jfacta(ipr)=jfacta(ipr)/100.d0
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C  Additional code to read reaction names and branch numbers
C  (ppm, 6/98, bmy, 9/99)     
         read (cline(7:10),"(a4)") rnames(ipr)
         rnames(ipr) = trim(rnames(ipr))
         branch(ipr) = 1
         do i=1,ipr-1
            if (rnames(ipr) == rnames(i)) branch(ipr) = branch(i) + 1
         enddo
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         go to 10
      endif
 20   close(nj1)
      if(ipr.ne.jppj) then
         write(6,1000) ipr,jppj
         stop
      endif
c
c Print details to standard output
      write(6,1100) ipr
      write(6,1200) (i, jlabel(i), jfacta(i),i=1,ipr)
c
      return
 1000 format(' Error: ',i3,' photolysis labels but ',i3,' reactions')
 1100 format(' Fast-J Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,10(3(i2,': ',a7,' (Q.Y. ',f5.3,') '),/,3x))
      end
