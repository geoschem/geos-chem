!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: rd_js
!
! !DESCRIPTION: Rereads the "ratj.d" file to map photolysis rate to reaction
!  Read in quantum yield 'jfacta' and fastj labels 'jlabel'
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE RD_JS( NJ1, NAMFIL, am_I_Root )
!
! !USES:
!
      USE CMN_FJ_MOD, ONLY : JPMAX,  JPPJ
      USE CMN_FJ_MOD, ONLY : RNAMES, BRANCH
      USE JV_CMN_MOD

      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      INTEGER,          INTENT(INOUT) :: nj1         ! Logical unit #
      CHARACTER(LEN=*), INTENT(IN)    :: namfil      ! File name
      LOGICAL,          INTENT(IN)    :: am_I_Root   ! Is this the root CPU? 
!
! !INPUT PARAMETERS: 
!
!
! !REMARKS:
!     jfacta    Quantum yield (or multiplication factor) for photolysis
!     jlabel    Reference label identifying appropriate J-value to use
!     ipr       Photolysis reaction counter - should total 'jppj'
! 
! !REVISION HISTORY: 
!  01 Jun 1998 - P. Murti    - Initial version
!  30 Jul 2012 - R. Yantosca - Now accept am_I_Root as an argument when
!                              running with the traditional driver main.F
!  30 Jul 2012 - R. Yantosca - Add reference to findFreeLUN under DEVEL tag
!  03 Aug 2012 - R. Yantosca - Restore NJ1 to INTENT(IN) status
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
      integer       ipr, i
      character*120 cline
c
c Reread the ratj.d file to map photolysis rate to reaction
c                     Read in quantum yield jfacta and fastj label jlabel
      ipr=0
      ! open file
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
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Additional code to read reaction names and branch numbers
!  (ppm, 6/98, bmy, 9/99)     
!         read (cline(7:10),"(a4)") rnames(ipr)
! FP increase length of name (6/2009)
         read (cline(7:13),"(a7)") rnames(ipr)
         rnames(ipr) = trim(rnames(ipr))
         branch(ipr) = 1
         do i=1,ipr-1
            if (rnames(ipr) == rnames(i)) branch(ipr) = branch(i) + 1
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         go to 10
      endif
 20   close(nj1)
      if(ipr.ne.jppj) then
         IF ( am_I_Root ) write(6,1000) ipr,jppj
         stop
      endif
c
c Print details to standard output
      IF ( am_I_Root ) THEN
         write(6,1100) ipr
         write(6,1200) (i, jlabel(i), jfacta(i),i=1,ipr)
      ENDIF
c
      return
 1000 format(' Error: ',i3,' photolysis labels but ',i3,' reactions')
 1100 format(' Fast-J Photolysis Scheme: considering ',i2,' reactions')
 1200 format(3x,10(3(i2,': ',a7,' (Q.Y. ',f5.3,') '),/,3x))
      END SUBROUTINE RD_JS
!EOC
