! $Id: GAUSSP.f,v 1.1 2009/09/16 14:06:48 bmy Exp $
      SUBROUTINE GAUSSP (N,XPT,XWT)
!-----------------------------------------------------------------------
!  Loads in pre-set Gauss points for 4 angles from 0 to +1 in cos(theta)=mu
!-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER N,I
      REAL*8  XPT(N),XWT(N)
      REAL*8 GPT4(4),GWT4(4)
      DATA GPT4/.06943184420297D0,.33000947820757D0,.66999052179243D0,
     G          .93056815579703D0/
      DATA GWT4/.17392742256873D0,.32607257743127D0,.32607257743127D0,
     W          .17392742256873D0/
      N = 4
      DO I=1,N
        XPT(I) = GPT4(I)
        XWT(I) = GWT4(I)
      ENDDO
      RETURN
      END
