       subroutine write49
     x     (namediag, xtrac, xtau1, xtau2, xskip,
     x      ni, nj, nl, xi0, xj0, xl0,
     x      xarray,iut)

       
c write the tsYYMMDD.bpch binary file (diagnostic 49)

       implicit none

       character*40  namediag
       integer xtrac, xskip
       integer ni, nj, nl, xi0, xj0, xl0, iut
       integer i, j, l
       real*4 xarray(ni, nj, nl)
       real*8 xtau1, xtau2

       write(iut) namediag, xtrac, xtau1, xtau2, xskip
       write(iut) ni, nj, nl, xi0, xj0, xl0
       write(iut) (((xarray(i,j,l), i=1,ni), j=1,nj), l=1,nl)


!       write(74,100) namediag, xtau1, xtau2, xskip
!       write(74,101) ni, nj, nl, xi0, xj0, xl0
!       write(74,102) (((xarray(i,j,l), i=1,ni), j=1,nj), l=1,nl)
! 100   format(a40,f7.1,f7.1,i8)
! 101   format(6i5)
! 102   format(5(1pe11.2))
       return
       end
