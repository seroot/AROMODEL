!  
!  Written by Leandro Martínez, 2009-2011.
!  Copyright (c) 2009-2011, Leandro Martínez, Jose Mario Martinez,
!  Ernesto G. Birgin.
!  
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU General Public License
!  as published by the Free Software Foundation; either version 2
!  of the License, or (at your option) any later version.
!  
! subroutine restmol: either compute the restraint function
!                     value for a single molecule or solve
!                     the problem of puting this molecule
!                     in the restraint region
!  

subroutine restmol(itype,ilubar,n,x,fx,solve)

  use sizes
  use compute_data
  use usegencan
  implicit none

  integer :: n, nsafe, ntotsafe, itype, i, ilubar, nmoltype, ip1, ip2
  double precision :: x(n), fx
  logical :: solve, initsafe
 
  ! Saving global problem variables

  nsafe = n
  ntotsafe = ntotmol
  nmoltype = nmols(itype)
  do i = 1, ntype
    compsafe(i) = comptype(i)
  end do
  initsafe = init1

  ! Preparing system to solve for this molecule

  n = 6
  ntotmol = 1      
  nmols(itype) = 1
  xmol(1) = x(ilubar+1)
  xmol(2) = x(ilubar+2)
  xmol(3) = x(ilubar+3)
  xmol(4) = x(ilubar+ntotsafe*3+1)
  xmol(5) = x(ilubar+ntotsafe*3+2)
  xmol(6) = x(ilubar+ntotsafe*3+3)
  do i = 1, ntype
    if(i.eq.itype) then
      comptype(i) = .true.
    else
      comptype(i) = .false.
    end if
  end do
  init1 = .true.
      
  ! If not going to solve the problem, compute energy and return

  if(.not.solve) then
    call computef(n,xmol,fx)
    ! Otherwise, put this molecule in its constraints
  else
    ip1 = iprint1
    ip2 = iprint2
    iprint1 = 0
    iprint2 = 0
    call pgencan(n,xmol,fx)
    iprint1 = ip1
    iprint2 = ip2
  end if       

  ! Restoring original problem data

  ntotmol = ntotsafe
  n = nsafe
  nmols(itype) = nmoltype
  x(ilubar+1) = xmol(1) 
  x(ilubar+2) = xmol(2) 
  x(ilubar+3) = xmol(3) 
  x(ilubar+ntotmol*3+1) = xmol(4) 
  x(ilubar+ntotmol*3+2) = xmol(5) 
  x(ilubar+ntotmol*3+3) = xmol(6) 
  do i = 1, ntype
    comptype(i) = compsafe(i)
  end do
  init1 = initsafe

  return
end subroutine restmol

