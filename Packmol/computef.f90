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
!
! Subroutine that computes the function value
!
  
subroutine computef(n,x,f) 
      
  use sizes
  use compute_data
  implicit none

  integer :: n, i, j, k
  integer :: ilugan, ilubar, icart, itype, imol, iatom, idatom, &
             iboxx, iboxy, iboxz 

  double precision :: v1(3), v2(3), v3(3) 
  double precision :: x(n)
  double precision :: f,fparc,fplus
  double precision :: xtemp, ytemp, ztemp
  double precision :: xbar, ybar, zbar
  double precision :: beta, gama, teta
  double precision :: flast

  ! Reset function value

  f = 0.d0 
  frest = 0.d0
  fdist = 0.d0

  ! Reset boxes

  if(.not.init1) call resetboxes()

  ! Transform baricenter and angles into cartesian coordinates 
  ! Computes cartesian coordinates from vector x and coor 

  ilubar = 0 
  ilugan = ntotmol*3 
  icart = natfix

  do itype = 1, ntype 
    if(.not.comptype(itype)) then
      icart = icart + nmols(itype)*natoms(itype)
    else
    do imol = 1, nmols(itype) 

      xbar = x(ilubar+1) 
      ybar = x(ilubar+2) 
      zbar = x(ilubar+3) 
  
      ! Computing the rotation matrix

      beta = x(ilugan+1)
      gama = x(ilugan+2)
      teta = x(ilugan+3)

      call eulerrmat(beta,gama,teta,v1,v2,v3)  

      ! Looping over the atoms of this molecule
  
      idatom = idfirst(itype) - 1
      do iatom = 1, natoms(itype) 

        icart = icart + 1
        idatom = idatom + 1

        ! Computing the cartesian coordinates for this atom

        call compcart(icart,xbar,ybar,zbar, &
                      coor(idatom,1),coor(idatom,2),coor(idatom,3), &
                      v1,v2,v3)

        ! Adding to f the value relative to constraints for this atom

        call comprest(icart,fplus)
        f = f + fplus
        frest = dmax1(frest,fplus)
        if(move) fatom(icart) = fatom(icart) + fplus

        ! Putting atoms in their boxes

        if(.not.init1) then

          xtemp = xcart(icart,1) - sizemin(1) 
          ytemp = xcart(icart,2) - sizemin(2) 
          ztemp = xcart(icart,3) - sizemin(3) 

          iboxx = int(xtemp/boxl(1)) + 1
          iboxy = int(ytemp/boxl(2)) + 1
          iboxz = int(ztemp/boxl(3)) + 1

          if(xtemp.le.0) iboxx = 1
          if(ytemp.le.0) iboxy = 1
          if(ztemp.le.0) iboxz = 1 
          if(iboxx.gt.nboxes(1)) iboxx = nboxes(1)
          if(iboxy.gt.nboxes(2)) iboxy = nboxes(2)
          if(iboxz.gt.nboxes(3)) iboxz = nboxes(3)

          latomnext(icart) = latomfirst(iboxx,iboxy,iboxz)
          latomfirst(iboxx,iboxy,iboxz) = icart

          ibtype(icart) = itype
          ibmol(icart) = imol

        end if

      end do 
 
      ilugan = ilugan + 3 
      ilubar = ilubar + 3 

    end do
    end if
  end do            

  if(init1) return

  ! Minimum distance function evaluation

  do i = 1, nboxes(1)
    do j = 1, nboxes(2)
      do k = 1, nboxes(3)

        icart = latomfirst(i,j,k)
        do while ( icart .ne. 0 ) 

          if(comptype(ibtype(icart))) then

            ! Vector that keeps the value for this atom

            if(move) flast = f 

            ! Interactions inside box

            f = f + fparc(icart,latomnext(icart))

            ! Interactions of boxes that share faces

            f = f + fparc(icart,latomfirst(i+1,j,k))
            f = f + fparc(icart,latomfirst(i,j+1,k))
            f = f + fparc(icart,latomfirst(i,j,k+1))

            ! Interactions of boxes that share axes

            f = f + fparc(icart,latomfirst(i+1,j+1,k))
            f = f + fparc(icart,latomfirst(i+1,j,k+1))
            f = f + fparc(icart,latomfirst(i+1,j-1,k))
            f = f + fparc(icart,latomfirst(i+1,j,k-1))
            f = f + fparc(icart,latomfirst(i,j+1,k+1))
            f = f + fparc(icart,latomfirst(i,j+1,k-1))

            ! Interactions of boxes that share vertices

            f = f + fparc(icart,latomfirst(i+1,j+1,k+1))
            f = f + fparc(icart,latomfirst(i+1,j+1,k-1))
            f = f + fparc(icart,latomfirst(i+1,j-1,k+1))
            f = f + fparc(icart,latomfirst(i+1,j-1,k-1))

            ! If going to move bad molecules, update fatom

            if(move) fatom(icart) = fatom(icart) + f - flast

          end if

          icart = latomnext(icart)
        end do
  
      end do
    end do
  end do

  return
end subroutine computef

