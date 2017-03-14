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
! Subroutine that sets the sizes of all allocatable arrays
!

subroutine setsizes()

  use sizes
  use compute_data
  use input
  use usegencan
  use flashsort
 
  implicit none
  integer :: i, ival, ilast, iline, itype
  integer :: ioerr
  integer :: strlength
  character(len=200) :: record, word, blank

  ! Instructions on how to run packmol

  write(*,*) ' Packmol must be run with: packmol < inputfile.inp '
  write(*,*)
  write(*,*) ' Userguide at: www.ime.unicamp.br/~martinez/packmol '
  write(*,*)
      
  ! Getting input lines from the input file

  write(*,*) ' Reading input file... (Control-C aborts)'   

  do i = 1, 200
    blank(i:i) = ' '
  end do
  nlines = 0
  maxkeywords = 0
  ntype = 0
  do
    read(5,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit

    ! Remove comments

    i = 0
    do while( i < 200 ) 
      i = i + 1
      if ( record(i:i) == '#' ) exit
    end do
    i = i - 1
    if ( i > 0 ) then
      record = record(1:i)//blank(i+1:200)
    else
      cycle
    end if
    if ( strlength(record) < 1 ) cycle
      
    ! Number of lines of the input file

    nlines = nlines + 1

    ! Check the number of keywords in this line

    i = 0
    ival = 0
    do while(i < 200)
      i = i + 1
      ilast = i
      do while(record(i:i) > ' ' .and. i < 200)
        i = i + 1
      end do
      if(i > ilast) then
        ival = ival + 1
        maxkeywords = max(maxkeywords,ival)
      end if
    end do  
  end do
  rewind(5)

  allocate(inputfile(nlines),keyword(nlines,maxkeywords))

  ! Read input to inputfile array

  iline = 0
  do
    read(5,"(a200)",iostat=ioerr) record
    if ( ioerr /= 0 ) exit

    ! Remove comments

    i = 0
    do while( i < 200 ) 
      i = i + 1
      if ( record(i:i) == '#' ) exit
    end do
    i = i - 1
    if ( i > 0 ) then
      record = record(1:i)//blank(i+1:200)
    else
      cycle
    end if
    if ( strlength(record) < 1 ) cycle

    iline = iline + 1
    inputfile(iline) = record
  end do

  ! Read all keywods into keyword array

  call getkeywords()

  ! Checking the filetype of coordinate files (default is pdb)

  tinker = .false.
  pdb = .false.
  xyz = .false.
  moldy = .false.
  fbins = dsqrt(3.d0)
  do i = 1, nlines
    if(keyword(i,1).eq.'filetype') then
      if(keyword(i,2).eq.'tinker') tinker = .true.
      if(keyword(i,2).eq.'pdb') pdb = .true.
      if(keyword(i,2).eq.'xyz') xyz = .true.
      if(keyword(i,2).eq.'moldy') moldy = .true.
    end if
    if(keyword(i,1).eq.'fbins') then
      record = keyword(i,2)
      read(record,*,iostat=ioerr) fbins
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Invalid value for fbins. '
        stop
      end if
    end if
  end do
  if(.not.pdb.and..not.tinker.and..not.xyz.and..not.moldy) then
    pdb = .true.
    write(*,*)
    write(*,*)' WARNING: File type not (correctly?) specified, using PDB'
  end if  

  ! Getting the number of different types of molecules

  ntype = 0
  do iline = 1, nlines
    if ( keyword(iline,1) == "structure" ) then
      ntype = ntype + 1
      if ( keyword(iline,2) == "none" ) then
        write(*,*) ' ERROR: structure without filename. '
        write(*,*) ' The syntax must be, for example: structure water.pdb '
        stop 
      end if
    end if
  end do

  allocate(nmols(ntype),natoms(ntype),idfirst(ntype),constrain_rot(ntype,3),&
           rot_bound(ntype,3,2),dmax(ntype),&
           cmxmin(ntype),cmymin(ntype),cmzmin(ntype),&
           cmxmax(ntype),cmymax(ntype),cmzmax(ntype),&
           comptype(ntype),compsafe(ntype),&
           restart_from(0:ntype),restart_to(0:ntype))

  ! Reading the number of molecules of each type, and the number of atoms
  ! of each molecule type

  itype = 0
  do iline = 1, nlines
    if ( keyword(iline,1) == "structure" ) then
      itype = itype + 1
      natoms(itype) = 0
      nmols(itype) = 0

      ! Read the number of atoms of this type of molecule

      open(10,file=keyword(iline,2),status='old',iostat=ioerr)
      if( ioerr /= 0 ) call failopen(keyword(iline,2))
      if ( pdb ) then
        do
          read(10,"(a200)",iostat=ioerr) record
          if ( ioerr /= 0 ) exit
          if ( record(1:4) == "ATOM" .or. record(1:6) == "HETATM" ) then
            natoms(itype) = natoms(itype) + 1
          end if
        end do
      end if
      if ( tinker ) then
        do
          read(10,*,iostat=ioerr) i
          if ( ioerr /= 0 ) cycle
          natoms(itype) = i
          exit
        end do
      end if
      if ( xyz ) then
        read(10,*,iostat=ioerr) i
        if ( ioerr == 0 ) natoms(itype) = i
      end if
      if ( moldy ) then
        read(10,*,iostat=ioerr) word, i
        if ( ioerr == 0 ) natoms(itype) = i
      end if
      close(10)
      if ( natoms(itype) == 0 ) then
        write(*,*) ' ERROR: Could not read any atom from file: ', &
                   trim(adjustl(keyword(iline,2)))
      end if

    end if

    ! Read number of molecules for each type 

    if ( keyword(iline,1) == "number" ) then
      read(keyword(iline,2),*,iostat=ioerr) nmols(itype)
      if ( ioerr /= 0 ) then
        write(*,*) ' ERROR: Error reading number of molecules of type ', itype
        stop  
      end if
    end if

  end do
  do itype = 1, ntype
    if ( nmols(itype) == 0 ) then  
      write(*,*) ' Warning: Number of molecules not set for type '&
                 ,itype,': assuming 1 '
      nmols(itype) = 1 
    end if
  end do

  ! Total number of atoms and molecules

  ntotat = 0
  ntotmol = 0
  do itype = 1, ntype
    ntotat = ntotat + nmols(itype)*natoms(itype)
    ntotmol = ntotmol + nmols(itype)
  end do

  ! The number of variables of the problem

  nn = ntotmol*6

  ! The number of bins of the linked cell method in each direction

  nbp = int((fbins*dble(ntotat))**(1.d0/3.d0)) + 1

  ! Checking the total number of restrictions defined

  i = 0
  do iline = 1, nlines
    if ( keyword(iline,1) == 'fixed' .or. &
         keyword(iline,1) == 'inside' .or. &
         keyword(iline,1) == 'outside' .or. &
         keyword(iline,1) == 'over' .or. &
         keyword(iline,1) == 'below' .or. &
         keyword(iline,1) == 'constrain_rotation' ) then
      i = i + 1 
    end if
  end do
  maxrest = i
  mrperatom = i

  ! Allocate arrays depending on ntotat, nn, nbp, maxrest and mrperatom

  allocate(nratom(ntotat),iratom(ntotat,mrperatom),ibmol(ntotat),&
           ibtype(ntotat),xcart(ntotat,3),coor(ntotat,3),&
           radius(ntotat),radius_ini(ntotat),&
           gxcar(ntotat,3),&
           fatom(ntotat),latomnext(ntotat),&
           fmol(ntotat),radiuswork(ntotat))
  allocate(ityperest(maxrest),restpars(maxrest,9))
  allocate(latomfirst(0:nbp+1,0:nbp+1,0:nbp+1),&
           latomfix(0:nbp+1,0:nbp+1,0:nbp+1))
  allocate(xmol(nn))

  ! Allocate other arrays used for input and output data

  allocate(nconnect(ntotat,9),maxcon(ntotat),&
           amass(ntotat),charge(ntotat),ele(ntotat))

  allocate(irestline(maxrest),linestrut(ntype,2),resnumbers(ntype),&
           input_itype(ntype),changechains(ntype),&
           fixedoninput(ntype),pdbfile(ntype),name(ntype))

  ! Allocate vectors for flashsort

  allocate(indflash(ntotat),lflash(ntotat))

  ! Allocate arrays for GENCAN

  allocate(l(nn),u(nn),wd(8*nn),wi(nn),g(nn))

end subroutine setsizes

