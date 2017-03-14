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
! Module that carries the input parameters read from the input file
!

module input

  use sizes
  implicit none

  integer :: nlines
  integer :: nrest
  integer :: seed
  integer :: nloop
  integer :: writeout
  integer :: ntfix
  integer :: ntcon(9) 

  integer, allocatable :: nconnect(:,:)  ! (ntotat,9)
  integer, allocatable :: irestline(:) ! (maxrest)
  integer, allocatable :: linestrut(:,:) ! (ntype,2)
  integer, allocatable :: resnumbers(:) ! (ntype)
  integer, allocatable :: maxcon(:) ! (ntotat)
  integer, allocatable :: input_itype(:) ! (ntype)

  double precision :: dism
  double precision :: precison
  double precision :: sidemax
  double precision :: discale
  double precision :: movefrac
  double precision :: add_sides_fix    
  double precision :: precision
  double precision :: fbins

  double precision, allocatable :: amass(:) ! (ntotat)
  double precision, allocatable :: charge(:) ! (ntotat)
  
  logical :: writebad
  logical :: tinker
  logical :: pdb
  logical :: xyz
  logical :: moldy
  logical :: check
  logical :: chkgrad
  logical :: randini
  logical :: movebadrandom
  logical :: add_amber_ter
  logical :: add_box_sides
  logical :: fix
  logical :: avoidoverlap

  logical, allocatable :: changechains(:) ! (ntype)
  logical, allocatable :: fixedoninput(:) ! (ntype)
  
  character(len=200) :: xyzout

  character(len=3), allocatable :: ele(:) ! (ntotat)
  character(len=80), allocatable :: pdbfile(:) ! (ntype)
  character(len=200), allocatable :: name(:) ! (ntype)
  character(len=200), allocatable :: keyword(:,:) ! (nlines,maxkeywords)
  character(len=200), allocatable :: inputfile(:) ! (nlines)
  character(len=200), allocatable :: restart_from(:), restart_to(:) ! (0:ntype)

end module input
