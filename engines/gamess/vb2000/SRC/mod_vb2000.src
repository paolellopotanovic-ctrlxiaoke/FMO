! This file contains modules needed by VB2000 code
! MODULE mod_vb2000
!> @brief   VB2000 common interface to GAMESS 
!> @details Replaces the following common block:
!>          COMMON /VBINTF/ VBENGY,LPROP,NOSYMV,MAXOCC,LNOS,JSTEP,MGMS
!>          It is used in vb2000.src, vb2gms.src,
!>          gamess.src, inputa.src,parley.src, and statpt.src
!
!> @author  David Sousa
!> @date    Sep, 2021
!
!> @params  vbengy - final energy of the VB2000 calculation
!> @params  lprop  - flag that specifies property calculation (0 if $NOVBPROP is present)
!> @params  nosymv - flag that stores value of NOSYM in $CONTRL
!> @params  maxocc - maximum number of occupied orbitals
!> @params  lnos   - flag that toggles natural orbital printing in AIMPAC ($PRINT_NOS)
!> @params  jstep  - number of current iteration (OPTIMIZE,HESSIAN,IRC,etc.)
!> @params  mgms   - flag that skips GAMESS HF calculation ($SKIPGMSHF)
      module mod_vb2000
      implicit none
      double precision :: vbengy
      integer          :: lprop, nosymv, maxocc, lnos, jstep, mgms
      end module mod_vb2000
!
! MODULE genwrk
!> @brief   Dynamic memory allocation in VB2000
!> @details Replaces the following common block:
!>          COMMON /GENWRK/RWRK(MAXWK)
!>          used only in vb2000.src
!
!> @author  David Sousa
!> @date    Sep, 2022
!
!> @params  RWRK - main array of the program
      module genwrk
      implicit none 
      double precision, allocatable :: RWRK(:)
      end module genwrk
!
! MODULE vb2mem
!> @brief   Variables related to dynamic memory allocation in VB2000
!> @details Replaces the following common block:
!>          COMMON /VB2MEM/LOCX,LOCMEM,LOFFS,MEMWRK,MSCRATCH,MHIGH
!>          used only in vb2000.src
!
!> @author  David Sousa
!> @date    Sep, 2022
      module vb2mem
      implicit none 
      integer :: memwrk,mscratch,mhigh
      end module vb2mem
