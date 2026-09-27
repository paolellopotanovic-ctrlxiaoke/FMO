!> @brief  The old common block
!>       COMMON /IOFILE/ ir,iw,ip,is,ipk,idaf,nav,ioda(950)
!>
module comm_IOFILE
  IMPLICIT NONE
  INTEGER :: IR         ! input file
  INTEGER :: IW         ! output file (normal log file)
  INTEGER :: IP         ! the supplemental output file, file
  ! ~/scr/xxx.dat perhaps
  INTEGER :: IJK        ! repulsion integrals in the AO basis
  INTEGER :: IJKT       ! IJKT is the repulsion integrals,
  ! transformed (the T) to the MO basis.
  INTEGER :: IDAF       ! unit number for the direct access
  ! file (DAREAD/DAWRIT).
  INTEGER :: NAV        ! obsolete (associated variable for Fortran IO)
  INTEGER, DIMENSION(950) :: IODA ! a pointer to the physical
  ! record number where the "logical record" is stored.

  ! For extracting values from common block
  COMMON /IOFILE/ IR, IW, IP, IJK, IJKT, IDAF, NAV, IODA
end module comm_IOFILE

!> @brief  The old common block
!>       COMMON /PAR   / me,master,nproc,ibtyp,iptim,goparr,dskwrk,maswrk
module comm_PAR
  IMPLICIT NONE
  INTEGER :: ME             ! rank of the calling process in the DDI group
  INTEGER :: MASTER         ! rank 0, the rank that is supposed to do
  ! output, eg to .log and .dat.
  INTEGER :: NPROC          ! number of compute processes in the DDI group
  INTEGER :: IBTYP          ! 1 = dynamic load balancing, 0 = static load balancing
  INTEGER :: IPTIM          !
  LOGICAL :: GOPARR         ! flag to denote whether or not it's a parallel run
  ! = NPROC.GT.1   or if   PARALL=.TRUE. in $SYSTEM
  LOGICAL :: DSKWRK         ! flag to control how files are read
  ! disk_worker = .true.:
  ! A really big file like IJK (AOINTS) is chopped into subfiles, with
  ! each rank doing part of the AO integrals, and saving its own part
  ! into a separate file.  This is disk_worker being .TRUE. and leads
  ! to files
  !   xxxx.F08       (since IJK=8)
  !   xxxx.F08.001
  !   xxxx.F08.002
  !   xxxx.F08.003
  ! for a four rank job.  Note that 000 is never appended to a file name,
  ! but if the four CPUs are inside the same node, and sharing the same
  ! disk, we have to add the rank number to get a unique file name.
  !
  ! disk_worker = .false.:
  ! Other files such are smaller, and kept on only one rank.
  ! These files are stored only on the MASTER rank, as just
  !   xxx.F25
  ! and shared with other ranks later by DDI_BCAST.
  !
  ! A call to SQREAD or SQWRIT therefore considers disk_worker, in deciding
  ! how many ranks do the I/O and if there's any broadcasting.
  ! DAREAD always does a broadcast to the other ranks.
  LOGICAL :: MASWRK         ! = ME.EQ.MASTER (true if the calling rank
  ! is the master rank)

  ! For extracting values from common block
  COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK, MASWRK
end module comm_PAR

!> @brief  The old common block
!>       COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN, NAUXSHL
module comm_DFTPAR
  use prec, only: fp
  implicit none
  real(kind=fp) :: DFTTYP(20)
  !      DFTTYP(1)=0 no DFT
  !      DFTTYP(1)=1 local DFT
  !      DFTTYP(1)=2 gradient corrected DFT (e.g. GGA)
  !      DFTTYP(2-19)= fractional weights of grid-free functionals
  !      2  X-alpha with a=2/3
  !      3  Hartree-Fock exact exchange scale
  !      4  Becke88 exchange
  !      5  De Pristo-Kress exchange
  !      6  CAMA exchange
  !      7  CAMB exchange (and correlation)
  !      8  WIGNER exchange-correlation
  !      9  WIGNER scaled form
  !     10  WIGNER exponential form
  !     11  VWN5 local correlation
  !     12  Perdew-Wang local correlation
  !     13  LYP gradient corrected correlation (not spin polarized)
  !     14  non-local MP2 correlation scale
  !     15  remaining are unused
  !     16
  !     17
  !     18
  !     19
  !     20  type of auxiliary basis set used
  real(kind=fp) :: EXENA ! alpha DFT exchange + correlation energy
  real(kind=fp) :: EXENB ! beta DFT exchange + correlation energy
  real(kind=fp) :: EXENC ! ???
  integer :: IDFT34      ! Parameter of Grid-free DFT, controlled by THREE key in $DFT group
  !                           if it is set to 3, will use 3 center integrals (RI)
  !                           if it is set to 4, will use 4 center integrals
  integer :: NAUXFUN     ! number of auxiliary functions
  integer :: NAUXSHL     ! number of auxiliary shells

  ! For extracting values from common block
  COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN, NAUXSHL
!$omp threadprivate(/DFTPAR/)
end module comm_DFTPAR
!*MODULE MCPDAT
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:35:43 UTC
!> @params COEBK
!> @params OEMP
!> @params CEMP
!> @params ALPHA
!> @params COEFA
!> @params ALPHB
!> @params COEFB
!> @params AA
!> @params BB
!> @params AGI
!> @params AGJ
!> @params NTERMA
!> @params NTERMB
!> @params NSIT
!> @params NPOS
!> @params NEFF
!> @params NQNDMP
!> @params IPOS
!> @params MPOS
!> @params NTC
!> @params NGC
!> @params IPSET
!> @params IOSET
!> @params NCO
!> @params NATM
module comm_MCPDAT
  use mx_limits, only: MXELM, MXNT, MXPOS, MXIO, MXT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(MXIO) :: COEBK
  DOUBLE PRECISION, DIMENSION(MXNT,MXIO) :: OEMP
  DOUBLE PRECISION, DIMENSION(MXNT,MXIO) :: CEMP
  DOUBLE PRECISION, DIMENSION(MXT,MXELM) :: ALPHA
  DOUBLE PRECISION, DIMENSION(MXT,MXELM) :: COEFA
  DOUBLE PRECISION, DIMENSION(MXT,MXELM) :: ALPHB
  DOUBLE PRECISION, DIMENSION(MXT,MXELM) :: COEFB
  DOUBLE PRECISION, DIMENSION(3) :: AA
  DOUBLE PRECISION, DIMENSION(3) :: BB
  DOUBLE PRECISION :: AGI
  DOUBLE PRECISION :: AGJ
  INTEGER, DIMENSION(MXELM) :: NTERMA
  INTEGER, DIMENSION(MXELM) :: NTERMB
  INTEGER, DIMENSION(MXPOS,MXELM) :: NSIT
  INTEGER, DIMENSION(MXELM) :: NPOS
  INTEGER, DIMENSION(MXELM) :: NEFF
  INTEGER, DIMENSION(4,MXIO) :: NQNDMP
  INTEGER, DIMENSION(MXPOS,MXIO) :: IPOS
  INTEGER, DIMENSION(MXIO) :: MPOS
  INTEGER, DIMENSION(MXIO) :: NTC
  INTEGER, DIMENSION(MXIO) :: NGC
  INTEGER :: IPSET
  INTEGER :: IOSET
  INTEGER :: NCO
  INTEGER :: NATM
end module comm_MCPDAT
!*MODULE MCPITM
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:36:10 UTC
!> @params ZETA
!> @params XI
!> @params AB2
!> @params CB2
!> @params PC2
!> @params XIC
!> @params AD2
!> @params QC2
!> @params ZEFF
!> @params EXA
!> @params EXB
!> @params COEFA
!> @params COEFB
!> @params ETA
!> @params ZETACABC
!> @params ETAD
!> @params ZETADABC
!> @params ZETAD
!> @params XYZCB
!> @params XYZQA
!> @params XYZQD
!> @params XYZQC
!> @params XYZGA
!> @params XYZGB
!> @params SSSC
!> @params SSSD
!> @params IELEM
!> @params IMCP
!> @params IATOM
!> @params NTERMA
!> @params NTERMB
!> @params INTTYP
!> @params MMAX
module comm_MCPITM
  use mx_limits, only: MXM, MXT
  IMPLICIT NONE
  DOUBLE PRECISION :: ZETA
  DOUBLE PRECISION :: XI
  DOUBLE PRECISION :: AB2
  DOUBLE PRECISION :: CB2
  DOUBLE PRECISION :: PC2
  DOUBLE PRECISION, DIMENSION(MXT) :: XIC
  DOUBLE PRECISION, DIMENSION(MXT) :: AD2
  DOUBLE PRECISION, DIMENSION(MXT) :: QC2
  DOUBLE PRECISION :: ZEFF
  DOUBLE PRECISION, DIMENSION(MXT) :: EXA
  DOUBLE PRECISION, DIMENSION(MXT) :: EXB
  DOUBLE PRECISION, DIMENSION(MXT) :: COEFA
  DOUBLE PRECISION, DIMENSION(MXT) :: COEFB
  DOUBLE PRECISION, DIMENSION(MXT) :: ETA
  DOUBLE PRECISION, DIMENSION(MXT) :: ZETACABC
  DOUBLE PRECISION, DIMENSION(MXT) :: ETAD
  DOUBLE PRECISION, DIMENSION(MXT) :: ZETADABC
  DOUBLE PRECISION, DIMENSION(MXT) :: ZETAD
  DOUBLE PRECISION, DIMENSION(3) :: XYZCB
  DOUBLE PRECISION, DIMENSION(MXT,3) :: XYZQA
  DOUBLE PRECISION, DIMENSION(MXT,3) :: XYZQD
  DOUBLE PRECISION, DIMENSION(MXT,3) :: XYZQC
  DOUBLE PRECISION, DIMENSION(MXT,3) :: XYZGA
  DOUBLE PRECISION, DIMENSION(MXT,3) :: XYZGB
  DOUBLE PRECISION, DIMENSION(0:MXM,MXT) :: SSSC
  DOUBLE PRECISION, DIMENSION(MXT) :: SSSD
  INTEGER :: IELEM
  INTEGER :: IMCP
  INTEGER :: IATOM
  INTEGER :: NTERMA
  INTEGER :: NTERMB
  INTEGER :: INTTYP
  INTEGER :: MMAX
end module comm_MCPITM
!*MODULE MCPITMG
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:36:40 UTC
!> @params XYZGC
module comm_MCPITMG
  use mx_limits, only: MXT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(MXT,3) :: XYZGC
end module comm_MCPITMG
!*MODULE MCPOUT
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:37:10 UTC
!> @params GMCP0C
!> @params GMCP0D
module comm_MCPOUT
  use mx_limits, only: MXC
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(MXC,MXC) :: GMCP0C
  DOUBLE PRECISION, DIMENSION(MXC,MXC) :: GMCP0D
end module comm_MCPOUT
!*MODULE MCPOUTG
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:37:39 UTC
!> @params GMCP1C
!> @params GMCP1D
module comm_MCPOUTG
  use mx_limits, only: MXC
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(MXC,3,MXC) :: GMCP1C
  DOUBLE PRECISION, DIMENSION(MXC,3,MXC) :: GMCP1D
end module comm_MCPOUTG
!*MODULE MCPWRKC
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:38:09 UTC
!> @params PSSC
!> @params PSPC
!> @params DSSC
!> @params DSPC
!> @params DSDC
!> @params FSSC
!> @params FSPC
!> @params FSDC
!> @params FSFC
!> @params GSSC
!> @params GSPC
!> @params GSDC
!> @params GSFC
!> @params GSGC
!> @params HSSC
!> @params HSPC
!> @params HSDC
!> @params HSFC
!> @params HSGC
!> @params HSHC
!> @params ISSC
!> @params ISPC
!> @params ISDC
!> @params ISFC
!> @params ISGC
!> @params ISHC
!> @params ISIC
module comm_MCPWRKC
  use mx_limits, only: MXT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(0:12,MXT,3,1) :: PSSC
  DOUBLE PRECISION, DIMENSION(0:11,MXT,3,3) :: PSPC
  DOUBLE PRECISION, DIMENSION(0:11,MXT,6,1) :: DSSC
  DOUBLE PRECISION, DIMENSION(0:10,MXT,6,3) :: DSPC
  DOUBLE PRECISION, DIMENSION(0: 9,MXT,6,6) :: DSDC
  DOUBLE PRECISION, DIMENSION(0:10,MXT,10,1) :: FSSC
  DOUBLE PRECISION, DIMENSION(0: 9,MXT,10,3) :: FSPC
  DOUBLE PRECISION, DIMENSION(0: 8,MXT,10,6) :: FSDC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,10,10) :: FSFC
  DOUBLE PRECISION, DIMENSION(0: 9,MXT,15,1) :: GSSC
  DOUBLE PRECISION, DIMENSION(0: 8,MXT,15,3) :: GSPC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,15,6) :: GSDC
  DOUBLE PRECISION, DIMENSION(0: 6,MXT,15,10) :: GSFC
  DOUBLE PRECISION, DIMENSION(0: 5,MXT,15,15) :: GSGC
  DOUBLE PRECISION, DIMENSION(0: 8,MXT,21,1) :: HSSC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,21,3) :: HSPC
  DOUBLE PRECISION, DIMENSION(0: 6,MXT,21,6) :: HSDC
  DOUBLE PRECISION, DIMENSION(0: 5,MXT,21,10) :: HSFC
  DOUBLE PRECISION, DIMENSION(0: 4,MXT,21,15) :: HSGC
  DOUBLE PRECISION, DIMENSION(0: 3,MXT,21,21) :: HSHC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,28,1) :: ISSC
  DOUBLE PRECISION, DIMENSION(0: 6,MXT,28,3) :: ISPC
  DOUBLE PRECISION, DIMENSION(0: 5,MXT,28,6) :: ISDC
  DOUBLE PRECISION, DIMENSION(0: 4,MXT,28,10) :: ISFC
  DOUBLE PRECISION, DIMENSION(0: 3,MXT,28,15) :: ISGC
  DOUBLE PRECISION, DIMENSION(0: 2,MXT,28,21) :: ISHC
  DOUBLE PRECISION, DIMENSION(0: 1,MXT,28,28) :: ISIC
end module comm_MCPWRKC
!*MODULE MCPWRKCG
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:38:40 UTC
!> @params SPSC
!> @params PPSC
!> @params PPPC
!> @params DPSC
!> @params DPPC
!> @params DPDC
!> @params FPSC
!> @params FPPC
!> @params FPDC
!> @params FPFC
!> @params GPSC
!> @params GPPC
!> @params GPDC
!> @params GPFC
!> @params GPGC
!> @params HPSC
!> @params HPPC
!> @params HPDC
!> @params HPFC
!> @params HPGC
!> @params HPHC
!> @params SSPC
!> @params PSDC
!> @params DSFC
!> @params FSGC
!> @params GSHC
!> @params HSIC
module comm_MCPWRKCG
  use mx_limits, only: MXT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(0:11,MXT,1,3,1) :: SPSC
  DOUBLE PRECISION, DIMENSION(0:10,MXT,3,3,1) :: PPSC
  DOUBLE PRECISION, DIMENSION(0: 9,MXT,3,3,3) :: PPPC
  DOUBLE PRECISION, DIMENSION(0: 9,MXT,6,3,1) :: DPSC
  DOUBLE PRECISION, DIMENSION(0: 8,MXT,6,3,3) :: DPPC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,6,3,6) :: DPDC
  DOUBLE PRECISION, DIMENSION(0: 8,MXT,10,3,1) :: FPSC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,10,3,3) :: FPPC
  DOUBLE PRECISION, DIMENSION(0: 6,MXT,10,3,6) :: FPDC
  DOUBLE PRECISION, DIMENSION(0: 5,MXT,10,3,10) :: FPFC
  DOUBLE PRECISION, DIMENSION(0: 7,MXT,15,3,1) :: GPSC
  DOUBLE PRECISION, DIMENSION(0: 6,MXT,15,3,3) :: GPPC
  DOUBLE PRECISION, DIMENSION(0: 5,MXT,15,3,6) :: GPDC
  DOUBLE PRECISION, DIMENSION(0: 4,MXT,15,3,10) :: GPFC
  DOUBLE PRECISION, DIMENSION(0: 3,MXT,15,3,15) :: GPGC
  DOUBLE PRECISION, DIMENSION(0: 6,MXT,21,3,1) :: HPSC
  DOUBLE PRECISION, DIMENSION(0: 5,MXT,21,3,3) :: HPPC
  DOUBLE PRECISION, DIMENSION(0: 4,MXT,21,3,6) :: HPDC
  DOUBLE PRECISION, DIMENSION(0: 3,MXT,21,3,10) :: HPFC
  DOUBLE PRECISION, DIMENSION(0: 2,MXT,21,3,15) :: HPGC
  DOUBLE PRECISION, DIMENSION(0: 1,MXT,21,3,21) :: HPHC
  DOUBLE PRECISION, DIMENSION(0:5,MXT,1,3) :: SSPC
  DOUBLE PRECISION, DIMENSION(0:5,MXT,3,6) :: PSDC
  DOUBLE PRECISION, DIMENSION(0:5,MXT,6,10) :: DSFC
  DOUBLE PRECISION, DIMENSION(0:5,MXT,10,15) :: FSGC
  DOUBLE PRECISION, DIMENSION(0:5,MXT,15,21) :: GSHC
  DOUBLE PRECISION, DIMENSION(0:5,MXT,21,28) :: HSIC
end module comm_MCPWRKCG
!*MODULE MCPWRKD
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:39:10 UTC
!> @params PSSD
!> @params PSPD
!> @params DSSD
!> @params DSPD
!> @params DSDD
!> @params FSSD
!> @params FSPD
!> @params FSDD
!> @params FSFD
!> @params GSSD
!> @params GSPD
!> @params GSDD
!> @params GSFD
!> @params GSGD
!> @params HSSD
!> @params HSPD
!> @params HSDD
!> @params HSFD
!> @params HSGD
!> @params HSHD
!> @params ISSD
!> @params ISPD
!> @params ISDD
!> @params ISFD
!> @params ISGD
!> @params ISHD
!> @params ISID
module comm_MCPWRKD
  use mx_limits, only: MXT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(MXT,3,1) :: PSSD
  DOUBLE PRECISION, DIMENSION(MXT,3,3) :: PSPD
  DOUBLE PRECISION, DIMENSION(MXT,6,1) :: DSSD
  DOUBLE PRECISION, DIMENSION(MXT,6,3) :: DSPD
  DOUBLE PRECISION, DIMENSION(MXT,6,6) :: DSDD
  DOUBLE PRECISION, DIMENSION(MXT,10,1) :: FSSD
  DOUBLE PRECISION, DIMENSION(MXT,10,3) :: FSPD
  DOUBLE PRECISION, DIMENSION(MXT,10,6) :: FSDD
  DOUBLE PRECISION, DIMENSION(MXT,10,10) :: FSFD
  DOUBLE PRECISION, DIMENSION(MXT,15,1) :: GSSD
  DOUBLE PRECISION, DIMENSION(MXT,15,3) :: GSPD
  DOUBLE PRECISION, DIMENSION(MXT,15,6) :: GSDD
  DOUBLE PRECISION, DIMENSION(MXT,15,10) :: GSFD
  DOUBLE PRECISION, DIMENSION(MXT,15,15) :: GSGD
  DOUBLE PRECISION, DIMENSION(MXT,21,1) :: HSSD
  DOUBLE PRECISION, DIMENSION(MXT,21,3) :: HSPD
  DOUBLE PRECISION, DIMENSION(MXT,21,6) :: HSDD
  DOUBLE PRECISION, DIMENSION(MXT,21,10) :: HSFD
  DOUBLE PRECISION, DIMENSION(MXT,21,15) :: HSGD
  DOUBLE PRECISION, DIMENSION(MXT,21,21) :: HSHD
  DOUBLE PRECISION, DIMENSION(MXT,28,1) :: ISSD
  DOUBLE PRECISION, DIMENSION(MXT,28,3) :: ISPD
  DOUBLE PRECISION, DIMENSION(MXT,28,6) :: ISDD
  DOUBLE PRECISION, DIMENSION(MXT,28,10) :: ISFD
  DOUBLE PRECISION, DIMENSION(MXT,28,15) :: ISGD
  DOUBLE PRECISION, DIMENSION(MXT,28,21) :: ISHD
  DOUBLE PRECISION, DIMENSION(MXT,28,28) :: ISID
end module comm_MCPWRKD
!*MODULE MCPWRKDG
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2020-07-26 10:39:40 UTC
!> @params SPSD
!> @params PPSD
!> @params PPPD
!> @params DPSD
!> @params DPPD
!> @params DPDD
!> @params FPSD
!> @params FPPD
!> @params FPDD
!> @params FPFD
!> @params GPSD
!> @params GPPD
!> @params GPDD
!> @params GPFD
!> @params GPGD
!> @params HPSD
!> @params HPPD
!> @params HPDD
!> @params HPFD
!> @params HPGD
!> @params HPHD
module comm_MCPWRKDG
  use mx_limits, only: MXT
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(MXT,1,3,1) :: SPSD
  DOUBLE PRECISION, DIMENSION(MXT,3,3,1) :: PPSD
  DOUBLE PRECISION, DIMENSION(MXT,3,3,3) :: PPPD
  DOUBLE PRECISION, DIMENSION(MXT,6,3,1) :: DPSD
  DOUBLE PRECISION, DIMENSION(MXT,6,3,3) :: DPPD
  DOUBLE PRECISION, DIMENSION(MXT,6,3,6) :: DPDD
  DOUBLE PRECISION, DIMENSION(MXT,10,3,1) :: FPSD
  DOUBLE PRECISION, DIMENSION(MXT,10,3,3) :: FPPD
  DOUBLE PRECISION, DIMENSION(MXT,10,3,6) :: FPDD
  DOUBLE PRECISION, DIMENSION(MXT,10,3,10) :: FPFD
  DOUBLE PRECISION, DIMENSION(MXT,15,3,1) :: GPSD
  DOUBLE PRECISION, DIMENSION(MXT,15,3,3) :: GPPD
  DOUBLE PRECISION, DIMENSION(MXT,15,3,6) :: GPDD
  DOUBLE PRECISION, DIMENSION(MXT,15,3,10) :: GPFD
  DOUBLE PRECISION, DIMENSION(MXT,15,3,15) :: GPGD
  DOUBLE PRECISION, DIMENSION(MXT,21,3,1) :: HPSD
  DOUBLE PRECISION, DIMENSION(MXT,21,3,3) :: HPPD
  DOUBLE PRECISION, DIMENSION(MXT,21,3,6) :: HPDD
  DOUBLE PRECISION, DIMENSION(MXT,21,3,10) :: HPFD
  DOUBLE PRECISION, DIMENSION(MXT,21,3,15) :: HPGD
  DOUBLE PRECISION, DIMENSION(MXT,21,3,21) :: HPHD
end module comm_MCPWRKDG
!*MODULE MRCOM
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2022-04-04 09:59:46 UTC
!> @params ZVTOL
module comm_MRCOM
    IMPLICIT NONE
    DOUBLE PRECISION :: ZVTOL
end module comm_MRCOM
!*MODULE MRSYM
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2022-04-04 10:00:03 UTC
!> @params STSYM
module comm_MRSYM
    IMPLICIT NONE
    DOUBLE PRECISION :: STSYM
end module comm_MRSYM
!*MODULE MRREW
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Konstantin Komarov
!> @date   Generated at 2022-05-26
!> @params IXCORE
module comm_INFOXR
    use mx_limits, only: mxao
    IMPLICIT NONE
    INTEGER, DIMENSION(mxao) :: IXCORE
end module comm_INFOXR
!*MODULE NONAD
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2022-04-04 10:00:20 UTC
!> @params NAMD
!> @params NDSWCH
!> @params NDRST
!> @params NDTLF
!> @params COLD
!> @params THRSHE
module comm_NONAD
    USE MX_LIMITS, ONLY: MXATM
    IMPLICIT NONE
    LOGICAL :: NAMD
    LOGICAL :: NDSWCH
    LOGICAL :: NDRST
    INTEGER :: NDTLF
    DOUBLE PRECISION :: COLD(3,MXATM)
    DOUBLE PRECISION :: THRSHE
end module comm_NONAD
!*MODULE REXOPT
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2022-04-04 10:00:37 UTC
!> @params REXTYPE
!> @params REXTARGET
!> @params REXSHIFT
!> @params REXDIIS
!> @params REXLDL
!> @params RLXDEN
!> @params REXEKT
!> @params EKTEA
!> @params REXCG
!> @params RXCGIT
!> @params RXCGTH
module comm_REXOPT
    IMPLICIT NONE
    INTEGER :: REXTYPE
    INTEGER :: REXTARGET
    DOUBLE PRECISION :: REXSHIFT
    LOGICAL :: REXDIIS
    INTEGER :: REXLDL
    LOGICAL :: RLXDEN
    LOGICAL :: REXEKT
    LOGICAL :: EKTEA
    LOGICAL :: REXCG
    INTEGER :: RXCGIT
    DOUBLE PRECISION :: RXCGTH
end module comm_REXOPT
!*MODULE REKSCM
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2022-04-04 10:00:55 UTC
!> @params NMICRO
!> @params MTTYP
!> @params WPPS
!> @params WOSS
!> @params G1
!> @params DNR
!> @params DNS
!> @params DELTA
!> @params FR
!> @params FS
module comm_REKSCM
    IMPLICIT NONE
    INTEGER :: NMICRO(2,2,4)
    INTEGER :: MTTYP
    DOUBLE PRECISION :: WPPS
    DOUBLE PRECISION :: WOSS
    DOUBLE PRECISION :: G1
    DOUBLE PRECISION :: DNR
    DOUBLE PRECISION :: DNS
    DOUBLE PRECISION :: DELTA
    DOUBLE PRECISION :: FR
    DOUBLE PRECISION :: FS
end module comm_REKSCM
!*MODULE SSR
!> @brief  NA
!> @detail this module is analog of old common block
!> @author Igor S. Gerasimov
!> @date   Generated at 2022-04-04 10:01:12 UTC
!> @params S2SARE
!> @params WS2SA
!> @params S3SARE
!> @params WS3SA
!> @params CLXGRD
module comm_SSR
    IMPLICIT NONE
    DOUBLE PRECISION :: S2SARE(2,2)
    DOUBLE PRECISION :: WS2SA(2)
    DOUBLE PRECISION :: S3SARE(3,3)
    DOUBLE PRECISION :: WS3SA(3)
    DOUBLE PRECISION :: CLXGRD(4)
end module comm_SSR
