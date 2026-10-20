C*MODULE MRSF   *DECK MRSFDFTCALC
C>
C>    @brief   Mixed-Reference Spin-Flip (MRSF)-TDDFT calculation
C>
C>    @details Driver for MRSF-TDDFT based on
C>             JCP 149 104101 (2018), JCP 150 184111 (2019)
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE MRSFDFTCALC
      USE comm_MRSYM, ONLY: STSYM
C
      USE lrcdft, ONLY: LCFLAG
      USE MX_LIMITS, ONLY: mxrt, mxgrid, mxatm
      USE constants, only: one
      USE comm_EFPTDG
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      LOGICAL :: ALPHKWD, BETAKWD, DSKWRK, GOPARR, MASWRK, MRDEA, MREKT,&
     &           SG1, SG1T, TAMMD, TPA, TRIPLET
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp) :: CNVTOL, DFTGTHR, DFTTHR, EXENA, EXENB, EXENC,     &
     &                EXETYP, RUNTYP, SW0, SWOFF
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      INTEGER :: IBTYP, ICUT, IDAF, IDFT34, IFMODIM, IP, IPK, IPTIM, IR,&
     &           IRECTD, IS, ITDFG, ITDPRP, ITOL, IW, JANS, JANST,      &
     &           MASTER, MAXVEC, ME, MODTD, MONOC1, MONOC2, MONVR1,     &
     &           MONVR2, MTHST, MULTD, NAUXFUN, NAUXSHL, NAV, NDFTFG,   &
     &           NEVALS, NGLEVL, NHLEVL, NLEBT, NONEQR, NOPK, NORMF,    &
     &           NORMP, NPHI, NPHI0, NPHIT, NPRINT, NPROC, NRAD, NRAD0, &
     &           NRADT, NSTAT, NTHE, NTHE0, NTHET, NTHST, NTRIAL
      INTEGER, DIMENSION(4) :: IFEDAT
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(MXRT) :: MOCC, MVIR
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0, NLEB, NLEB0
      REAL(KIND=dp), DIMENSION(2) :: PFREQ
      REAL(KIND=dp), DIMENSION(3) :: SPCP
      REAL(KIND=dp), DIMENSION(3,MXRT) :: TDM
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /DFLEB / NLEB, NLEB0
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /INFOEX/ TDM, MOCC, MVIR, MONOC1, MONVR1, IFMODIM, MONOC2, &
     &                MONVR2
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      COMMON /RUNOPT/ RUNTYP, EXETYP, NEVALS, NGLEVL, NHLEVL
C
      LOGICAL :: GOPSAV, LMRSF, MRSFQ, MRSFS, MRSFT, SOME
      REAL(KIND=dp) :: HFSCAL, SPCP1, SPCP2, SPCP3
      INTEGER :: NCONV, NDER, NLEBS, NPHIS, NPROCS, NRADS, NTHES
      INTEGER, SAVE :: NTHST_SAVE
C
      LOGICAL, SAVE :: FIRST = .TRUE.
      INTEGER, SAVE :: MXITER = 200
C
      REAL(KIND=dp), PARAMETER :: BLANK = transfer('        ',1.0d0)
      REAL(KIND=dp), PARAMETER :: CHECK = transfer('CHECK   ',1.0d0)
      REAL(KIND=dp), PARAMETER :: ROHF = transfer('ROHF    ',1.0d0)
      REAL(KIND=dp), PARAMETER :: SFDFT = transfer('SFDFT   ',1.0d0)
C
      SOME = MASWRK .AND. NPRINT.NE.-5
C
      LMRSF = .TRUE.
      MRSFS = .FALSE.
      MRSFT = .FALSE.
      MRSFQ = .FALSE.
      IF(MULTD.EQ.1) MRSFS = .TRUE.
      IF(MULTD.EQ.3) MRSFT = .TRUE.
      IF(MULTD.EQ.5) MRSFQ = .TRUE.
C
C     ----- EXTRACT SOME INPUT PARAMETERS FROM $TDDFT -----
C
      IF(SOME) THEN
        IF(MRSFQ)
     *     WRITE(IW,9000) NSTAT,NTHST,MAXVEC,NTRIAL,CNVTOL
        IF(MRSFS .OR. MRSFT) THEN
           HFSCAL = DFTTYP(3)
           SPCP1 = SPCP(1)
           SPCP2 = SPCP(2)
           SPCP3 = SPCP(3)
           IF(SPCP1.EQ.-ONE) SPCP1 = HFSCAL
           IF(SPCP2.EQ.-ONE) SPCP2 = HFSCAL
           IF(SPCP3.EQ.-ONE) SPCP3 = HFSCAL
C
           WRITE(IW,9001) NSTAT,NTHST,MAXVEC,NTRIAL,CNVTol,
     *                    HFSCAL,SPCP1,SPCP2,SPCP3
        ENDIF
        WRITE(IW,9011)
      END IF
C
C****                             ****************************************
C**** MRSF-DFT ENERGY CALCULATION ****************************************
C****                             ****************************************
      CALL VCLR(TDM,1,3*MXRT)

      IF(FIRST) THEN
         FIRST = .FALSE.
      ELSE IF(.NOT.FIRST .and. STSYM .ne. blank) THEN
         NTHST = NTHST_SAVE
      ENDIF
      CALL SPNFLPCL(SFDFT,CNVTOL,MOCC,MVIR,NSTAT,
     *              NTRIAL,MAXVEC,NTHST,MXITER,NCONV,
     *              LMRSF,MRSFS,MRSFT,MRSFQ,SPCP,
     *              NTHST_SAVE)
C
C        IF NO EXCITED STATE DENSITY WILL BE GENERATED, JUST QUIT.
C
      CALL DERCHK(NDER)
c     IF(NDER.EQ.0  .AND.  ITDPRP.EQ.0) RETURN
      if((nder.gt.0 .and. nglevl.eq.0) .or. itdprp.ne.0) then
         continue
      else
         return
      end if
C
C****                             **************************************
C**** SF-DFT GRADIENT CALCULATION **************************************
C****                             **************************************
      IF(EXETYP.EQ.CHECK) NCONV = 0
      IF(NCONV.NE.0) THEN
         IF(MASWRK) WRITE(IW,*)
     *     ' NO GRADIENT. MRSF-DFT ENERGY CALCULATION IS NOT CONVERGED.'
         CALL ABRT
      END IF
C
C     INITIALIZE SF-DFT/EFP1 GRADIENT
C
      IF(ITDEFG.NE.0) CALL TDEFPINIT
C
C     PARALLEL VERSION OF THE GRADIENT CALCULATION ONLY WORKS FOR
C     DYANAMIC LOAD BALANCING
C
      GOPSAV = GOPARR
      NPROCS = NPROC
C
      IF(SOME) WRITE(IW,9101)
C
      IF(NTHST.GT.MXRT) THEN
         IF(MASWRK .AND. .not.LMRSF)
     *      WRITE(IW,*) 'SFDFT-GRADIENT CALCULATION IS ALLOWED',
     *                  'ONLY FOR IROOT.LE.',MXRT
         IF(MASWRK .AND. LMRSF)
     *      WRITE(IW,*) 'MRSFDFT-GRADIENT CALCULATION IS ALLOWED',
     *                  'ONLY FOR IROOT.LE.',MXRT
         CALL ABRT
      END IF
C
C     -- SWITCH TO THE SMALL GRID SIZE FOR TDDFT
C
      NRADS   = NRAD
      NTHES   = NTHE
      NPHIS   = NPHI
      NLEBS   = NLEB(1)
      NRAD    = NRADT
      NTHE    = NTHET
      NPHI    = NPHIT
      NLEB(1) = NLEBT
C
C     POINT GROUP SYMMETRY CANNOT BE USED DURING THE RESPONSES.
C     THE SELECTED STATE MAY NOT HAVE A SYMMETRIC DENSITY, SO
C     WE ALSO HAVE TO DO THE GRADIENT INTEGRAL TERMS (LATER)
C     W/O SYMMETRY.
C
      CALL SYMOFF
      CALL SFROGRAD(SFDFT,NTHST,LMRSF,MRSFS,MRSFT,MRSFQ,SPCP,
     *              MREKT,MRDEA)
      CALL SYMON
C
      IF(MASWRK)      WRITE(IW,9501)
      CALL TIMIT(1)
C
C     ----- CALCULATE PROPERTY -----
C
      IF(NPRINT.NE.-5) CALL PROPTY('TDDF')
C     NATURAL ORBITAL
      IF (ITDPRP.NE.0) THEN
         CALL SYMOFF
         CALL TDDNOS
         CALL SYMON
      ENDIF
C
C     -- SWITCH BACK TO THE REGULAR GRID SIZE, E.G. FOR GRADIENT STEP
C
      NRAD = NRADS
      NTHE = NTHES
      NPHI = NPHIS
      NLEB(1) = NLEBS
C
      NPROC  = NPROCS
      GOPARR = GOPSAV
C
      RETURN
 9000 FORMAT(/10X,47(1H-)/10X,'SF-DFT INPUT PARAMETERS ',
     *            '(EXTRACTED FROM $TDDFT)'/10X,47(1H-)//
     *        5X,'NSTATE=',I8,'  IROOT=',I8,' MAXVEC=',I8/
     *        5X,'NTRIAL=',I8,' CNVTOL=',1P,E10.2)
 9001 FORMAT(/8X,50(1H-)/8X,'MRSF-DFT INPUT PARAMETERS ',
     *            '(EXTRACTED FROM $TDDFT)'/8X,50(1H-)//
     *        5X,'NSTATE=',I8,'  IROOT=',I8,' MAXVEC=',I8/
     *        5X,'NTRIAL=',I8,' CNVTOL=',1P,E10.2//
     *        5X,'FITTING PARAMETERS OF SPIN-PAIRING COUPLINGS'/
     *        7X,'   cHF    ','  CO-CO   ','  OV-OV   ','  CO-OV   '/
     *        5X,0P,F10.5,F10.5,F10.5,F10.5)
 9011 FORMAT(/5X,'********* NOTE: THE CURRENT STATUS OF MRSFDFT',
     *           ' *********'//
     *        5X,'(1) ROHF REFERENCE IS AVAILABLE.'/
     *        5X,'(2) COLLINEAR APPROXIMATION FOR XC FUNCTIONALS.'/
     *       /5X,55(1H*)/)
 9101 FORMAT(/1X,71(1H-)/
     *  17X,'MRSF-TDDFT ENERGY GRADIENT CALCULATION'/
     *       1X,71(1H-)/)
 9501 FORMAT(1X,'..... DONE WITH MRSF-DFT EXCITED STATE RESPONSE',
     *          ' AND DENSITY MATRIX .....')
      END SUBROUTINE MRSFDFTCALC
c*MODULE MRSF   *deck mrsf2e
C>
C>    @brief   clone of SF2E
C>             for singlet and triplet MRSF
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Terms co2v, co1v, cco1 and cco2 was replaced
C>             by o21v, co12. mrsfcbc is initial subroutine of this.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      subroutine mrsf2e(bo2v,bo1v,bco1,bco2,ball,o21v,co12,
     &                  agdlr,ao21v,aco12,
     &                  ado2v,ado1v,adco1,
     &                  adco2,buf,ibuf,nbf,nv,nmax)
C
      USE mx_limits, only: mxsh, mxgtot, mxao
      USE constants, only: one
      USE lrcdft, ONLY: LCFLAG
      USE camdft, ONLY: ALPHAC => cam_alpha, BETAC => cam_beta, CAMFLAG
      USE prec, ONLY: dp
      implicit none
c
      real(kind=dp), dimension(mxgtot) :: cd, cf, cg, ch, ci, cp, cs, ex
      logical :: dirscf, dskwrk, fdiff, goparr, lrint, maswrk, pack2e
      integer, dimension(mxao) :: ia
      integer :: ibtyp, icut, idaf, iecp, iefld, igrdtyp, ijk, inttyp,  &
     &           ip, ipk, iptim, ir, ischwz, itol, iw, master, me, nav, &
     &           necp, nhex, nintmx, nopk, normf, normp, nprint, nproc, &
     &           nshell, ntupl
      integer, dimension(950) :: ioda
      integer, dimension(mxsh) :: katom, kloc, kmax, kmin, kng, kstart, &
     &                            ktype
      real(kind=dp), dimension(1) :: xx
      common /fmcom / xx
      common /ijpair/ ia
      common /intfil/ nintmx, nhex, ntupl, pack2e, inttyp, igrdtyp
      common /intopt/ ischwz, iecp, necp, iefld
      common /iofile/ ir, iw, ip, ijk, ipk, idaf, nav, ioda
      common /nlrcf / lrint
      common /nshel / ex, cs, cp, cd, cf, cg, ch, ci, kstart, katom,    &
     &                ktype, kng, kloc, kmin, kmax, nshell
      common /optscf/ dirscf, fdiff
      common /output/ nprint, itol, icut, normf, normp, nopk
      common /par   / me, master, nproc, ibtyp, iptim, goparr, dskwrk,  &
     &                maswrk
c
      integer :: nbf, nmax, nv
      real(kind=dp), dimension(*) :: agdlr, adco1, adco2, ado1v, ado2v, &
     &                               aco12, ao21v
      real(kind=dp), dimension(*) :: ball, bco1, bco2, bo1v, bo2v,      &
     &                               o21v, co12
      real(kind=dp), dimension(nintmx) :: buf
      integer, dimension(nintmx) :: ibuf
c
      real(kind=dp) :: dummy, scalfac
      integer :: last, ldsh, lghond, lmax, loadfm, lxints, maxg, nangm, &
     &           nbf3, need, nint, nschwz, nsh2
      logical :: schwrz, tdskwrk
C
      real(kind=dp), parameter :: rhf = transfer('RHF     ',1.0d0)
c
c     --- form square non-symmetric fock-like matrix ---
c         direct method = recompute 2e- ao integrals
c         standard method = process integrals from disk
c
      nint = 0
      nschwz = 0
      nbf3 = nbf*nbf

      if(dirscf) then

         schwrz=ischwz.gt.0
         nsh2 = (nshell*nshell+nshell)/2
c
         call baschk(lmax)
                       nangm =  4
         if(lmax.eq.2) nangm =  6
         if(lmax.eq.3) nangm = 10
         if(lmax.eq.4) nangm = 15
         if(lmax.eq.5) nangm = 21
         if(lmax.eq.6) nangm = 28
         maxg = nangm**4
c
         call valfm(loadfm)
         lghond = loadfm + 1
         lxints = lghond + maxg
         ldsh   = lxints + nsh2
         last   = ldsh     + nsh2
         need   = last - loadfm - 1
         call getfm(need)
c
         if(schwrz) then
            dummy = 0.0d+00
            call shltd(rhf,ball,dummy,xx(ldsh),nbf,nsh2,nv)
            call daread(idaf,ioda,xx(lxints),nsh2,54,0)
C
         end if
         call mrsftwoei(schwrz,nint,nschwz,nbf,xx(lxints),nsh2,
     *                  xx(lghond),maxg,ia,ball,agdlr,xx(ldsh),
     *                  nv,.true.,nmax,bo2v,bo1v,bco1,bco2,ball,
     *                  o21v,co12,
     *                  agdlr,ao21v,aco12,
     *                  ado2v,ado1v,adco1,adco2)
         call retfm(need)
      else
         tdskwrk = dskwrk
         dskwrk = .true.
         call seqrew(ijk)
         if(lcflag) then
           lrint = .true.
           call mrsfadisk(bo2v,bo1v,bco1,bco2,ball,
     *                    o21v,co12,
     *                    nbf,agdlr,ao21v,aco12,
     *                    ado2v,ado1v,adco1,adco2,
     *                    buf,ibuf,nintmx,nopk,nv)
           lrint = .false.
         endif
         if(CAMFLAG) then
          LRINT = .TRUE.
          SCALFAC = BETAC/ALPHAC
          call mrsfadisk(bo2v,bo1v,bco1,bco2,ball,
     *                   o21v,co12,
     *                   nbf,agdlr,ao21v,aco12,
     *                   ado2v,ado1v,adco1,adco2,
     *                   buf,ibuf,nintmx,nopk,nv)
          if(scalfac.ne.one) call dscal(nbf3*nv,scalfac,agdlr,1)
          lrint = .false.
         endif
         call mrsfadisk(bo2v,bo1v,bco1,bco2,ball,
     *                  o21v,co12,
     *                  nbf,agdlr,ao21v,aco12,
     *                  ado2v,ado1v,adco1,adco2,
     *                  buf,ibuf,nintmx,nopk,nv)
         dskwrk = tdskwrk
      end if
c
c     --- sum up partial fock-like matrices ---
c
      if(goparr) then
         call ddi_gsumi(2311,nint,1)
         call ddi_gsumi(2312,nschwz,1)
         call ddi_gsumf(2313,agdlr,7*nmax*nbf3)
      end if
c
      return
      end subroutine mrsf2e
c*MODULE MRSF   *deck mrsfadisk
C>
C>    @brief   clone of SFADISK
C>             for singlet and triplet MRSF
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Terms co2v, co1v, cco1 and cco2 was replaced
C>             by o21v, co12. mrsfcbc is initial subroutine of this.
C>             The integral cycle was reduced by 16 lines.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      subroutine mrsfadisk(bo2v,bo1v,bco1,bco2,ball,
     &                     o21v,co12,
     &                     l1,agdlr,ao21v,aco12,
     &                     ado2v,ado1v,adco1,adco2,
     &                     xx,ix,nintmx,nopk,nv)
C
      USE mx_limits, ONLY: mxao
      USE lrcdft, ONLY: LCFLAG, LRFILE
      USE prec, ONLY: dp
      USE constants, ONLY: half
      IMPLICIT NONE
C
      INTEGER, DIMENSION(mxao) :: IA
      INTEGER :: IDAF, IP, IPK, IR, IS, IW, LABSIZ, NAV
      INTEGER, DIMENSION(950) :: IODA
      LOGICAL :: LRINT
      INTEGER, DIMENSION(7) :: NORDER
      COMMON /IJPAIR/ IA
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      COMMON /NLRCF / LRINT
      COMMON /ORDOPT/ NORDER
      COMMON /PCKLAB/ LABSIZ
C
      INTEGER :: L1, NINTMX, NOPK, NV
      REAL(KIND=dp), DIMENSION(L1,L1,*) :: ADCO1, ADCO2, ADO1V, ADO2V,
     &                                     AGDLR, AO21V, ACO12
      REAL(KIND=dp), DIMENSION(L1,L1,*) :: BCO1, BCO2, BO1V, BO2V,
     &                                     BALL, O21V, CO12
      INTEGER, DIMENSION(*) :: IX
      REAL(KIND=dp), DIMENSION(nintmx) :: XX
C
      INTEGER :: I, IMO, IPACK, J, JPACK, K, KPACK, L, LABEL, LPACK,    &
     &           M, NIJ, NINT, NKL, NPACK, NXX
      REAL(KIND=dp) :: VAL
c
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
      i = 0
      j = 0
      k = 0
      l = 0
      nxx = 1 ! for entering to the cycle
      if(nopk.ne.1) then
         write(iw,*) 'nopk.ne.1'
         write(iw,*) 'sfdft calculation does not support ',
     *               'this integral type'
         call abrt
      end if
c
c     ----- integrals are not in supermatrix form (nopk=.true.) -----
c
      do while(nxx .gt. 0)
      if(.not.LRINT) call pread(is,xx,ix,nxx,nintmx)
      if(     LRINT) call pread(lrfile,xx,ix,nxx,nintmx)
      if(nxx .ne. 0) then
         nint = iabs(nxx)
         if(nint .gt. nintmx) call abrt
         do imo = 1,nv
            do m = 1,nint
               npack = m
               if (labsiz .eq. 2) then
                  label = ix(npack)
                  ipack = ishft( label, -48 )
                  jpack = iand( ishft( label, -32 ), 65535 )
                  kpack = iand( ishft( label, -16 ), 65535 )
                  lpack = iand( label, 65535 )
               else if (labsiz .eq. 1) then
                  if ( mod(npack,2) .eq. 0 ) then
                     label = ix( npack/2 )
                     ipack = iand( ishft( label, -24 ), 255 )
                     jpack = iand( ishft( label, -16 ), 255 )
                     kpack = iand( ishft( label,  -8 ), 255 )
                     lpack = iand( label, 255 )
                  else
                     label = ix( (npack/2)+1 )
                     ipack = ishft( label, -56 )
                     jpack = iand( ishft( label, -48 ), 255 )
                     kpack = iand( ishft( label, -40 ), 255 )
                     lpack = iand( ishft( label, -32 ), 255 )
                  end if
               end if
               i = ipack
               j = jpack
               k = kpack
               l = lpack

               val = xx(m)
               nij = i*(i-1)+j
               nkl = k*(k-1)+l
c
c              using square canonical integral file
c
               if(norder(7) .eq. 1) then
                  if(nkl .gt. nij) CYCLE
                  if(i .eq. j) val = val*half
                  if(k .eq. l) val = val*half
                  if(nij .eq. nkl) val = val*half
               end if

               if(.not.LCFLAG .or. LRINT) then
                  agdlr(i,k,imo) = agdlr(i,k,imo) - val*ball(j,l,imo)
                  agdlr(k,i,imo) = agdlr(k,i,imo) - val*ball(l,j,imo)
                  agdlr(i,l,imo) = agdlr(i,l,imo) - val*ball(j,k,imo)
                  agdlr(l,i,imo) = agdlr(l,i,imo) - val*ball(k,j,imo)
                  agdlr(j,k,imo) = agdlr(j,k,imo) - val*ball(i,l,imo)
                  agdlr(k,j,imo) = agdlr(k,j,imo) - val*ball(l,i,imo)
                  agdlr(j,l,imo) = agdlr(j,l,imo) - val*ball(i,k,imo)
                  agdlr(l,j,imo) = agdlr(l,j,imo) - val*ball(k,i,imo)
               endif

               if(.not.LRINT) then
                  ao21v(i,l,imo) = ao21v(i,l,imo) - val*o21v(k,j,imo)! (ij|lk)
                  ao21v(k,j,imo) = ao21v(k,j,imo) - val*o21v(i,l,imo)! (kl|ji)
                  ao21v(i,k,imo) = ao21v(i,k,imo) - val*o21v(l,j,imo)! (ij|kl)
                  ao21v(l,j,imo) = ao21v(l,j,imo) - val*o21v(i,k,imo)! (lk|ji)
                  ao21v(j,l,imo) = ao21v(j,l,imo) - val*o21v(k,i,imo)! (ji|lk)
                  ao21v(k,i,imo) = ao21v(k,i,imo) - val*o21v(j,l,imo)! (kl|ij)
                  ao21v(j,k,imo) = ao21v(j,k,imo) - val*o21v(l,i,imo)! (ji|kl)
                  ao21v(l,i,imo) = ao21v(l,i,imo) - val*o21v(j,k,imo)! (lk|ij)

                  aco12(i,l,imo) = aco12(i,l,imo) - val*co12(k,j,imo)! (ij|lk)
                  aco12(k,j,imo) = aco12(k,j,imo) - val*co12(i,l,imo)! (kl|ji)
                  aco12(i,k,imo) = aco12(i,k,imo) - val*co12(l,j,imo)! (ij|kl)
                  aco12(l,j,imo) = aco12(l,j,imo) - val*co12(i,k,imo)! (lk|ji)
                  aco12(j,l,imo) = aco12(j,l,imo) - val*co12(k,i,imo)! (ji|lk)
                  aco12(k,i,imo) = aco12(k,i,imo) - val*co12(j,l,imo)! (kl|ij)
                  aco12(j,k,imo) = aco12(j,k,imo) - val*co12(l,i,imo)! (ji|kl)
                  aco12(l,i,imo) = aco12(l,i,imo) - val*co12(j,k,imo)! (lk|ij)

                  ado2v(i,k,imo) = ado2v(i,k,imo) - val*bo2v(j,l,imo)! (ij|lk)
                  ado2v(k,i,imo) = ado2v(k,i,imo) - val*bo2v(l,j,imo)! (kl|ji)
                  ado2v(i,l,imo) = ado2v(i,l,imo) - val*bo2v(j,k,imo)! (ij|kl)
                  ado2v(l,i,imo) = ado2v(l,i,imo) - val*bo2v(k,j,imo)! (lk|ji)
                  ado2v(j,k,imo) = ado2v(j,k,imo) - val*bo2v(i,l,imo)! (ji|lk)
                  ado2v(k,j,imo) = ado2v(k,j,imo) - val*bo2v(l,i,imo)! (kl|ij)
                  ado2v(j,l,imo) = ado2v(j,l,imo) - val*bo2v(i,k,imo)! (ji|kl)
                  ado2v(l,j,imo) = ado2v(l,j,imo) - val*bo2v(k,i,imo)! (lk|ij)
                  ado2v(i,j,imo) = ado2v(i,j,imo) + val*bo2v(k,l,imo)! (ij|lk)
                  ado2v(k,l,imo) = ado2v(k,l,imo) + val*bo2v(i,j,imo)! (kl|ji)
                  ado2v(i,j,imo) = ado2v(i,j,imo) + val*bo2v(l,k,imo)! (ij|kl)
                  ado2v(l,k,imo) = ado2v(l,k,imo) + val*bo2v(i,j,imo)! (lk|ji)
                  ado2v(j,i,imo) = ado2v(j,i,imo) + val*bo2v(k,l,imo)! (ji|lk)
                  ado2v(k,l,imo) = ado2v(k,l,imo) + val*bo2v(j,i,imo)! (kl|ij)
                  ado2v(j,i,imo) = ado2v(j,i,imo) + val*bo2v(l,k,imo)! (ji|kl)
                  ado2v(l,k,imo) = ado2v(l,k,imo) + val*bo2v(j,i,imo)! (lk|ij)

                  ado1v(i,k,imo) = ado1v(i,k,imo) - val*bo1v(j,l,imo)
                  ado1v(k,i,imo) = ado1v(k,i,imo) - val*bo1v(l,j,imo)
                  ado1v(i,l,imo) = ado1v(i,l,imo) - val*bo1v(j,k,imo)
                  ado1v(l,i,imo) = ado1v(l,i,imo) - val*bo1v(k,j,imo)
                  ado1v(j,k,imo) = ado1v(j,k,imo) - val*bo1v(i,l,imo)
                  ado1v(k,j,imo) = ado1v(k,j,imo) - val*bo1v(l,i,imo)
                  ado1v(j,l,imo) = ado1v(j,l,imo) - val*bo1v(i,k,imo)
                  ado1v(l,j,imo) = ado1v(l,j,imo) - val*bo1v(k,i,imo)
                  ado1v(i,j,imo) = ado1v(i,j,imo) + val*bo1v(k,l,imo)! (ij|lk)
                  ado1v(k,l,imo) = ado1v(k,l,imo) + val*bo1v(i,j,imo)! (kl|ji)
                  ado1v(i,j,imo) = ado1v(i,j,imo) + val*bo1v(l,k,imo)! (ij|kl)
                  ado1v(l,k,imo) = ado1v(l,k,imo) + val*bo1v(i,j,imo)! (lk|ji)
                  ado1v(j,i,imo) = ado1v(j,i,imo) + val*bo1v(k,l,imo)! (ji|lk)
                  ado1v(k,l,imo) = ado1v(k,l,imo) + val*bo1v(j,i,imo)! (kl|ij)
                  ado1v(j,i,imo) = ado1v(j,i,imo) + val*bo1v(l,k,imo)! (ji|kl)
                  ado1v(l,k,imo) = ado1v(l,k,imo) + val*bo1v(j,i,imo)! (lk|ij)

                  adco1(i,k,imo) = adco1(i,k,imo) - val*bco1(j,l,imo)
                  adco1(k,i,imo) = adco1(k,i,imo) - val*bco1(l,j,imo)
                  adco1(i,l,imo) = adco1(i,l,imo) - val*bco1(j,k,imo)
                  adco1(l,i,imo) = adco1(l,i,imo) - val*bco1(k,j,imo)
                  adco1(j,k,imo) = adco1(j,k,imo) - val*bco1(i,l,imo)
                  adco1(k,j,imo) = adco1(k,j,imo) - val*bco1(l,i,imo)
                  adco1(j,l,imo) = adco1(j,l,imo) - val*bco1(i,k,imo)
                  adco1(l,j,imo) = adco1(l,j,imo) - val*bco1(k,i,imo)
                  adco1(i,j,imo) = adco1(i,j,imo) + val*bco1(k,l,imo)! (ij|lk)
                  adco1(k,l,imo) = adco1(k,l,imo) + val*bco1(i,j,imo)! (kl|ji)
                  adco1(i,j,imo) = adco1(i,j,imo) + val*bco1(l,k,imo)! (ij|kl)
                  adco1(l,k,imo) = adco1(l,k,imo) + val*bco1(i,j,imo)! (lk|ji)
                  adco1(j,i,imo) = adco1(j,i,imo) + val*bco1(k,l,imo)! (ji|lk)
                  adco1(k,l,imo) = adco1(k,l,imo) + val*bco1(j,i,imo)! (kl|ij)
                  adco1(j,i,imo) = adco1(j,i,imo) + val*bco1(l,k,imo)! (ji|kl)
                  adco1(l,k,imo) = adco1(l,k,imo) + val*bco1(j,i,imo)! (lk|ij)

                  adco2(i,k,imo) = adco2(i,k,imo) - val*bco2(j,l,imo)
                  adco2(k,i,imo) = adco2(k,i,imo) - val*bco2(l,j,imo)
                  adco2(i,l,imo) = adco2(i,l,imo) - val*bco2(j,k,imo)
                  adco2(l,i,imo) = adco2(l,i,imo) - val*bco2(k,j,imo)
                  adco2(j,k,imo) = adco2(j,k,imo) - val*bco2(i,l,imo)
                  adco2(k,j,imo) = adco2(k,j,imo) - val*bco2(l,i,imo)
                  adco2(j,l,imo) = adco2(j,l,imo) - val*bco2(i,k,imo)
                  adco2(l,j,imo) = adco2(l,j,imo) - val*bco2(k,i,imo)
                  adco2(i,j,imo) = adco2(i,j,imo) + val*bco2(k,l,imo)! (ij|lk)
                  adco2(k,l,imo) = adco2(k,l,imo) + val*bco2(i,j,imo)! (kl|ji)
                  adco2(i,j,imo) = adco2(i,j,imo) + val*bco2(l,k,imo)! (ij|kl)
                  adco2(l,k,imo) = adco2(l,k,imo) + val*bco2(i,j,imo)! (lk|ji)
                  adco2(j,i,imo) = adco2(j,i,imo) + val*bco2(k,l,imo)! (ji|lk)
                  adco2(k,l,imo) = adco2(k,l,imo) + val*bco2(j,i,imo)! (kl|ij)
                  adco2(j,i,imo) = adco2(j,i,imo) + val*bco2(l,k,imo)! (ji|kl)
                  adco2(l,k,imo) = adco2(l,k,imo) + val*bco2(j,i,imo)! (lk|ij)
               endif
            end do
         end do
c
      END IF
      end do
c
      if(.not.LRINT) call seqrew(is)
      if(     LRINT) call seqrew(LRFILE)
c
      return
      end subroutine mrsfadisk
C
C*MODULE MRSF   *DECK MRSFTWOEI
C>
C>    @brief  clone of SFTWOEI
C>            for singlet and triplet MRSF
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Terms co2v, co1v, cco1 and cco2 was replaced
C>             by o21v, co12. mrsfcbc is initial subroutine of this.
C>             Updated dirfck procedure for mrsf. ii,jj,kk,ll and val are written
C>             to the buffer (call mrsfaddtobuf) until it is full
C>             (MXBUF = 1E+6 ~ 8Mb of cash memory) then the flush buffer
C>             (mrsfflushbuf) is called to calculate the integrals and update
C>             the fock matrix. In case the buffer is not filled after the loop
C>             to the end, mrsfflushbuf is called once again to finish.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE MRSFTWOEI(SCHWRZ,NINT,NSCHWZ,L1,XINTS,NSH2,GHONDO,MAXG,
     &                     IA,DA,F,DSH,NV,lmrsf,nmax,
     &                     bo2v,bo1v,bco1,bco2,ball,
     &                     o21v,co12,
     &                     agdlr,ao21v,aco12,
     &                     ado2v,ado1v,adco1,adco2)
C
C     THIS IS A DUMMY CODE OF 'SUBROUTINE TWOEI' @ INT2A.SRC
C
      USE MX_LIMITS, ONLY: MXSH, MXGTOT, MXATM
      USE prec, ONLY: dp
      USE constants, ONLY: zero
      IMPLICIT NONE
C
      LOGICAL :: BLOCK, DSKWRK, GOPARR, LTRMST, MASWRK, OUT, PANDK, PK
      REAL(KIND=dp), DIMENSION(MXGTOT) :: CD, CF, CG, CH, CI, CP, CS, EX
      REAL(KIND=dp) :: CUTOFF, QQ4, TIMLIM, TOL
      INTEGER :: IBTYP, ICOUNT, ICUT, IDAF, IEXCH, IJ, IJKL, ININTIC,   &
     &           INTLOC, IP, IPK, IPTIM, IR, IREST, IS, ISH, IST, ITOL, &
     &           IW, JSH, JST, KL, KSH, KST, LABSIX, LBUFPIC, LIT,      &
     &           LIXIC, LJT, LKT, LLT, LOCI, LOCJ, LOCK, LOCL, LSH, LST,&
     &           MASTER, MAXI, MAXJ, MAXK, MAXL, ME, MINI, MINJ, MINK,  &
     &           MINL, NANGM, NAV, NFILE, NFLTRM, NHTSHL, NIJ, NINTIC,  &
     &           NINTIX, NOPK, NORMF, NORMP, NPRINT, NPROC, NPSTRM,     &
     &           NRCTRM, NREC, NSHELL, NT, NXXIC
      INTEGER, DIMENSION(4) :: IJKLXX, NGTH
      CHARACTER(8) :: INAOFL
      INTEGER, DIMENSION(48) :: INVT
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG, KSTART, &
     &                            KTYPE
      INTEGER, DIMENSION(MXATM,48) :: MAPCTR
      INTEGER, DIMENSION(MXSH,48) :: MAPSHL
      INTEGER, DIMENSION(3) :: NORGSH, NORGSP
      REAL(KIND=dp), DIMENSION(432) :: T
      COMMON /ELGFIL/ NFILE, INAOFL
      COMMON /ELGTRM/ LTRMST, NFLTRM, NRCTRM, NPSTRM, NHTSHL
      COMMON /INT2IC/ NINTIC, ININTIC, NXXIC, LBUFPIC, LIXIC, LABSIX,   &
     &                NINTIX
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
      COMMON /PKFIL / PK, PANDK, BLOCK
      COMMON /RESTAR/ TIMLIM, IREST, NREC, INTLOC, IST, JST, KST, LST
      COMMON /SHLEXC/ NORGSH, NORGSP, IEXCH, NANGM, NGTH
      COMMON /SHLG70/ ISH, JSH, KSH, LSH, IJKLXX
      COMMON /SHLNOS/ QQ4, LIT, LJT, LKT, LLT, LOCI, LOCJ, LOCK, LOCL,  &
     &                MINI, MINJ, MINK, MINL, MAXI, MAXJ, MAXK, MAXL,   &
     &                NIJ, IJ, KL, IJKL
      COMMON /SHLT  / TOL, CUTOFF, ICOUNT, OUT
      COMMON /SYMTRY/ MAPSHL, MAPCTR, T, INVT, NT
C
      INTEGER :: L1, MAXG, MXBUF, NCUR, NINT, NMAX, NSCHWZ, NSH2, NV
      LOGICAL :: LMRSF, SCHWRZ
      REAL(KIND=dp), DIMENSION(*) :: AGDLR, ADCO1, ADCO2, ADO1V, ADO2V,
     &                               ACO12, AO21V
      REAL(KIND=dp), DIMENSION(*) :: BALL, BCO1, BCO2, BO1V, BO2V,
     &                               O21V, CO12
      REAL(KIND=dp), DIMENSION(*) :: DA, F
      REAL(KIND=dp), DIMENSION(NSH2) :: DSH, XINTS
      REAL(KIND=dp), DIMENSION(MAXG) :: GHONDO
      INTEGER, DIMENSION(L1) :: IA
C
      LOGICAL :: C1GRP, CMBDIR, DLB, GPSAVE, NPSYM, SCHSKP, SKIPA,      &
     &           SKIPB, SKIPC, SLB
      REAL(KIND=dp) :: DENMAX, Q4, TEST, TIM
      INTEGER :: I, ID, IH, II, IJIJ, IPCOUNT, IT, J0, JD, JH, JJ, K0,  &
     &           KD, KH, KK, KLKL, L0, LD, LL, LMAX, M, MINE, N4, ND,   &
     &           NEXT
      INTEGER, DIMENSION(48) :: M0, MI, MJ, MK
      REAL(KIND=dp) :: SCHWDN

      INTEGER, ALLOCATABLE :: IIB(:), KKB(:), JJB(:), LLB(:)
      REAL(KIND=dp), ALLOCATABLE :: VALB(:)

C          ----- TWO-ELECTRON INTEGRALS -----
C     THIS VERSION CAN HANDLE S,P,D,F,G AND L SHELLS
C
      TIM = ZERO
      CALL TSECND(TIM)
C
      CMBDIR= .TRUE.
C
C        THE OLD FASHIONED PARALLEL INTEGRAL TRANSFORMATIONS DO NOT
C        ALLOW THE AO INTEGRAL WORK TO BE RUN IN PARALLEL.  THE MODERN
C        DISTRIBUTED MEMORY TRANSFORMATIONS DO NOT PASS THROUGH HERE.
C
      GPSAVE = GOPARR
C
C     ----- INITIALIZATION FOR PARALLEL WORK -----
C     BOTH STATIC AND DYNAMIC LOAD BALANCING ARE IMPLEMENTED BELOW
C
      SLB = GOPARR  .AND.  IBTYP.EQ.0
      DLB = GOPARR  .AND.  IBTYP.EQ.1
      NEXT = -1
      MINE = -1
      IPCOUNT = ME - 1
C
      C1GRP = NT.EQ.1
C
      CALL BASCHK(LMAX)
                    NANGM =  4
      IF(LMAX.EQ.2) NANGM =  6
      IF(LMAX.EQ.3) NANGM = 10
      IF(LMAX.EQ.4) NANGM = 15
      IF(LMAX.EQ.5) NANGM = 21
      IF(LMAX.EQ.6) NANGM = 28
      NGTH(4) = 1
      NGTH(3) = NGTH(4) * NANGM
      NGTH(2) = NGTH(3) * NANGM
      NGTH(1) = NGTH(2) * NANGM
      IF(NOPK.EQ.0) THEN
         NORGSH(1) = 0
         NORGSH(2) = NORGSH(1) + NANGM**4
         NORGSH(3) = NORGSH(2) + NANGM**4
         NORGSP(1) = 0
         NORGSP(2) = 256
         NORGSP(3) = 512
      ELSE
         DO I=1,3
            NORGSH(I) = 0
            NORGSP(I) = 0
         END DO
      END IF
C
      NINT  = 0
      NCUR  = 0
      NSCHWZ= 0
      SCHSKP=.FALSE.
      DENMAX = ZERO
C
C        NOW WE ARE READY TO LOOP OVER ALL NSHELL**4 SHELL QUARTETS

      MXBUF = 1000000        !SIZE OF BUFFER FOR INTEGRALS
      ALLOCATE(IIB(MXBUF))
      ALLOCATE(JJB(MXBUF))
      ALLOCATE(KKB(MXBUF))
      ALLOCATE(LLB(MXBUF))
      ALLOCATE(VALB(MXBUF))
C
C     ----- I SHELL -----
C
      DO II = IST,NSHELL
C
C     ----- CHECK CPU TIME -----
C
      CALL TSECND(TIM)
      IF(TIM.GE.TIMLIM) THEN
C        NOTHING CAN BE DONE FOR IN-CORE INTEGRALS: JUST FORGET THEM
         IF(MASWRK) WRITE(IW,9030)
         CALL ABRT
      END IF
C
C         ELONGATION METHOD
      IF(LTRMST) THEN
         IF(II.EQ.NHTSHL) THEN
            NFLTRM = NFILE
            NRCTRM = NREC
            NPSTRM = ICOUNT
         END IF
      END IF
C
C     ----- PRINT INTERMEDIATE RESTART DATA -----
C
      IF(NPRINT.NE.-5  .AND.  .NOT.CMBDIR .AND. MASWRK) THEN
         IF(ICOUNT.LE.NINTIC) THEN
            WRITE(IW,9015) II,JST,KST,LST,ICOUNT
         ELSE
            WRITE(IW,9010) II,JST,KST,LST,NREC,ICOUNT-NINTIC
         END IF
      END IF
C
C     ----- SKIP I SHELL IF NOT SYMMETRY UNIQUE -----
C     THIS, AND THE SIMILAR BRANCHINGS FOR THE J, K, AND L LOOPS IS
C     WHAT GENERATES THE "PETITE" RATHER THAN "GRANDE" INTEGRAL LIST.
C
      IF(C1GRP) THEN
         MI(1)=II
      ELSE
         DO IT = 1,NT
            ID = MAPSHL(II,IT)
            IF (ID .GT. II) GO TO 920
            MI(IT) = ID
         END DO
      END IF
C
C     ----- J SHELL -----
C
      J0 = JST
      DO JJ = J0,II
      JST = 1
C
      IF(C1GRP) THEN
         MJ(1)=JJ
      ELSE
         DO IT = 1,NT
            ID = MI(IT)
            JD = MAPSHL(JJ,IT)
            MJ(IT) = JD
            IF (ID .LT. JD) THEN
               ND = ID
               ID = JD
               JD = ND
            END IF
            IF (ID .LE. II) THEN
               IF (JD .GT. JJ) GO TO 900
            ELSE
               GO TO 900
            END IF
         END DO
      END IF
C
C     ----- GO PARALLEL! -----
C
      IF (DLB) THEN
         MINE = MINE + 1
         IF (MINE.GT.NEXT) CALL DDI_DLBNEXT(NEXT)
         IF (NEXT.NE.MINE) CYCLE
      END IF
C
C     ----- K SHELL -----
C
      K0 = KST
      DO KK = K0,JJ
         KST = 1
C
         IF(C1GRP) THEN
            MK(1)=KK
         ELSE
            DO IT = 1,NT
               ID = MI(IT)
               JD = MJ(IT)
               KD = MAPSHL(KK,IT)
               MK(IT) = KD
  240          IF (ID .LT. JD) THEN
                  ND = ID
                  ID = JD
                  JD = ND
               END IF
               IF (JD .LT. KD) THEN
                  ND = JD
                  JD = KD
                  KD = ND
                  GO TO 240
               ELSE IF (ID .LE. II) THEN
                  IF (JD .LE. JJ) THEN
                     IF (KD .GT. KK) GO TO 880
                  ELSE
                     GO TO 880
                  END IF
               ELSE
                  GO TO 880
               END IF
            END DO
         END IF
C
C     ----- GO PARALLEL! -----
C
      IF(SLB) THEN
         IPCOUNT = IPCOUNT + 1
         IF (MOD(IPCOUNT,NPROC).NE.0) CYCLE
      END IF
C
C     ----- L SHELL ----
C
      L0 = LST
      DO LL = L0,KK
         LST = 1
C
         IF(C1GRP) THEN
            M0(1)=1
            N4=1
         ELSE
            N4 = 0
            DO IT = 1,NT
               ID = MI(IT)
               JD = MJ(IT)
               KD = MK(IT)
               LD = MAPSHL(LL,IT)
  380          IF (ID .LT. JD) THEN
                  ND = ID
                  ID = JD
                  JD = ND
               END IF
  400          IF (JD .LT. KD) THEN
                  ND = JD
                  JD = KD
                  KD = ND
                  GO TO 380
               ELSE IF (KD .LT. LD) THEN
                  ND = KD
                  KD = LD
                  LD = ND
                  GO TO 400
               ELSE IF (ID .LE. II) THEN
                  IF (JD .LE. JJ) THEN
                     IF (KD .LE. KK) THEN
                        IF (LD .LE. LL) THEN
                           N4 = N4+1
                           M0(N4) = IT
                        ELSE
                           CYCLE
                        END IF
                     ELSE
                        CYCLE
                     END IF
                  ELSE
                     CYCLE
                  END IF
               ELSE
                  CYCLE
               END IF
            END DO
         END IF
C
C         THE LOOP STRUCTURE IN THIS ROUTINE IS DESIGNED TO FACILITATE
C         SUPERMATRIX CONSTRUCTION BY HAVING UP TO THREE "EXCHANGED"
C         QUARTETS AVAILABLE AT ONCE.  THE LOOP STRUCTURE TO GENERATE
C         A MORE NORMAL CANONICAL ORDERING OF THE QUARTETS HITS THE
C         SAME QUARTETS IN A SLIGHTLY DIFFERENT ORDER, BUT BOTH LOOPS
C         WILL DO EXACTLY THE SAME QUARTETS.
C
C             CANONICAL                      SUPERMATRIX
C         DO ISH=1,NSHELL                 DO II=1,NSHELL
C           DO JSH=1,ISH                    DO JJ=1,II
C             IJSH = IA(ISH)+JSH
C             DO KSH=1,ISH                    DO KK=1,JJ
C               DO LSH=1,KSH                    DO LL=1,KK
C                 KLSH=IA(KSH)+LSH
C                 IF(IJSH.LT.KLSH),               [II JJ|KK LL],
C                    CYCLE KSH LOOP               [II KK|JJ LL],
C                 [ISH JSH|KSH LSH]               [II LL|JJ KK]
C               ENDDO                           ENDDO
C             ENDDO                           ENDDO
C           ENDDO                           ENDDO
C         ENDDO                           ENDDO
C
C     ----- CHECK FOR REDUNDANIES BETWEEN THE 3 COMBINATIONS -----
C            (II,JJ//KK,LL), (II,KK//JJ,LL), (II,LL//JJ,KK)
C
      SKIPA =  JJ.EQ.KK
      SKIPB = (II.EQ.KK) .OR. (JJ.EQ.LL)
      SKIPC = (II.EQ.JJ) .OR. (KK.EQ.LL)
      NPSYM = .FALSE.
      IF(.NOT.(SKIPA .OR. SKIPB .OR. SKIPC)) THEN
         NPSYM = .TRUE.
         DO M = 1,N4
            IT = M0(M)
            IH = MI(IT)
            JH = MJ(IT)
            IF(JH.LE.IH) THEN
               ID = IH
               JD = JH
            ELSE
               ID = JH
               JD = IH
            END IF
            IF(.NOT.SKIPA)
     *         SKIPA = (ID.EQ.II .AND. JD.EQ.KK) .OR.
     *                 (ID.EQ.JJ .AND. JD.EQ.LL)
            IF(.NOT.SKIPB)
     *         SKIPB = (ID.EQ.II .AND. JD.EQ.LL) .OR.
     *                 (ID.EQ.JJ .AND. JD.EQ.KK)
            IF(SKIPA .AND. SKIPB) THEN
               SKIPC = .TRUE.
               GO TO 720
            ELSE
               KH = MK(IT)
               IF(KH.LE.IH) THEN
                  ID = IH
                  KD = KH
               ELSE
                  ID = KH
                  KD = IH
               END IF
               IF(.NOT.SKIPC)
     *            SKIPC = (ID.EQ.II .AND. KD.EQ.LL) .OR.
     *                    (ID.EQ.JJ .AND. KD.EQ.KK)
               IF(SKIPA .AND. SKIPC) THEN
                  SKIPB = .TRUE.
                  GO TO 720
               ELSE IF(SKIPB .AND. SKIPC) THEN
                  SKIPA = .TRUE.
                  GO TO 720
               END IF
            END IF
         END DO
      END IF
C
C        GENERATE SYMMETRY FACTOR -Q4- FOR THIS QUARTET IN PETITE LIST
C
  720 Q4 = NT
      Q4 = Q4 / N4
C
C     ----- (II,JJ//KK,LL) -----
C
      IEXCH = 1
      ISH = II
      JSH = JJ
      KSH = KK
      LSH = LL
      QQ4 = Q4
      IF(SKIPA .AND. NPSYM) QQ4 = QQ4+Q4
      IF(SKIPB .AND. NPSYM) QQ4 = QQ4+Q4
      GO TO 780
C
C     ----- (II,LL//JJ,KK) -----
C
  760 IF (SKIPB .OR. SKIPC) CYCLE
      IEXCH = 3
      ISH = II
      JSH = LL
      KSH = JJ
      LSH = KK
      QQ4 = Q4
C
C        ----- COMPUTE TWO-ELECTRON INTEGRALS ----
C
  780 IF(PK .AND. IEXCH.EQ.1) CALL ZPKOUT(ISH,JSH,KSH,LSH,GHONDO,       &
     *                                    SKIPA,SKIPB,SKIPC,NPSYM)
C
C     APPLY THE SCHWARZ INEQUALITY TO SCREEN OUT SMALL INTEGRALS,
C      (II,JJ//KK,LL) .LE.  SQRT( (II,JJ//II,JJ)*(KK,LL//KK,LL) )
C     SEE, FOR EXAMPLE, J.L.WHITTEN, J.CHEM.PHYS. 58,4496-4501(1973)
C
      IF(SCHWRZ) THEN
         IJIJ = (ISH*ISH-ISH)/2 + JSH
         KLKL = (KSH*KSH-KSH)/2 + LSH
         TEST = QQ4*XINTS(IJIJ)*XINTS(KLKL)
         DENMAX = SCHWDN(DSH,ISH,JSH,KSH,LSH,IA)
c        write(6,*) "denmax =",denmax,test
         TEST = TEST*DENMAX
         SCHSKP = TEST.LT.CUTOFF
         IF(SCHSKP) NSCHWZ = NSCHWZ + 1
      END IF
C
C        ----- ELECTRON REPULSION INTEGRAL CALCULATION -----
C     THIS MAY USE ROTATED AXIS, ERIC, OR RYS QUADRATURE METHODS
C

      CALL SHELLQUART(ISH,JSH,KSH,LSH,GHONDO)
C
      CALL MRSFADDTOBUF(NCUR,IIB,KKB,JJB,LLB,VALB,GHONDO,NINT,L1,NV,
     &                  MXBUF,BO2V,BO1V,BCO1,BCO2,BALL,
     &                  O21V,CO12,
     &                  AGDLR,AO21V,ACO12,
     &                  ADO2V,ADO1V,ADCO1,ADCO2)

      IF (IEXCH .EQ. 1) THEN
C
C     ----- (II,KK//JJ,LL) -----
C
         IF (SKIPA) GO TO 760
         IEXCH = 2
         ISH = II
         JSH = KK
         KSH = JJ
         LSH = LL
         QQ4 = Q4
         IF (SKIPC .AND. NPSYM) QQ4 = QQ4+Q4
         GO TO 780
      ELSE IF (IEXCH .EQ. 2) THEN
         GO TO 760
      END IF
C
C     ----- END OF SHELL LOOPS -----
C
      END DO
  880 END DO
  900 END DO
  920 END DO
      CALL MRSFFLUSHBUF(NCUR,IIB,KKB,JJB,LLB,VALB,
     &                  BO2V,BO1V,BCO1,BCO2,BALL,
     &                  O21V,CO12,AGDLR,AO21V,ACO12,
     &                  ADO2V,ADO1V,ADCO1,ADCO2,L1,NV)
C
      IF(DLB) CALL DDI_DLBRESET
C
      GOPARR = GPSAVE
      RETURN
C
 9010 FORMAT(1X,'II,JST,KST,LST =',4I3,' NREC =',I10,' INTLOC =',I5)
 9015 FORMAT(1X,'II,JST,KST,LST =',4I3,' IN CORE, INTLOC =',I12)
 9030 FORMAT(//1X,'*** THIS JOB HAS EXHAUSTED ITS CPU TIME ***'/
     *         1X,'     (WHILE COMPUTING 2E- INTEGRALS)'///)
      END SUBROUTINE MRSFTWOEI
C
C*MODULE MRSF   *DECK MRSFADDTOBUF
C>
C>    @brief  clone of SFDIRFCK for MRSF-TDDFT
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details MRSFDIRFCK was splitted into MRSFADDTOBUF and MRSFFLUSHBUF.
C>             MRSFTWOEI is initial subroutine for this.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE MRSFADDTOBUF(NCUR,IIB,KKB,JJB,LLB,VALB,
     &                        GHONDO,NINT,NBF,NV,MXBUF,
     &                        BO2V,BO1V,BCO1,BCO2,BALL,
     &                        O21V,CO12,AGDLR,AO21V,ACO12,
     &                        ADO2V,ADO1V,ADCO1,ADCO2)
C
      USE mx_limits, ONLY: mxsh, mxgtot
      USE prec, ONLY: dp
      USE constants, ONLY: half
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(MXGTOT) :: CD, CF, CG, CH, CI, CP, CS, EX
      REAL(KIND=dp) :: CUTOFF, QQ4, TOL
      LOGICAL :: IANDJ, KANDL, OUT, SAME
      INTEGER :: ICOUNT, IJ, IJKL, ISH, JSH, KL, KSH, LIT, LJT, LKT,    &
     &           LLT, LOCI, LOCJ, LOCK, LOCL, LSH, LSTRI, LSTRJ, LSTRK, &
     &           LSTRL, MAXI, MAXJ, MAXK, MAXL, MINI, MINJ, MINK, MINL, &
     &           NIJX, NSHELL
      INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG, KSTART, &
     &                            KTYPE
      COMMON /ERIOUT/ ISH, JSH, KSH, LSH, LSTRI, LSTRJ, LSTRK, LSTRL
      COMMON /MISC  / IANDJ, KANDL, SAME
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
      COMMON /SHLNOS/ QQ4, LIT, LJT, LKT, LLT, LOCI, LOCJ, LOCK, LOCL,  &
     &                MINI, MINJ, MINK, MINL, MAXI, MAXJ, MAXK, MAXL,   &
     &                NIJX, IJ, KL, IJKL
      COMMON /SHLT  / TOL, CUTOFF, ICOUNT, OUT
C
      INTEGER :: NINT
      REAL(KIND=dp), DIMENSION(*) :: GHONDO
C
      REAL(KIND=dp) :: CUTINT, VAL
      INTEGER :: I, II, IJKL_ERIC, IJK_ERIC, IJ_ERIC, IMO, I_ERIC, J,   &
     &           JJ, JMAX, K, KK, L, LL, LMAX, N
      REAL(KIND=dp), DIMENSION(*) :: AGDLR, ADCO1, ADCO2, ADO1V, ADO2V, &
     &                               ACO12, AO21V
      REAL(KIND=dp), DIMENSION(*) :: BALL, BCO1, BCO2, BO1V, BO2V,      &
     &                               O21V, CO12

      INTEGER :: NBF, NCUR, NV, MXBUF
      INTEGER, DIMENSION(*) :: IIB, KKB, JJB, LLB
      REAL(KIND=dp), DIMENSION(*) :: VALB
C
C-NEXT STATEMENT IS FOR VARIOUS IBM XLF 3.X AND 5.X COMPILERS-
C
      INTEGER, SAVE :: IJN, KLN
C
C     --- FORM A BATCH OF SQUARE NON-SYMMETRIC FOCK-LIKE MATRICES
C     DIRECTLY FROM INTEGRALS ---
C
      CUTINT = CUTOFF
C
      SAME  = ISH .EQ. KSH .AND. JSH .EQ. LSH
      IANDJ = ISH .EQ. JSH
      KANDL = KSH .EQ. LSH
C
      MINI = KMIN(ISH)
      MINJ = KMIN(JSH)
      MINK = KMIN(KSH)
      MINL = KMIN(LSH)
      MAXI = KMAX(ISH)
      MAXJ = KMAX(JSH)
      MAXK = KMAX(KSH)
      MAXL = KMAX(LSH)
      LOCI = KLOC(ISH)-MINI
      LOCJ = KLOC(JSH)-MINJ
      LOCK = KLOC(KSH)-MINK
      LOCL = KLOC(LSH)-MINL
C
      IJN = 0
      JMAX = MAXJ
      DO I = MINI,MAXI
         I_ERIC = (I-MINI)*LSTRI + 1
         IF (IANDJ) JMAX = I
         DO J = MINJ,JMAX
            IJ_ERIC = (J-MINJ)*LSTRJ + I_ERIC
            IJN = IJN+1
            LMAX = MAXL
            KLN = 0
            DO K =  MINK,MAXK
               IJK_ERIC = (K-MINK)*LSTRK + IJ_ERIC
               IF (KANDL) LMAX = K
               DO L = MINL,LMAX
                  KLN = KLN+1
                  IF(SAME .AND. KLN.GT.IJN) GO TO 340
C
                  IJKL_ERIC = (L-MINL)*LSTRL + IJK_ERIC
                  VAL = GHONDO( IJKL_ERIC )
C
                  IF(ABS(VAL).LT.CUTINT) CYCLE
                  NINT = NINT + 1
C
                  II = LOCI+I
                  JJ = LOCJ+J
                  KK = LOCK+K
                  LL = LOCL+L
                  IF (II .LT. JJ) THEN
                     N = II
                     II = JJ
                     JJ = N
                  END IF
                  IF (KK .LT. LL) THEN
                     N = KK
                     KK = LL
                     LL = N
                  END IF
                  IF (II .LE. KK) THEN
                     IF (JJ .GE. LL) GO TO 180
                  ELSE
                     GO TO 180
                  END IF
                  N = II
                  II = KK
                  KK = N
                  N = JJ
                  JJ = LL
                  LL = N
C
  180             IF(II.EQ.JJ) VAL = VAL*HALF
                  IF(KK.EQ.LL) VAL = VAL*HALF
                  IF(II.EQ.KK  .AND.  JJ.EQ.LL) VAL = VAL*HALF
                  NCUR = NCUR+1
                  IIB(NCUR) = II
                  JJB(NCUR) = JJ
                  KKB(NCUR) = KK
                  LLB(NCUR) = LL
                  VALB(NCUR) = VAL
                  IF (NCUR == MXBUF) THEN
                   CALL MRSFFLUSHBUF(NCUR,IIB,KKB,JJB,LLB,VALB,
     &                               BO2V,BO1V,BCO1,BCO2,BALL,
     &                               O21V,CO12,
     &                               AGDLR,AO21V,ACO12,
     &                               ADO2V,ADO1V,ADCO1,ADCO2,NBF,NV)
                   NCUR = 0
                  END IF
               END DO
            END DO
 340     END DO
      END DO
C
      RETURN
      END SUBROUTINE MRSFADDTOBUF
C
C*MODULE MRSF   *DECK MRSFFLUSHBUF
C>
C>    @brief  clone of SFDIRFCK for MRSF-TDDFT
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details MRSFDIRFCK was splitted into MRSFADDTOBUF and MRSFFLUSHBUF.
C>             MRSFTWOEI is initial subroutine for this.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE MRSFFLUSHBUF(NCUR,IIB,KKB,JJB,LLB,VALB,
     &                        BO2V,BO1V,BCO1,BCO2,BALL,
     &                        O21V,CO12,
     &                        AGDLR,AO21V,ACO12,
     &                        ADO2V,ADO1V,ADCO1,ADCO2,NBF,NV)

      USE prec, ONLY: dp
      IMPLICIT NONE

      INTEGER :: NBF, NCUR, NV
      REAL(KIND=dp), DIMENSION(NBF,NBF,*) :: ADCO1, ADCO2, ADO1V, ADO2V,
     &                               AGDLR, ACO12, AO21V
      REAL(KIND=dp), DIMENSION(NBF,NBF,*) :: BCO1, BCO2, BO1V, BO2V,
     &                               BALL, O21V, CO12

      INTEGER, DIMENSION(*) :: IIB, KKB, JJB, LLB
      REAL(KIND=dp), DIMENSION(*) :: VALB

      INTEGER :: II, IMO, JJ, KK, LL, M
      REAL(KIND=dp) :: VAL

      DO IMO = 1,NV
         DO M = 1,NCUR
            II = IIB(M)
            KK = KKB(M)
            JJ = JJB(M)
            LL = LLB(M)
            VAL = VALB(M)
C
C   --- LOOP OVER BATCH OF MOS AND MAKE CONTRIB. TO
C        --- FOCK-LIKE MATRIX ---
C
            agdlr(ii,kk,imo) = agdlr(ii,kk,imo) - val*ball(jj,ll,imo)
            agdlr(kk,ii,imo) = agdlr(kk,ii,imo) - val*ball(ll,jj,imo)
            agdlr(ii,ll,imo) = agdlr(ii,ll,imo) - val*ball(jj,kk,imo)
            agdlr(ll,ii,imo) = agdlr(ll,ii,imo) - val*ball(kk,jj,imo)
            agdlr(jj,kk,imo) = agdlr(jj,kk,imo) - val*ball(ii,ll,imo)
            agdlr(kk,jj,imo) = agdlr(kk,jj,imo) - val*ball(ll,ii,imo)
            agdlr(jj,ll,imo) = agdlr(jj,ll,imo) - val*ball(ii,kk,imo)
            agdlr(ll,jj,imo) = agdlr(ll,jj,imo) - val*ball(kk,ii,imo)

            ao21v(ii,ll,imo) = ao21v(ii,ll,imo) - val*o21v(kk,jj,imo)! (ij|lk)
            ao21v(kk,jj,imo) = ao21v(kk,jj,imo) - val*o21v(ii,ll,imo)! (kl|ji)
            ao21v(ii,kk,imo) = ao21v(ii,kk,imo) - val*o21v(ll,jj,imo)! (ij|kl)
            ao21v(ll,jj,imo) = ao21v(ll,jj,imo) - val*o21v(ii,kk,imo)! (lk|ji)
            ao21v(jj,ll,imo) = ao21v(jj,ll,imo) - val*o21v(kk,ii,imo)! (ji|lk)
            ao21v(kk,ii,imo) = ao21v(kk,ii,imo) - val*o21v(jj,ll,imo)! (kl|ij)
            ao21v(jj,kk,imo) = ao21v(jj,kk,imo) - val*o21v(ll,ii,imo)! (ji|kl)
            ao21v(ll,ii,imo) = ao21v(ll,ii,imo) - val*o21v(jj,kk,imo)! (lk|ij)

            aco12(ii,ll,imo) = aco12(ii,ll,imo) - val*co12(kk,jj,imo)! (ij|lk)
            aco12(kk,jj,imo) = aco12(kk,jj,imo) - val*co12(ii,ll,imo)! (kl|ji)
            aco12(ii,kk,imo) = aco12(ii,kk,imo) - val*co12(ll,jj,imo)! (ij|kl)
            aco12(ll,jj,imo) = aco12(ll,jj,imo) - val*co12(ii,kk,imo)! (lk|ji)
            aco12(jj,ll,imo) = aco12(jj,ll,imo) - val*co12(kk,ii,imo)! (ji|lk)
            aco12(kk,ii,imo) = aco12(kk,ii,imo) - val*co12(jj,ll,imo)! (kl|ij)
            aco12(jj,kk,imo) = aco12(jj,kk,imo) - val*co12(ll,ii,imo)! (ji|kl)
            aco12(ll,ii,imo) = aco12(ll,ii,imo) - val*co12(jj,kk,imo)! (lk|ij)

            ado2v(ii,kk,imo) = ado2v(ii,kk,imo) - val*bo2v(jj,ll,imo)! (ij|lk)
            ado2v(kk,ii,imo) = ado2v(kk,ii,imo) - val*bo2v(ll,jj,imo)! (kl|ji)
            ado2v(ii,ll,imo) = ado2v(ii,ll,imo) - val*bo2v(jj,kk,imo)! (ij|kl)
            ado2v(ll,ii,imo) = ado2v(ll,ii,imo) - val*bo2v(kk,jj,imo)! (lk|ji)
            ado2v(jj,kk,imo) = ado2v(jj,kk,imo) - val*bo2v(ii,ll,imo)! (ji|lk)
            ado2v(kk,jj,imo) = ado2v(kk,jj,imo) - val*bo2v(ll,ii,imo)! (kl|ij)
            ado2v(jj,ll,imo) = ado2v(jj,ll,imo) - val*bo2v(ii,kk,imo)! (ji|kl)
            ado2v(ll,jj,imo) = ado2v(ll,jj,imo) - val*bo2v(kk,ii,imo)! (lk|ij)
            ado2v(ii,jj,imo) = ado2v(ii,jj,imo) + val*bo2v(kk,ll,imo)! (ij|lk)
            ado2v(kk,ll,imo) = ado2v(kk,ll,imo) + val*bo2v(ii,jj,imo)! (kl|ji)
            ado2v(ii,jj,imo) = ado2v(ii,jj,imo) + val*bo2v(ll,kk,imo)! (ij|kl)
            ado2v(ll,kk,imo) = ado2v(ll,kk,imo) + val*bo2v(ii,jj,imo)! (lk|ji)
            ado2v(jj,ii,imo) = ado2v(jj,ii,imo) + val*bo2v(kk,ll,imo)! (ji|lk)
            ado2v(kk,ll,imo) = ado2v(kk,ll,imo) + val*bo2v(jj,ii,imo)! (kl|ij)
            ado2v(jj,ii,imo) = ado2v(jj,ii,imo) + val*bo2v(ll,kk,imo)! (ji|kl)
            ado2v(ll,kk,imo) = ado2v(ll,kk,imo) + val*bo2v(jj,ii,imo)! (lk|ij)

            ado1v(ii,kk,imo) = ado1v(ii,kk,imo) - val*bo1v(jj,ll,imo)
            ado1v(kk,ii,imo) = ado1v(kk,ii,imo) - val*bo1v(ll,jj,imo)
            ado1v(ii,ll,imo) = ado1v(ii,ll,imo) - val*bo1v(jj,kk,imo)
            ado1v(ll,ii,imo) = ado1v(ll,ii,imo) - val*bo1v(kk,jj,imo)
            ado1v(jj,kk,imo) = ado1v(jj,kk,imo) - val*bo1v(ii,ll,imo)
            ado1v(kk,jj,imo) = ado1v(kk,jj,imo) - val*bo1v(ll,ii,imo)
            ado1v(jj,ll,imo) = ado1v(jj,ll,imo) - val*bo1v(ii,kk,imo)
            ado1v(ll,jj,imo) = ado1v(ll,jj,imo) - val*bo1v(kk,ii,imo)
            ado1v(ii,jj,imo) = ado1v(ii,jj,imo) + val*bo1v(kk,ll,imo)! (ij|lk)
            ado1v(kk,ll,imo) = ado1v(kk,ll,imo) + val*bo1v(ii,jj,imo)! (kl|ji)
            ado1v(ii,jj,imo) = ado1v(ii,jj,imo) + val*bo1v(ll,kk,imo)! (ij|kl)
            ado1v(ll,kk,imo) = ado1v(ll,kk,imo) + val*bo1v(ii,jj,imo)! (lk|ji)
            ado1v(jj,ii,imo) = ado1v(jj,ii,imo) + val*bo1v(kk,ll,imo)! (ji|lk)
            ado1v(kk,ll,imo) = ado1v(kk,ll,imo) + val*bo1v(jj,ii,imo)! (kl|ij)
            ado1v(jj,ii,imo) = ado1v(jj,ii,imo) + val*bo1v(ll,kk,imo)! (ji|kl)
            ado1v(ll,kk,imo) = ado1v(ll,kk,imo) + val*bo1v(jj,ii,imo)! (lk|ij)

            adco1(ii,kk,imo) = adco1(ii,kk,imo) - val*bco1(jj,ll,imo)
            adco1(kk,ii,imo) = adco1(kk,ii,imo) - val*bco1(ll,jj,imo)
            adco1(ii,ll,imo) = adco1(ii,ll,imo) - val*bco1(jj,kk,imo)
            adco1(ll,ii,imo) = adco1(ll,ii,imo) - val*bco1(kk,jj,imo)
            adco1(jj,kk,imo) = adco1(jj,kk,imo) - val*bco1(ii,ll,imo)
            adco1(kk,jj,imo) = adco1(kk,jj,imo) - val*bco1(ll,ii,imo)
            adco1(jj,ll,imo) = adco1(jj,ll,imo) - val*bco1(ii,kk,imo)
            adco1(ll,jj,imo) = adco1(ll,jj,imo) - val*bco1(kk,ii,imo)
            adco1(ii,jj,imo) = adco1(ii,jj,imo) + val*bco1(kk,ll,imo)! (ij|lk)
            adco1(kk,ll,imo) = adco1(kk,ll,imo) + val*bco1(ii,jj,imo)! (kl|ji)
            adco1(ii,jj,imo) = adco1(ii,jj,imo) + val*bco1(ll,kk,imo)! (ij|kl)
            adco1(ll,kk,imo) = adco1(ll,kk,imo) + val*bco1(ii,jj,imo)! (lk|ji)
            adco1(jj,ii,imo) = adco1(jj,ii,imo) + val*bco1(kk,ll,imo)! (ji|lk)
            adco1(kk,ll,imo) = adco1(kk,ll,imo) + val*bco1(jj,ii,imo)! (kl|ij)
            adco1(jj,ii,imo) = adco1(jj,ii,imo) + val*bco1(ll,kk,imo)! (ji|kl)
            adco1(ll,kk,imo) = adco1(ll,kk,imo) + val*bco1(jj,ii,imo)! (lk|ij)

            adco2(ii,kk,imo) = adco2(ii,kk,imo) - val*bco2(jj,ll,imo)
            adco2(kk,ii,imo) = adco2(kk,ii,imo) - val*bco2(ll,jj,imo)
            adco2(ii,ll,imo) = adco2(ii,ll,imo) - val*bco2(jj,kk,imo)
            adco2(ll,ii,imo) = adco2(ll,ii,imo) - val*bco2(kk,jj,imo)
            adco2(jj,kk,imo) = adco2(jj,kk,imo) - val*bco2(ii,ll,imo)
            adco2(kk,jj,imo) = adco2(kk,jj,imo) - val*bco2(ll,ii,imo)
            adco2(jj,ll,imo) = adco2(jj,ll,imo) - val*bco2(ii,kk,imo)
            adco2(ll,jj,imo) = adco2(ll,jj,imo) - val*bco2(kk,ii,imo)
            adco2(ii,jj,imo) = adco2(ii,jj,imo) + val*bco2(kk,ll,imo)! (ij|lk)
            adco2(kk,ll,imo) = adco2(kk,ll,imo) + val*bco2(ii,jj,imo)! (kl|ji)
            adco2(ii,jj,imo) = adco2(ii,jj,imo) + val*bco2(ll,kk,imo)! (ij|kl)
            adco2(ll,kk,imo) = adco2(ll,kk,imo) + val*bco2(ii,jj,imo)! (lk|ji)
            adco2(jj,ii,imo) = adco2(jj,ii,imo) + val*bco2(kk,ll,imo)! (ji|lk)
            adco2(kk,ll,imo) = adco2(kk,ll,imo) + val*bco2(jj,ii,imo)! (kl|ij)
            adco2(jj,ii,imo) = adco2(jj,ii,imo) + val*bco2(ll,kk,imo)! (ji|kl)
            adco2(ll,kk,imo) = adco2(ll,kk,imo) + val*bco2(jj,ii,imo)! (lk|ij)
         END DO
      END DO
C
      RETURN
      END SUBROUTINE MRSFFLUSHBUF
C*MODULE MRSF   *DECK MRSFESUM
C>
C>    @brief  clone of SFESUM
C>            for singlet and triplet UKS/MRSF
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE MRSFESUM(EA,EB,PMO,Z,LX,L7,NOCA,NOCB,IVEC)
C
      USE constants, ONLY: two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IVEC, L7, LX, NOCA, NOCB
      REAL(KIND=dp), DIMENSION(*) :: EA, EB
      REAL(KIND=dp), DIMENSION(L7,*) :: PMO, Z
C
      INTEGER :: I, IJ, J
C
C     ----- ADD (EA-EI)*ZAI -----
C
      IJ = 0
      DO J=NOCB+1,LX
         DO I=1,NOCA
            if(i.eq.noca-1 .and. j.eq.noca-1) then
               IJ = IJ + 1
               PMO(IJ,IVEC) = PMO(IJ,IVEC) +
     *         (eb(noca-1)+eb(noca)-ea(noca-1)-ea(noca))*Z(IJ,IVEC)/two

            else if(i.eq.noca   .and. j.eq.noca ) then
            else
               IJ = IJ + 1
               PMO(IJ,IVEC) = PMO(IJ,IVEC) + (EB(J)-EA(I))*Z(IJ,IVEC)
            endif
         END DO
      END DO
c
      RETURN
      END SUBROUTINE MRSFESUM
C
c*MODULE MRSF   *deck mrsfmntoia
C>
C>    @brief  clone of SFMNTOIA
C>            for singlet and triplet MRSF
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Terms co2v, co1v, cco1 and cco2 was replaced
C>             by o21v, co12. mrsfcbc is initial subroutine of this.
C>             Optimized cycles (3 lvl -> 2 lvl)
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE MRSFMNTOIA(agdlr,ao21v,aco12,
     &                      ado2v,ado1v,adco1,adco2,
     &                      pmo,va,vb,scr,wrk,l1,lx,l7,noca,
     &                      nocb,ivec,lmrsf,mrsfs,mrsft)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, one, two
      IMPLICIT NONE
C
      INTEGER :: IVEC, L1, L7, LX, NOCA, NOCB
      LOGICAL :: LMRSF, MRSFS, MRSFT
      REAL(KIND=dp), DIMENSION(L1,*) :: ADCO1, ADCO2, ADO1V, ADO2V,
     &                            AGDLR, ACO12, AO21V
      REAL(KIND=dp), DIMENSION(L1) :: TMP1, TMP2
      REAL(KIND=dp), DIMENSION(l7,*) :: PMO
      REAL(KIND=dp), DIMENSION(lx,*) :: SCR, WRK
      REAL(KIND=dp), DIMENSION(l1,*) :: VA, VB
C
      REAL(KIND=dp) :: DUMN
      INTEGER :: I, IJ, IJD, IJG, IJLR2, J, M, N, LR1, LR2
C
C
      CALL DGEMM('T','N',LX,L1,L1,ONE,VA,L1,agdlr,L1,ZERO,WRK,LX)
      CALL DGEMM('N','N',LX,LX,L1,ONE,WRK,LX,VB,L1,ZERO,scr,LX)

      CALL DGEMM('T','N',LX,L1,L1,ONE,VA,L1,agdlr,L1,ZERO,WRK,LX)
      CALL DGEMM('N','N',LX,LX,L1,ONE,WRK,LX,VB,L1,ZERO,scr,LX)
c 1
      LR1 = NOCA-1
      LR2 = NOCA
c
c     ----- (m,n) to (i+,n) -----
c
      ij = 0
      if (mrsfs) then
         do j = nocb+1,lx
         do i = 1,noca
            if(i.eq.lr2 .and. j.eq.lr2) cycle
            if(i.eq.lr1 .and. j.eq.lr1) then
               ij = ij+1
               pmo(ij,ivec) = (scr(lr1,lr1)-scr(lr2,lr2))/sqrt(two)
            else
               ij = ij+1
               pmo(ij,ivec) = scr(i,j)
            endif
         end do
         end do
      elseif (mrsft) then
         do j = nocb+1,lx
         do i = 1,noca
            if(i.eq.lr2 .and. j.eq.lr2) cycle
            if(i.eq.lr2 .and. j.eq.lr1) cycle
            if(i.eq.lr1 .and. j.eq.lr2) cycle
            if(i.eq.lr1 .and. j.eq.lr1) then
               ij = ij+1
               pmo(ij,ivec) = (scr(lr1,lr1)+scr(lr2,lr2))/sqrt(two)
            else
               ij = ij+1
               pmo(ij,ivec) = scr(i,j)
            endif
         end do
         end do
      endif

      ijg   = (noca-1-nocb-1)*noca+noca
      ijd   = (noca  -nocb-1)*noca+noca-1
      ijlr2 = (noca  -nocb-1)*noca+noca
c 3
      tmp1 = zero
      tmp2 = zero
      call dgemv('n',l1,l1,one,ado1v,l1,vb(:,lr2),1,one,tmp1,1)
      call dgemv('n',l1,l1,one,aco12,l1,vb(:,lr1),1,one,tmp2,1)
      tmp1 = tmp1+tmp2
      do i = 1,noca-2
         ij = (noca-nocb-1)*noca+i
         if(ij.gt.ijlr2) ij = ij-1
         if(mrsft .and. ij.gt.ijg) ij = ij-1
         if(mrsft .and. ij.gt.ijd) ij = ij-1
         dumn = zero
         do n = 1,l1
             dumn = dumn+va(n,i)*tmp1(n)
         end do
         pmo(ij,ivec) = pmo(ij,ivec)+dumn
      end do
c 4
      tmp1 = zero
      tmp2 = zero
      call dgemv('n',l1,l1,one,ado2v,l1,vb(:,lr1),1,one,tmp1,1)
      call dgemv('n',l1,l1,one,aco12,l1,vb(:,lr2),1,one,tmp2,1)
      tmp1 = tmp1-tmp2
      do i = 1,noca-2
         ij = (lr1-nocb-1)*noca+i
         if(ij.gt.ijlr2) ij = ij-1
         if(mrsft .and. ij.gt.ijg) ij = ij-1
         if(mrsft .and. ij.gt.ijd) ij = ij-1
         dumn = zero
         do n = 1,l1
            dumn = dumn+va(n,i)*tmp1(n)
         end do
         pmo(ij,ivec) = pmo(ij,ivec)+dumn
      end do
c 5
      tmp1 = zero
      tmp2 = zero
      call dgemv('t',l1,l1,one,adco2,l1,va(:,lr1),1,one,tmp1,1)
      call dgemv('t',l1,l1,one,ao21v,l1,va(:,lr2),1,one,tmp2,1)
      tmp1 = tmp1+tmp2
      do j = noca+1,lx
         ij = (j-nocb-1)*noca+lr1
         if(ij.gt.ijlr2) ij = ij-1
         if(mrsft .and. ij.gt.ijg) ij = ij-1
         if(mrsft .and. ij.gt.ijd) ij = ij-1
         dumn = zero
         do m = 1,l1
            dumn = dumn+vb(m,j)*tmp1(m)
         end do
         pmo(ij,ivec) = pmo(ij,ivec)+dumn
      end do
c 6
      tmp1 = zero
      tmp2 = zero
      call dgemv('t',l1,l1,one,adco1,l1,va(:,lr2),1,one,tmp1,1)
      call dgemv('t',l1,l1,one,ao21v,l1,va(:,lr1),1,one,tmp2,1)
      tmp1 = tmp1-tmp2
      do j = noca+1,lx
         ij = (j-nocb-1)*noca+lr2
         if(ij.gt.ijlr2) ij = ij-1
         if(mrsft .and. ij.gt.ijg) ij = ij-1
         if(mrsft .and. ij.gt.ijd) ij = ij-1
         dumn = zero
         do m = 1,l1
            dumn = dumn+vb(m,j)*tmp1(m)
         end do
         pmo(ij,ivec) = pmo(ij,ivec)+dumn
      end do

      return
      end subroutine mrsfmntoia
C*MODULE MRSF   *DECK MRSFQMNTOIA
C>
C>    @brief  clone of SFMNTOIA
C>            for quintet MRSF
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE MRSFQMNTOIA(PAO,PMO,VA,VB,SCR,L1,LX,L7,NOCA,NOCB,IVEC)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, one
      IMPLICIT NONE
C
      INTEGER :: IVEC, L1, L7, LX, NOCA, NOCB
      REAL(KIND=dp), DIMENSION(L1,*) :: PAO
      REAL(KIND=dp), DIMENSION(L7,*) :: PMO
      REAL(KIND=dp), DIMENSION(NOCB,*) :: SCR
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA, VB
C
      REAL(KIND=dp) :: DUMN
      INTEGER :: I, IJ, J, N
C
C     ----- (M,N) TO (I+,N) -----
C
      CALL DGEMM('T','N',NOCB,L1,L1,ONE,VB,L1,PAO,L1,ZERO,SCR,NOCB)
C
C     ----- (I+,N) TO (I+A-) -----
C
      IJ = 0
      DO J=NOCA+1,LX
         DO I=1,NOCB
            IJ = IJ + 1
            DUMN = ZERO
            DO N=1,L1
               DUMN = DUMN+SCR(I,N)*VA(N,J)
            END DO
            PMO(IJ,IVEC) = DUMN
         END DO
      END DO
      RETURN
      END SUBROUTINE MRSFQMNTOIA
C
C*MODULE MRSF   *DECK MRSFROESUM
C>
C>    @brief  clone of SFROESUM
C>            for singlet and triplet ROKS/MRSF
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Optimized cycles (4 lvl -> 2 lvl)
C>
C>    @date    Mar, 2022 Performance improvements
C>
      subroutine mrsfroesum(tmp1,tmp2,xia,fij,fab,pmo,l1,lx,l7,noca,    &
     &                      nocb,ivec,lmrsf,mrsfs,mrsft)
C
      USE prec, ONLY: dp
      USE comm_INFOXR, ONLY: IXCORE
      USE mx_limits, only: mxao
      use constants, ONLY: zero, one, two
      implicit none
c
      logical :: dskwrk, goparr, maswrk
      integer :: ibtyp, iptim, master, me, nproc
      common /par   / me, master, nproc, ibtyp, iptim, goparr, dskwrk,  &
     &                maswrk
c
      integer :: iv, ivec, l1, l7, lx, noca, nocb
      logical :: lmrsf, mrsfs, mrsft
      real(kind=dp), dimension(lx,*) :: fab, fij, tmp1, tmp2, xia
      real(kind=dp), dimension(lx,lx) :: scr, sc1
      real(kind=dp), dimension(l7,*) :: pmo
c
      real(kind=dp) :: dumn, xlr
      integer :: i, i1, i2, ij, j, j1, j2, lr1, lr2, ijlr1, ijlr2

      integer idamax

      if(ixcore(idamax(mxao,ixcore,1)).ne.0) then
        do i1=1,nocb
           i=1
           do j=1,mxao
              if(ixcore(j).eq.i1) i=0
           end do
           if(i.eq.1) then
             fij(i1,i1)=-1.0d6
           end if
        end do
      end if 

      iv = ivec
      lr1 = noca-1
      lr2 = noca

      call dcopy(lx*lx,xia,1,scr,1)
      scr(lr1,lr1) = zero
      scr(lr2,lr2) = zero

! First normal terms
c contraction1
      do j1=nocb+1,lx
      do i2=1,noca
         dumn=zero
         do j2=nocb+1,lx
            dumn=dumn+fab(j1,j2)*scr(i2,j2)
         end do
         tmp1(i2,j1) = dumn
      end do
      end do
c     call dgemm('n','t',noca,lx,lx,one,
c    *           scr(1,nocb+1),lx,
c    *           fab(nocb-1,nocb-1),lx,
c    *           zero,tmp2,lx)
c     write(*,*) "kkk => tmp1"
c     call prsql(tmp1,lx,lx,lx)
c     write(*,*) "kkk => tmp2"
c     call prsql(tmp2,lx,lx,lx)

c contraction2
      do j2=nocb+1,lx
      do i1=1,noca
         dumn=zero
         do i2=1,noca
            dumn=dumn+fij(i1,i2)*scr(i2,j2)
         end do
         tmp2(i1,j2) = dumn
      end do
      end do

      ijlr2 = (noca  -nocb-1)*noca+noca
      ijlr1 = (noca-1-nocb-1)*noca+noca-1
      xlr = xia(lr1,lr1)
      if (mrsfs) then
         ij = ijlr1
         dumn = zero
         do i2 = 1,noca
            dumn = dumn-fij(lr1,i2)*scr(i2,lr1)
            dumn = dumn+fij(lr2,i2)*scr(i2,lr2)
         enddo
         do j2 = nocb+1,lx
            dumn = dumn+fab(lr1,j2)*scr(lr1,j2)
            dumn = dumn-fab(lr2,j2)*scr(lr2,j2)
         enddo
         pmo(ij,iv) = pmo(ij,iv)+dumn/sqrt(two)
     *              + xlr*(fab(lr1,lr1)+fab(lr2,lr2)
     *                    -fij(lr1,lr1)-fij(lr2,lr2))/two

         ij = 0
         do j = nocb+1,lx
         do i = 1,noca
            if(i.eq.lr2 .and. j.eq.lr2) cycle
            if(i.eq.lr1 .and. j.eq.lr1) then
               ij = ij+1
               cycle
            endif
            ij = ij+1
            pmo(ij,iv) = pmo(ij,iv)+tmp1(i,j)-tmp2(i,j)
            if(i.eq.lr1)pmo(ij,iv) = pmo(ij,iv)+fab(j,lr1)*xlr/sqrt(two)
            if(i.eq.lr2)pmo(ij,iv) = pmo(ij,iv)-fab(j,lr2)*xlr/sqrt(two)
            if(j.eq.lr1)pmo(ij,iv) = pmo(ij,iv)-fij(i,lr1)*xlr/sqrt(two)
            if(j.eq.lr2)pmo(ij,iv) = pmo(ij,iv)+fij(i,lr2)*xlr/sqrt(two)
         end do
         end do
      elseif (mrsft) then
         ij = ijlr1
         dumn = zero
         do i2 = 1,noca
            dumn = dumn-fij(lr1,i2)*scr(i2,lr1)
            dumn = dumn-fij(lr2,i2)*scr(i2,lr2)
         enddo
         do j2 = nocb+1,lx
            dumn = dumn+fab(lr1,j2)*scr(lr1,j2)
            dumn = dumn+fab(lr2,j2)*scr(lr2,j2)
         enddo
         pmo(ij,iv) = pmo(ij,iv)+dumn/sqrt(two)
     *              + xlr*(fab(lr1,lr1)+fab(lr2,lr2)
     *                    -fij(lr1,lr1)-fij(lr2,lr2))/two
         ij = 0
         do j = nocb+1,lx
         do i = 1,noca
            if(i.eq.lr1 .and. j.eq.lr2) cycle
            if(i.eq.lr2 .and. j.eq.lr1) cycle
            if(i.eq.lr2 .and. j.eq.lr2) cycle
            if(i.eq.lr1 .and. j.eq.lr1) then
               ij = ij + 1
               cycle
            endif
            ij = ij+1
            pmo(ij,iv) = pmo(ij,iv)+tmp1(i,j)-tmp2(i,j)
            if(i.eq.lr1)pmo(ij,iv) = pmo(ij,iv)+fab(j,lr1)*xlr/sqrt(two)
            if(i.eq.lr2)pmo(ij,iv) = pmo(ij,iv)+fab(j,lr2)*xlr/sqrt(two)
            if(j.eq.lr1)pmo(ij,iv) = pmo(ij,iv)-fij(i,lr1)*xlr/sqrt(two)
            if(j.eq.lr2)pmo(ij,iv) = pmo(ij,iv)-fij(i,lr2)*xlr/sqrt(two)
         end do
         end do
      endif
      return
      end subroutine mrsfroesum
C
C*MODULE MRSF   *DECK MRSFQROESUM
C>
C>    @brief  clone of SFROESUM
C>            for quintet ROKS/MRSF
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE MRSFQROESUM(FBZ,ZFA,PMO,LX,L7,NOCA,NOCB,IVEC)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IVEC, L7, LX, NOCA, NOCB
      REAL(KIND=dp), DIMENSION(NOCB,LX) :: FBZ, ZFA
      REAL(KIND=dp), DIMENSION(L7,*) :: PMO
C
      INTEGER :: I, IJ, J
C
      IJ = 0
      DO J=NOCA+1,LX
         DO I=1,NOCB
            IJ = IJ + 1
c           ZFA - FBZ
            PMO(IJ,IVEC) = PMO(IJ,IVEC) + ZFA(I,J) - FBZ(I,J)
         END DO
      END DO
      RETURN
      END SUBROUTINE MRSFQROESUM
c
C*MODULE MRSF   *DECK MRSFDMAT
C>
C>    @brief  clone of SFDMAT
C>            for singlet and triplet MRSF
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE mrsfdmat(r,rao,xv12,ca,cb,da,db,ta,tb,wrk1,wrk2,noca,  &
     &                    nocb,nvirb,l1,l2,lx,l7,mrsfs,mrsft,mrsfq)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, one, two
      implicit none
c
      integer :: idaf, ip, ipk, ir, is, iw, nav
      integer, dimension(950) :: ioda
      common /iofile/ ir, iw, ip, is, ipk, idaf, nav, ioda
c
      integer :: l1, l2, l7, lx, noca, nocb, nvirb
      logical :: mrsfq, mrsfs, mrsft
      real(kind=dp), dimension(l1,lx) :: ca, cb
      real(kind=dp), dimension(l2) :: da, db, ta, tb
      real(kind=dp), dimension(l7) :: r
      real(kind=dp), dimension(l1,l1) :: rao
      real(kind=dp), dimension(*) :: wrk1, wrk2
      real(kind=dp), dimension(noca,nvirb) :: xv12
c
      integer :: i, ij, ijd, ijg, ijlr1, ijlr2, j

      ijlr1=(noca-1-nocb-1)*noca+noca-1
      ijg  =(noca-1-nocb-1)*noca+noca
      ijd  =(noca  -nocb-1)*noca+noca-1
      ijlr2=(noca  -nocb-1)*noca+noca

      if (mrsfs) then
         do i=1,noca
            do j=nocb+1,lx
               ij  =(j-nocb-1)*noca+i
               if(ij.eq.ijlr1) then
                  xv12(i,j-nocb)= r(ijlr1)/sqrt(two)
                  cycle
               else if(ij.eq.ijlr2) then
                  xv12(i,j-nocb)=-r(ijlr1)/sqrt(two)
                  cycle
               endif
               if(ij.gt.ijlr2) ij=ij-1
               xv12(i,j-nocb)=r(ij)
            enddo
         enddo
      else if (mrsft) then
         do i=1,noca
            do j=nocb+1,lx
               ij  =(j-nocb-1)*noca+i
               if(ij.eq.ijlr1) then
                  xv12(i,j-nocb)= r(ijlr1)/sqrt(two)
                  cycle
               else if(ij.eq.ijlr2) then
                  xv12(i,j-nocb)= r(ijlr1)/sqrt(two)
                  cycle
               else if(ij.eq.ijg  ) then
                  xv12(i,j-nocb)=zero
                  cycle
               else if(ij.eq.ijd  ) then
                  xv12(i,j-nocb)=zero
                  cycle
               endif
               if(ij.gt.ijg .and. ij.lt.ijd) then
                  ij=ij-1
               else if(ij.gt.ijd .and. ij.lt.ijlr2) then
                  ij=ij-2
               else if(ij.gt.ijlr2) then
                  ij=ij-3
               endif
               xv12(i,j-nocb)=r(ij)
            enddo
         enddo
      endif
c
c mo(i+,a-) -> ao(m,n)
c
      CALL SFIATOGEN(XV12,WRK1,LX,L7,NOCA,NOCB,1,1,
     *               .false.,.false.,.false.,.false.)
      CALL DGEMM('N','N',L1,LX,LX,ONE,CA,L1,WRK1,LX,ZERO,WRK2,L1)
      CALL DGEMM('N','T',L1,L1,LX,ONE,WRK2,L1,CB,L1,ZERO,RAO,L1)
C
C     ----- UNRELAXED DIFFERENCE DENSITY MATRIX -----
C
C OCC(ALPHA)-OCC(ALPHA)
C
      CALL DGEMM('N','T',NOCA,NOCA,NVIRB,-ONE,XV12,NOCA,XV12,
     *           NOCA,ZERO,WRK1,NOCA)
C MO(I+,J+) -> AO(M,N)
      CALL DGEMM('N','N',L1,NOCA,NOCA,ONE,CA,L1,WRK1,NOCA,ZERO,
     *           WRK2,L1)
      CALL DGEMM('N','T',L1,L1,NOCA,ONE,WRK2,L1,CA,L1,ZERO,WRK1,L1)
      CALL GENTOCAN(WRK1,TA,L1)
C TOTAL (= REFERENCE + UNRELAXED) DENSITY FOR ALPHA
      CALL DAREAD(IDAF,IODA,DA,L2,418,0)
      CALL VADD(DA,1,TA,1,DA,1,L2)
C
C VIRT(BETA)-VIRT(BETA)
C
      CALL DGEMM('T','N',NVIRB,NVIRB,NOCA,ONE,XV12,NOCA,XV12,
     *           NOCA,ZERO,WRK1,NVIRB)
C MO(A-,B-) -> AO(M,N)
      CALL DGEMM('N','N',L1,NVIRB,NVIRB,ONE,CB(1,NOCB+1),L1,
     *           WRK1,NVIRB,ZERO,WRK2,L1)
      CALL DGEMM('N','T',L1,L1,NVIRB,ONE,WRK2,L1,CB(1,NOCB+1),
     *           L1,ZERO,WRK1,L1)
      CALL GENTOCAN(WRK1,TB,L1)
C TOTAL (= REFERENCE + UNRELAXED) DENSITY FOR BETA
      CALL DAREAD(IDAF,IODA,DB,L2,428,0)
      CALL VADD(DB,1,TB,1,DB,1,L2)
c
      return
      end subroutine mrsfdmat

C*MODULE MRSF   *DECK MRSFQDMAT
C>
C>    @brief  clone of SFDMAT
C>            for quintet MRSF
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE MRSFQDMAT(R,RAO,CA,CB,DA,DB,TA,TB,WRK1,WRK2,NOCA,NOCB, &
     &                     NVIRA,L1,L2,LX,L7)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, one
      IMPLICIT NONE
C
      INTEGER :: IDAF, IP, IPK, IR, IS, IW, NAV
      INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
C
      INTEGER :: L1, L2, L7, LX, NOCA, NOCB, NVIRA
      REAL(KIND=dp), DIMENSION(L1,LX) :: CA, CB
      REAL(KIND=dp), DIMENSION(L2) :: DA, DB, TA, TB
      REAL(KIND=dp), DIMENSION(NOCB,NVIRA) :: R
      REAL(KIND=dp), DIMENSION(L1,L1) :: RAO
      REAL(KIND=dp), DIMENSION(*) :: WRK1, WRK2
C
C MO(I+,A-) -> AO(M,N)
C
      CALL SFIATOGEN(R,WRK1,LX,L7,NOCA,NOCB,1,1,
     *               .true.,.false.,.false.,.true.)
      CALL DGEMM('N','N',L1,LX,LX,ONE,CB,L1,WRK1,LX,ZERO,WRK2,L1)
      CALL DGEMM('N','T',L1,L1,LX,ONE,WRK2,L1,CA,L1,ZERO,RAO,L1)
C
C     ----- UNRELAXED DIFFERENCE DENSITY MATRIX -----
C
C VIRT(ALPHA)-VIRT(ALPHA)
C
      CALL DGEMM('T','N',NVIRA,NVIRA,NOCB,ONE,R,NOCB,R,NOCB,ZERO,
     *           WRK1,NVIRA)
C MO(A-,B-) -> AO(M,N)
      CALL DGEMM('N','N',L1,NVIRA,NVIRA,ONE,CA(1,NOCA+1),L1,
     *           WRK1,NVIRA,ZERO,WRK2,L1)
      CALL DGEMM('N','T',L1,L1,NVIRA,ONE,WRK2,L1,CA(1,NOCA+1),
     *           L1,ZERO,WRK1,L1)
      CALL GENTOCAN(WRK1,TA,L1)
C TOTAL (= REFERENCE + UNRELAXED) DENSITY FOR ALPHA
      CALL DAREAD(IDAF,IODA,DA,L2,418,0)
      CALL VADD(DA,1,TA,1,DA,1,L2)
C
C OCC(BETA)-OCC(BETA)
C
      CALL DGEMM('N','T',NOCB,NOCB,NVIRA,-ONE,R,NOCB,R,NOCB,ZERO,
     *           WRK1,NOCB)
C MO(I+,J+) -> AO(M,N)
      CALL DGEMM('N','N',L1,NOCB,NOCB,ONE,CB,L1,WRK1,NOCB,ZERO,WRK2,L1)
      CALL DGEMM('N','T',L1,L1,NOCB,ONE,WRK2,L1,CB,L1,ZERO,WRK1,L1)
      CALL GENTOCAN(WRK1,TB,L1)
C TOTAL (= REFERENCE + UNRELAXED) DENSITY FOR BETA
      CALL DAREAD(IDAF,IODA,DB,L2,428,0)
      CALL VADD(DB,1,TB,1,DB,1,L2)
C
      RETURN
      END SUBROUTINE MRSFQDMAT
C*MODULE MRSF   *DECK MRSFCBC
C>
C>    @brief    trial vectors in ao basis for MRSF
C>
C>    @details  eight types of trial vectors
C>              defined in Eq. (III.27) of JCP 150 184111 (2019)
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Combined 19 cycles into 5 and 2 dgemm
C>             The substract of (co2v and co1v) and (cco1 and cco2)
C>             was replaced by o21v and co12. This led to changes in
C>             mrsf2e, mrsftwoei, mrsfadisk, mrsfmntoia, mrsfsp and
C>             dabmrsfdft subroutines.
C>             As a result, reduced memory usage up to 15 percent.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE MRSFCBC(va,vb,bvec,
     &                   bo2v,bo1v,bco1,bco2,ball,o21v,co12,
     &                   tmp,l1,lx,noca,nocb,mrsfs,mrsft)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, one, two
      IMPLICIT NONE
C
        INTEGER :: IDAF, IP, IPK, IR, IS, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
C
      INTEGER :: L1, L3, LX, NOCA, NOCB
      LOGICAL :: MRSFS, MRSFT
      REAL(KIND=dp), DIMENSION(L1,*) :: BCO1, BCO2, BO1V, BO2V,
     &                            BALL, O21V, CO12
      REAL(KIND=dp), DIMENSION(l1,noca-2) :: TMP
      REAL(KIND=dp), DIMENSION(l1,4) :: TMP1, TMP2
      REAL(KIND=dp), DIMENSION(lx,*) :: BVEC
      REAL(KIND=dp), DIMENSION(l1,*) :: VA, VB
C
      INTEGER :: I, J, M

      l3 = l1*l1

      tmp1 = zero
      do j = noca+1,lx
         tmp1(:,1) = tmp1(:,1)+vb(:,j)*bvec(noca,j)
         tmp1(:,2) = tmp1(:,2)+vb(:,j)*bvec(noca-1,j)
         tmp1(:,3) = tmp1(:,3)+vb(:,j)*bvec(noca,j)
         tmp1(:,4) = tmp1(:,4)+vb(:,j)*bvec(noca-1,j)
      end do

      tmp2 = zero
      do i = 1,noca-2
         tmp2(:,1) = tmp2(:,1)+va(:,i)*bvec(i,noca-1)
         tmp2(:,2) = tmp2(:,2)+va(:,i)*bvec(i,noca)
         tmp2(:,3) = tmp2(:,3)+va(:,i)*bvec(i,noca-1)
         tmp2(:,4) = tmp2(:,4)+va(:,i)*bvec(i,noca)
      end do

      do m = 1,l1
         bo2v(:,m) = bo2v(:,m)+va(:,noca)*tmp1(m,1)
         bo1v(:,m) = bo1v(:,m)+va(:,noca-1)*tmp1(m,2)
         o21v(:,m) = o21v(:,m)+va(:,noca-1)*tmp1(m,3)
     &                        -va(:,noca)*tmp1(m,4)
         bco1(:,m) = bco1(:,m)+vb(m,noca-1)*tmp2(:,1)
         bco2(:,m) = bco2(:,m)+vb(m,noca)*tmp2(:,2)
         co12(:,m) = co12(:,m)+vb(m,noca)*tmp2(:,3)
     &                        -vb(m,noca-1)*tmp2(:,4)
      end do

      call daxpy(l3,one,bo2v,1,ball,1)
      call daxpy(l3,one,bo1v,1,ball,1)
      call daxpy(l3,one,bco1,1,ball,1)
      call daxpy(l3,one,bco2,1,ball,1)

      tmp = zero
c      do i=1,noca-2
c      do m=1,l1
c      do j=noca+1,lx
c         tmp(m,i)=tmp(m,i)+vb(m,j)*bvec(i,j)
c      end do
c      end do
c      end do
      call dgemm('n','t',l1,noca-2,lx-noca,one,
     *           vb(1,noca+1),l1,
     *           bvec(1,noca+1),lx,
     *           zero,tmp,l1)

c      do n=1,l1
c      do m=1,l1
c      do i=1,noca-2
c         ball(n,m)=ball(n,m)+tmp(m,i)*va(n,i)
c      end do
c      end do
c      end do
      call dgemm('n','t',l1,l1,noca-2,one,
     *           va,l1,
     *           tmp,l1,
     *           one,ball,l1)

      if (mrsfs) then
         do m=1,l1
            ball(:,m) = ball(:,m)
     *         +va(:,noca)*bvec(noca,noca-1)*vb(m,noca-1)
     *         +va(:,noca-1)*bvec(noca-1,noca)*vb(m,noca)
     *         +(va(:,noca-1)*vb(m,noca-1)-va(:,noca)*vb(m,noca))
     *            *bvec(noca-1,noca-1)/sqrt(two)
         end do
      elseif (mrsft) then
         do m = 1,l1
            ball(:,m) = ball(:,m)
     *         +(va(:,noca-1)*vb(m,noca-1)+va(:,noca)*vb(m,noca))
     *            *bvec(noca-1,noca-1)/sqrt(two)
         end do
      endif

      return
      end subroutine mrsfcbc
C
c*MODULE MRSF   *deck mrbsym
C>
C>    @brief    irreducible representation
C>              of reference and excited states
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE MRBSYM(irmoa,irmob,labsym,inxov,ndsr,vro,noca,nocb,    &
     &                  nvirb,lx,l7,mxvec)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), PARAMETER :: THRESHOLD = 5.0D-02
C
      INTEGER, DIMENSION(14) :: IADDR1, IADDR2, IADDR3, IPA, IRPNAM,    &
     &                          LAMBD0, LAMBDA
      INTEGER, DIMENSION(2,14,14,14) :: IJREP
      INTEGER, DIMENSION(14,14) :: NIJREP
      COMMON /SYMMUL/ NIJREP, IJREP
      COMMON /SYMREP/ IRPNAM, IPA, LAMBDA, LAMBD0, IADDR1, IADDR2,      &
     &                IADDR3
C
      INTEGER :: L7, LX, MXVEC, NDSR, NOCA, NOCB, NVIRB
      INTEGER, DIMENSION(2,L7) :: INXOV
      INTEGER, DIMENSION(lx) :: IRMOA, IRMOB
      INTEGER, DIMENSION(0:ndsr) :: LABSYM
      REAL(KIND=dp), DIMENSION(l7,mxvec) :: VRO
C
      LOGICAL :: ABEL, CHECK, DSKWRK, GOPARR, MASWRK
      LOGICAL :: ABELPT
      INTEGER :: I, IINAM, IINAM0, IJ, IRP, ISTAT, J, JRP, NAMSAV
      REAL(KIND=dp) :: XDUM
C
      INTEGER, PARAMETER :: iqmark = transfer('????',0)
c
      abel = abelpt()
c
c      -- symmetry of reference state
c
      irp = irmoa(1)
      do i=2,noca
         jrp = irmoa(i)
         iinam = ijrep(2,1,irp,jrp)
         irp = iinam
      end do
c
      irp = iinam
      do i=1,nocb
         jrp = irmob(i)
         iinam = ijrep(2,1,irp,jrp)
         irp = iinam
      end do
      iinam0    = iinam
      labsym(0) = irpnam(iinam0)
c
      states: do istat=1,ndsr
c
c      -- symmetry of excited state
c     mike does not believe that this code always prints the correct
c     state symmetry, and has an n2 example to back this up.  he
c     originally used the lines that you see commented out, preferring
c     to print blanks for all states, rather than mistakes for some.
c
           check = .false.
           namsav = -179
           do ij=1,l7
              i = inxov(1,ij)
              j = inxov(2,ij)
              xdum = vro(ij,istat)
              if (abs(xdum).gt.threshold) then
                 irp = irmoa(i)
                 jrp = irmob(j)
                 iinam = ijrep(2,1,irp,jrp)
                 if(check) then
                    if(iinam.ne.namsav) then
                       labsym(istat) = iqmark
                       cycle
                    end if
                 else
                    check = .true.
                    namsav = iinam
                 end if
              end if
           end do
           iinam = ijrep(2,1,iinam,iinam0)
           labsym(istat) = irpnam(iinam)
      end do states
c
      return
c
 1000 FORMAT(1X,'SYMMETRY OF STATE ',I2,' =',4X,A4)
      END SUBROUTINE MRBSYM
C
C*MODULE MRSF  *DECK UTDOSCALCMR
C>
C>    @brief    clone of UTDOSCALC for MRSF
C>
C>    @details  computing transition dipole moment and
C>              oscillator strength of MRSF response states
C>
C>    @author   Hiroya Nakata
C>
      SUBROUTINE UTDOSCALCMR(OS,TXYZ,TDENA,TDENB,SCR,VRO,VLO,EE,VA,VB,  &
     &                       AX,AY,AZ,L1,LX,L2,L7,L7A,L7C,NDSR,INXOV,   &
     &                       lmrsf,mrsfs,mrsft,mrsfq,NTHST)
C
      USE mx_limits, ONLY: mxatm
      USE prec, ONLY: dp
      USE constants, ONLY: zero, one
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(3,MXATM) :: C
      REAL(KIND=dp) :: DMX, DMY, DMZ, OMXXX, OMXXY, OMXXZ, OMXYY,       &
     &                 OMXYZ, OMXZZ, OMYYY, OMYYZ, OMYZZ, OMZZZ, OXXX,  &
     &                 OXXY, OXXZ, OXYY, OXYZ, OXZZ, OYYY, OYYZ, OYZZ,  &
     &                 OZZZ, QMXX, QMXY, QMXZ, QMYY, QMYZ, QMZZ, QXX,   &
     &                 QXY, QXZ, QYY, QYZ, QZZ, XP, YP, ZP
      INTEGER, DIMENSION(MXATM) :: IAN
      INTEGER :: ICH, IDAF, IP, IPK, IR, IS, IW, MUL, NA, NAT, NAV, NB, &
     &           NE, NQMT, NUM
      INTEGER, DIMENSION(950) :: IODA
      REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
      COMMON /XYZPRP/ XP, YP, ZP, DMX, DMY, DMZ, QXX, QYY, QZZ, QXY,    &
     &                QXZ, QYZ, QMXX, QMYY, QMZZ, QMXY, QMXZ, QMYZ,     &
     &                OXXX, OXXY, OXXZ, OXYY, OYYY, OYYZ, OXZZ, OYZZ,   &
     &                OZZZ, OXYZ, OMXXX, OMXXY, OMXXZ, OMXYY, OMYYY,    &
     &                OMYYZ, OMXZZ, OMYZZ, OMZZZ, OMXYZ
C
      INTEGER :: L1, L2, L7, L7A, L7C, LX, NDSR, NTHST
      LOGICAL :: LMRSF, MRSFQ, MRSFS, MRSFT
      REAL(KIND=dp), DIMENSION(L2) :: AX, AY, AZ
      REAL(KIND=dp), DIMENSION(*) :: EE
      INTEGER, DIMENSION(2,*) :: INXOV
      REAL(KIND=dp), DIMENSION(NDSR) :: OS
      REAL(KIND=dp), DIMENSION(L1,L1) :: SCR, TDENA, TDENB
      REAL(KIND=dp), DIMENSION(3,NDSR,NDSR) :: TXYZ
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA, VB
      REAL(KIND=dp), DIMENSION(L7,NDSR,*) :: VLO, VRO
C
      INTEGER :: IST, JST, NOCA, NOCB, NVIRB
C
C     --- CALCULATE OSCILLATOR STRENGTH
C
C        GET DIPOLE INTEGRALS AT CENTER OF MASS
C
      CALL CALCOM(XP,YP,ZP)
      CALL DIPINT(XP,YP,ZP,.FALSE.)
      CALL DAREAD(IDAF,IODA,AX,L2,95,0)
      CALL DAREAD(IDAF,IODA,AY,L2,96,0)
      CALL DAREAD(IDAF,IODA,AZ,L2,97,0)
c
      NOCA  = NA
      NOCB  = NB
C
      CALL VCLR(TXYZ,1,3*NDSR*NDSR)
C
      DO IST=1,NDSR
         DO JST=1,NDSR
C
            if(ist.eq.jst) cycle
C
            nvirb=lx-nocb
            CALL TCONSTMRSF(TDENA,VRO,LX,L7,NOCA,NOCB,nvirb,IST,JST,
     *           TDENB,SCR,mrsfs,mrsft,mrsfq)
C
            CALL DGEMM('N','N',L1,LX,LX,ONE,VA,L1,TDENA,LX,ZERO,SCR,L1)
            CALL DGEMM('N','T',L1,L1,LX,ONE,SCR,L1,VB,L1,ZERO,TDENA,L1)
C
            CALL TRADIP2(TXYZ,TDENA,AX,AY,AZ,L1,L2,NDSR,IST,JST)
C
         ENDDO
      ENDDO
C
      CALL MROSCALC2(OS,TXYZ,EE,NDSR)
C
      RETURN
      END SUBROUTINE UTDOSCALCMR
C*MODULE MRSF  *DECK MROSCALC2
C>
C>    @brief    clone of OSCALC for MRSF
C>
C>    @details  computing oscillator strength
C>              btw MRSF excited states
C>
C>    @author   Hiroya Nakata
C>
      SUBROUTINE MROSCALC2(OS,TXYZ,EE,NDSR)
      USE prec, ONLY: dp
      USE constants, ONLY: zero, two, three
      IMPLICIT NONE
C
      INTEGER :: NDSR
      REAL(KIND=dp), DIMENSION(NDSR) :: EE
      REAL(KIND=dp), DIMENSION(NDSR,NDSR) :: OS
      REAL(KIND=dp), DIMENSION(3,NDSR,NDSR) :: TXYZ
C
      REAL(KIND=dp) :: DUM
      INTEGER :: IST, JST, K
C
      DO IST=1,NDSR
      DO JST=1,NDSR
         OS(IST,JST)=ZERO
         DO K=1,3
            DUM=TXYZ(K,IST,JST)**2*ABS(EE(IST)-EE(JST))*TWO/THREE
            OS(IST,JST)=OS(IST,JST)+DUM
         ENDDO
      ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE MROSCALC2
c*MODULE MRSF   *deck sffase
C>
C>    @brief    keep phase and order of MO for NAMD
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE SFFASE(veca,vecb,ea,eb,l1,l2,l3,lx,noca,nocb)
C
      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      USE mx_limits, only: mxatm, mxgtot, mxsh
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER, PARAMETER :: MXNDRT = 20
C
      REAL(KIND=dp), DIMENSION(3,mxatm) :: C
      REAL(KIND=dp), DIMENSION(mxgtot) :: CD, CF, CG, CH, CI, CP, CS, EX
      LOGICAL :: DSKWRK, GOPARR, MASWRK
      INTEGER, DIMENSION(mxatm) :: IAN
      INTEGER :: IBTYP, ICH, IDAF, IJKO, IJKT, IP, IPTIM, IR, IW,       &
     &           MASTER, ME, MUL, NA, NAT, NAV, NB, NE, NPROC,          &
     &           NQMT, NSHELL, NUM
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(mxsh) :: KATOM, KLOC, KMAX, KMIN, KNG, KSTART, &
     &                            KTYPE
      REAL(KIND=dp), DIMENSION(1) :: X
      REAL(KIND=dp), DIMENSION(mxatm) :: ZAN
      COMMON /FMCOM / X
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
      COMMON /IOFILE/ IR, IW, IP, IJKO, IJKT, IDAF, NAV, IODA
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: L1, L2, L3, LX
      REAL(KIND=dp), DIMENSION(lx) :: EA, EB
      REAL(KIND=dp), DIMENSION(l1*lx) :: VECA, VECB
C
      REAL(KIND=dp), DIMENSION(lx,lx) :: DKSMSA, DKSMSB
      REAL(KIND=dp) :: DMAX, DUMMY, SIGNA, SIGNB
      REAL(KIND=dp), DIMENSION(lx) :: ETMPN
      INTEGER :: I, I2N, I2O, IJ, INDAS, ITMP, ITMP1, ITN, ITO, IVECOA, &
     &           IVECOB, J, LAST, LOADFM, NDXA, NDXB,
     &           NEED, NEXCH, NOCA, NOCB
      INTEGER, DIMENSION(lx) :: IDXSEL, IDXNEW
      LOGICAL :: MOEXCH, NDBUG, NEWRND
      REAL(KIND=dp), DIMENSION(l1,lx) :: VECTMP
c
      moexch=.false.
c
c     l1 : for atomic orbital index
C     lx : for molecular orbital index
c
      call valfm(loadfm)
      ivecoa = loadfm + 1
      ivecob = ivecoa + l1*lx
      indas  = ivecob + l1*lx
      itmp1  = indas  + l1*l1
      last   = itmp1  + lx*l1
      need = last - loadfm - 1
      if(maswrk) write(iw,9100) need
      call getfm(need)
c
c   get ao 1electron integral matrix:
c   ndas : <old(A)|new(B')>
      call coovlp(1,dummy,x(indas),l1,l1,l2,nat,
     *        mxgtot,nshell,
     *        ex,cs,cp,cd,cf,cg,ch,ci,
     *        kstart,katom,ktype,kng,kloc,kmin,kmax,
     *        mxgtot,nshell,
     *        ex,cs,cp,cd,cf,cg,ch,ci,
     *        kstart,katom,ktype,kng,kloc,kmin,kmax,
     *        c,cold)
c
      call daread(idaf,ioda,x(ivecoa),l1*lx,700,0)
      call daread(idaf,ioda,x(ivecob),l1*lx,701,0)
c
c        mo transformation : matrix ndms
c         <old(I)|new(J')> = trpose[cold(AI)] . <old(A)|new(B')> . cnew(B'J')
c            alpha
      call mrtrbr(x(ivecoa),l1,l1,lx,x(indas),l1,l1,x(itmp1),lx)
      call mrarbr(x(itmp1),lx,lx,l1,veca,l1,lx,dksmsa,lx)
c            beta
      call mrtrbr(x(ivecob),l1,l1,lx,x(indas),l1,l1,x(itmp1),lx)
      call mrarbr(x(itmp1),lx,lx,l1,vecb,l1,lx,dksmsb,lx)

      if(moexch) then

C     mo exchange algorithm needed to be changed
c     alpha ks mo
      do itmp=1,lx
          idxsel(itmp) = itmp
      enddo

      do i2o=1,lx
         dmax=0.0D+00
         do i2n=1,lx
             if(dmax.lt.abs(dksmsa(i2o,i2n)) .and.
     *          idxsel(i2n).ne.0) then
                 if ( (i2o.le.nocb .and. i2n.le.nocb) .or.
     *                (i2o.gt.noca .and. i2n.gt.noca) .or.
     *                (i2o.gt.nocb .and. i2o.le.noca .and.
     *                 i2n.gt.nocb .and. i2n.le.noca) ) then
                      dmax=abs(dksmsa(i2o,i2n))
                      idxnew(i2o)=i2n
                 endif
             endif
         enddo
         idxsel(idxnew(i2o))=0
      enddo

c assert all the indices are selected
      do i=1,lx
         if (idxsel(i).ne.0) CALL ABRT
      enddo

      CALL DCOPY(lx*l1,veca,1,vectmp,1)

c exchange values
      do i=1,lx
         itmp = idxnew(i)
         do j=1,l1
            ij=l1*(i-1)+j
            veca(ij)=vectmp(j,itmp)
         enddo
      enddo

c     beta ks mo
      do itmp=1,lx
          idxsel(itmp) = itmp
      enddo

      do i2o=1,lx
         dmax=0.0D0
         do i2n=1,lx
             if(dmax.lt.abs(dksmsb(i2o,i2n)) .and.
     *          idxsel(i2n).ne.0) then
                  dmax=abs(dksmsb(i2o,i2n))
                  idxnew(i2o)=i2n
             endif
         enddo
         idxsel(idxnew(i2o))=0
      enddo

clsh test
      if (maswrk) then
         do i2o=1,lx
            write(iw,*) i2o, idxnew(i2o)
         enddo
      endif

c assert all the indices are selected
      do i=1,lx
         if (idxsel(i).ne.0) CALL ABRT
      enddo

      CALL DCOPY(lx*l1,vecb,1,vectmp,1)

c exchange values
      do i=1,lx
         itmp = idxnew(i)
         do j=1,l1
            ij=l1*(i-1)+j
            vecb(ij)=vectmp(j,itmp)
         enddo
      enddo

      call mrtrbr(x(ivecoa),l1,l1,lx,x(indas),l1,l1,x(itmp1),lx)
      call mrarbr(x(itmp1),lx,lx,l1,veca,l1,lx,dksmsa,lx)
c            beta
      call mrtrbr(x(ivecob),l1,l1,lx,x(indas),l1,l1,x(itmp1),lx)
      call mrarbr(x(itmp1),lx,lx,l1,vecb,l1,lx,dksmsb,lx)

      endif

C        eigenvector sign correction of different time
      do i=1,lx
        signa=dksmsa(i,i)
        signb=dksmsb(i,i)
        if(signa.lt.0.0d0) then
          do j=1,l1
            ij=l1*(i-1)+j
            veca(ij) = -1.d+00 * veca(ij)
          enddo
        endif
        if(signb.lt.0.d+00) then
          do j=1,l1
            ij=l1*(i-1)+j
            vecb(ij) = -1.d+00 * vecb(ij)
          enddo
        endif
      enddo
c
      call retfm(need)

      if(maswrk.and.     moexch) write(iw,9101)
      if(maswrk.and..not.moexch) write(iw,9102)
      call timit(1)

 9100 FORMAT(1X,'MEMORY REQUIRED FOR MO PHASE AND ORDER IS=',I10,
     *          ' WORDS.')
 9101 FORMAT(1X,'..... DONE WITH ALIGNING MO PHASE AND ORDER .....')
 9102 FORMAT(1X,'..... DONE WITH ALIGNING MO PHASE .....'/
     *          ' WARNING: MO ORDER IS NOT ALIGNED')
      END SUBROUTINE SFFASE

C*MODULE MRSF   *deck sfxfase
C>
C>    @brief    keep phase of X amplitudes for NAMD
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE SFXFASE(V,VO,L7,NDSR)
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: L7, NDSR
      REAL(KIND=dp), DIMENSION(l7,ndsr) :: V, VO
C
      INTEGER :: I, J
      REAL(KIND=dp), DIMENSION(ndsr) :: SIGNX
C
C OVERAL SIGN CORRECTION BEFORE STATE TRACKING
       DO I=1,NDSR
          SIGNX(I)=0.D+00
          DO J=1,L7
             SIGNX(I) = SIGNX(I) + VO(J,I)*V(J,I)
          ENDDO
          IF( SIGNX(I).LT.0.0D+00 ) THEN
             DO J=1,L7
                V(J,I) = -1.0D+00 * V(J,I)
             ENDDO
          ENDIF
       ENDDO
      END SUBROUTINE SFXFASE
C*MODULE MRSF   *DECK TSR2E
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE TSR2E(BSMA,BSMB,ABA,ABB,DJ2EA,DJ2EB,DK2EA,DK2EB,BUF,   &
     &                 IBUF,NBF,NV,NMIC,NMAX)
C
      USE mx_limits, ONLY: mxsh, mxgtot, mxao
      USE camdft, ONLY: ALPHAC => cam_alpha, BETAC => cam_beta,
     *  CAMMU => cam_mu, CAMFLAG
      USE lrcdft, ONLY: LCFLAG, EMU, EMU2
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(20) :: DFTTYP
        REAL(KIND=dp) :: EXENA, EXENB, EXENC
        INTEGER :: IDFT34, NAUXFUN, NAUXSHL
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
        REAL(KIND=dp), DIMENSION(1) :: XX
      COMMON /FMCOM / XX
        INTEGER, DIMENSION(MXAO) :: IA
      COMMON /IJPAIR/ IA
        INTEGER :: IGRDTYP, INTTYP, NHEX, NINTMX, NTUPL
        LOGICAL :: PACK2E
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
        INTEGER :: IECP, IEFLD, ISCHWZ, NECP
      COMMON /INTOPT/ ISCHWZ, IECP, NECP, IEFLD
        INTEGER :: IDAF, IJK, IPK, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
        LOGICAL ::  LRINT
      COMMON /NLRCF / LRINT
        REAL(KIND=dp), DIMENSION(MXGTOT) :: CD, CF, CG, CH, CI, CP, CS, &
     &                                      EX
        INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG,       &
     &                              KSTART, KTYPE
        INTEGER :: NSHELL
      COMMON /NSHEL / EX, CS, CP, CD, CF, CG, CH, CI, KSTART, KATOM,    &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
        LOGICAL :: DIRSCF, FDIFF
      COMMON /OPTSCF/ DIRSCF, FDIFF
        INTEGER :: ICUT, ITOL, NOPK, NORMF, NORMP, NPRINT
      COMMON /OUTPUT/ NPRINT, ITOL, ICUT, NORMF, NORMP, NOPK
        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: IBTYP, IPTIM, MASTER, ME, NPROC
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK
C
      INTEGER :: NBF, NMAX, NMIC, NV
      REAL(KIND=dp), DIMENSION(NBF,NBF,NMIC,NV) :: ABA, ABB, DJ2EA,     &
     &                                             DJ2EB, DK2EA, DK2EB
      REAL(KIND=dp), DIMENSION(NBF*NBF*NMIC*NV) :: BSMA, BSMB
      REAL(KIND=dp), DIMENSION(NINTMX) :: BUF
      INTEGER, DIMENSION(NINTMX) :: IBUF
C
      INTEGER :: I, IV, J, L, NBF3, NINT, NSCHWZ
      LOGICAL :: SCHWRZ, TDSKWRK
      REAL(KIND=dp) :: TMPA, TMPB
C
      REAL(KIND=dp), PARAMETER :: RHF = transfer('RHF     ',1.0d0)
C
C     --- FORM SQUARE NON-SYMMETRIC FOCK-LIKE MATRIX ---
C         DIRECT METHOD = RECOMPUTE 2E- AO INTEGRALS
C         STANDARD METHOD = PROCESS INTEGRALS FROM DISK
C
      NINT = 0
      NSCHWZ = 0
      NBF3 = NBF*NBF

      IF(DIRSCF) THEN
clsh I have to include direct scf
c         SCHWRZ=ISCHWZ.GT.0
c         NSH2 = (NSHELL*NSHELL+NSHELL)/2
cC
c         CALL BASCHK(LMAX)
c                       NANGM =  4
c         IF(LMAX.EQ.2) NANGM =  6
c         IF(LMAX.EQ.3) NANGM = 10
c         IF(LMAX.EQ.4) NANGM = 15
c         IF(LMAX.EQ.5) NANGM = 21
c         IF(LMAX.EQ.6) NANGM = 28
c         MAXG = NANGM**4
cC
c         CALL VALFM(LOADFM)
c         LGHOND = LOADFM + 1
c         LXINTS = LGHOND + MAXG
c         LDSH   = LXINTS + NSH2
c         LAST   = LDSH     + NSH2
c         NEED   = LAST - LOADFM - 1
c         CALL GETFM(NEED)
cC
c         IF(SCHWRZ) THEN
c            DUMMY = 0.0D+00
c            CALL SHLTD(RHF,BALL,DUMMY,XX(LDSH),NBF,NSH2,NV)
c            CALL DAREAD(IDAF,IODA,XX(LXINTS),NSH2,54,0)
cC
c         END IF
c         CALL MRSFTWOEI(SCHWRZ,NINT,NSCHWZ,NBF,XX(LXINTS),NSH2,
c     *       XX(LGHOND),MAXG,IA,BALL,AGDLR,XX(LDSH),NV,.TRUE.,NMAX,
c     *                  BO2V,BO1V,BCO1,BCO2,BALL,
c     *                  CO2V,CO1V,CCO1,CCO2,AGDLR,ACO2V,ACO1V,
c     *                  ACCO1,ACCO2,ADO2V,ADO1V,ADCO1,ADCO2)
c         CALL RETFM(NEED)
      ELSE
         TDSKWRK = DSKWRK
         DSKWRK  = .TRUE.
         if(LCFLAG) then
          LRINT=.TRUE.
          CALL TSRADISK(BSMA,BSMB,NBF,DJ2EA,DJ2EB,DK2EA,DK2EB,
     *                  BUF,IBUF,NINTMX,NOPK,NMIC,NV)
          LRINT=.FALSE.
         endif
         IF(CAMFLAG) THEN
           LRINT     = .TRUE.
           EMU       = CAMMU
           EMU2      = CAMMU*CAMMU
           DFTTYP(3) = BETAC
          CALL TSRADISK(BSMA,BSMB,NBF,DJ2EA,DJ2EB,DK2EA,DK2EB,
     *                  BUF,IBUF,NINTMX,NOPK,NMIC,NV)
           DFTTYP(3) = ALPHAC
           LRINT     = .FALSE.
         endif
         CALL SEQREW(IJK)
         CALL TSRADISK(BSMA,BSMB,NBF,DJ2EA,DJ2EB,DK2EA,DK2EB,
     *                 BUF,IBUF,NINTMX,NOPK,NMIC,NV)
         DSKWRK  = TDSKWRK
      END IF
C
C     --- SUM UP PARTIAL FOCK-LIKE MATRICES ---
C
      IF(GOPARR) THEN
         CALL DDI_GSUMI(2311,NINT,1)
         CALL DDI_GSUMI(2312,NSCHWZ,1)
         CALL DDI_GSUMF(2313,DJ2EA,NBF3*NMIC*NMAX)
         CALL DDI_GSUMF(2314,DJ2EB,NBF3*NMIC*NMAX)
         CALL DDI_GSUMF(2315,DK2EA,NBF3*NMIC*NMAX)
         CALL DDI_GSUMF(2316,DK2EB,NBF3*NMIC*NMAX)
      END IF
C
      ABA = DJ2EA + DK2EA
      ABB = DJ2EB + DK2EB
C
      RETURN
      END SUBROUTINE TSR2E
C*MODULE MRSF   *DECK TSRADISK
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
      SUBROUTINE TSRADISK(BSJMA,BSJMB,L1,DJ2EA,DJ2EB,DK2EA,DK2EB,XX,IX, &
     &                    NINTMX,NOPK,NMIC,NV)
C
      USE mx_limits, ONLY: mxao, mxatm
      USE camdft, ONLY: CAMFLAG
      USE lrcdft, ONLY: LCFLAG, LRFILE
      USE constants, ONLY: half
      USE prec, ONLY: dp
      IMPLICIT NONE
C
        REAL(KIND=dp), DIMENSION(20) :: DFTTYP
        REAL(KIND=dp) :: EXENA, EXENB, EXENC
        INTEGER :: IDFT34, NAUXFUN, NAUXSHL
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
        INTEGER, DIMENSION(MXAO) :: IA
      COMMON /IJPAIR/ IA
        REAL(KIND=dp), DIMENSION(3,MXATM) :: C
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, C, IAN
        INTEGER :: IDAF, IS,IPK, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IS, IPK, IDAF, NAV, IODA
        LOGICAL ::  LRINT
      COMMON /NLRCF / LRINT
        INTEGER, DIMENSION(7) :: NORDER
      COMMON /ORDOPT/ NORDER
        INTEGER :: LABSIZ
      COMMON /PCKLAB/ LABSIZ
C
      INTEGER :: L1, NINTMX, NMIC, NOPK, NV
      REAL(KIND=dp), DIMENSION(L1,L1,NMIC,NV) :: BSJMA, BSJMB, DJ2EA,   &
     &                                           DJ2EB, DK2EA, DK2EB
      INTEGER, DIMENSION(*) :: IX
      REAL(KIND=dp), DIMENSION(NINTMX) :: XX
C
      REAL(KIND=dp) :: CSCALT, DIJ, DIKA, DIKB, DILA, DILB, DJI, DJKA,  &
     &                 DJKB, DJLA, DJLB, DKIA, DKIB, DKJA, DKJB, DKL,   &
     &                 DLIA, DLIB, DLJA, DLJB, DLK, HFSCAL, VAL, VALES, &
     &                 VALEX
      INTEGER :: I, IMIC, IMO, IPACK, J, JPACK, K, KPACK, L, LPACK, M,  &
     &           NIJ, NINT, LABEL, NKL, NPACK, NXX
C
C
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
      HFSCAL=DFTTYP(3)
      CSCALT=1.0D+00
      IF(LCFLAG) THEN
         IF(LRINT) THEN
            HFSCAL=1.0D+00
            CSCALT=0.0D+00
         ELSE
            HFSCAL=0.0D+00
            CSCALT=1.0D+00
         ENDIF
      ENDIF
      IF(CAMFLAG.AND.LRINT) CSCALT=0.0D+00
      I = 0
      J = 0
      K = 0
      L = 0
      NXX = 1 ! for entering to the cycle
      IF(NOPK.NE.1) THEN
         WRITE(IW,*) 'NOPK.NE.1'
         WRITE(IW,*) 'REKS CALCULATION DOES NOT SUPPORT ',
     *               'THIS INTEGRAL TYPE'
         CALL ABRT
      END IF
C
C     ----- INTEGRALS ARE NOT IN SUPERMATRIX FORM (NOPK=.TRUE.) -----
C
      do while (nxx .gt. 0)
      if(.not.LRINT) CALL PREAD(IS,XX,IX,NXX,NINTMX)
      if(     LRINT) CALL PREAD(LRFILE,XX,IX,NXX,NINTMX)
      IF(NXX .NE. 0) THEN
         NINT = IABS(NXX)
         IF(NINT .GT. NINTMX) CALL ABRT
         DO IMO =1,NV                                                    !trial vector
            DO IMIC=1,NMIC                                               !microstate
               DO M = 1,NINT                                             !long loop over integrals; better to put it inside shorter loops
C
                  NPACK = M
                  IF (LABSIZ .EQ. 2) THEN
                     LABEL = IX(NPACK)
                     IPACK = ISHFT( LABEL, -48 )
                     JPACK = IAND( ISHFT( LABEL, -32 ), 65535 )
                     KPACK = IAND( ISHFT( LABEL, -16 ), 65535 )
                     LPACK = IAND( LABEL, 65535 )
                  ELSE IF (LABSIZ .EQ. 1) THEN
                     IF ( MOD(NPACK,2) .EQ. 0 ) THEN
                        LABEL = IX( NPACK/2 )
                        IPACK = IAND( ISHFT( LABEL, -24 ), 255 )
                        JPACK = IAND( ISHFT( LABEL, -16 ), 255 )
                        KPACK = IAND( ISHFT( LABEL,  -8 ), 255 )
                        LPACK = IAND( LABEL, 255 )
                     ELSE
                        LABEL = IX( (NPACK/2)+1 )
                        IPACK = ISHFT( LABEL, -56 )
                        JPACK = IAND( ISHFT( LABEL, -48 ), 255 )
                        KPACK = IAND( ISHFT( LABEL, -40 ), 255 )
                        LPACK = IAND( ISHFT( LABEL, -32 ), 255 )
                     END IF
                  END IF
                  I = IPACK
                  J = JPACK
                  K = KPACK
                  L = LPACK
C
                  VAL = XX(M)
                  NIJ = IA(I)+J
                  NKL = IA(K)+L
C
C           USING SQUARE CANONICAL INTEGRAL FILE
C
                  IF(NORDER(7) .EQ. 1) THEN
                     IF(NKL .GT. NIJ) CYCLE
                     IF(I .EQ. J) VAL=VAL*HALF
                     IF(K .EQ. L) VAL=VAL*HALF
                     IF(NIJ .EQ. NKL) VAL=VAL*HALF
                  END IF

                  VALEX  =  VAL * HFSCAL
                  VALES  =  VAL * CSCALT

                  DKJA=(BSJMA(K,J,IMIC,IMO)+BSJMA(J,K,IMIC,IMO))*HALF    !BSKMA(K,J,IMIC,IMO)
                  DLJA=(BSJMA(L,J,IMIC,IMO)+BSJMA(J,L,IMIC,IMO))*HALF    !BSKMA(L,J,IMIC,IMO)
                  DKIA=(BSJMA(K,I,IMIC,IMO)+BSJMA(I,K,IMIC,IMO))*HALF    !BSKMA(K,I,IMIC,IMO)
                  DLIA=(BSJMA(L,I,IMIC,IMO)+BSJMA(I,L,IMIC,IMO))*HALF    !BSKMA(L,I,IMIC,IMO)
                  DJLA=(BSJMA(J,L,IMIC,IMO)+BSJMA(L,J,IMIC,IMO))*HALF    !BSKMA(J,L,IMIC,IMO)
                  DILA=(BSJMA(I,L,IMIC,IMO)+BSJMA(L,I,IMIC,IMO))*HALF    !BSKMA(I,L,IMIC,IMO)
                  DJKA=(BSJMA(J,K,IMIC,IMO)+BSJMA(K,J,IMIC,IMO))*HALF    !BSKMA(J,K,IMIC,IMO)
                  DIKA=(BSJMA(I,K,IMIC,IMO)+BSJMA(K,I,IMIC,IMO))*HALF    !BSKMA(I,K,IMIC,IMO)

                  DK2EA(I,L,IMIC,IMO)=DK2EA(I,L,IMIC,IMO)-VALEX*DKJA     !(IJ|LK)
                  DK2EA(I,K,IMIC,IMO)=DK2EA(I,K,IMIC,IMO)-VALEX*DLJA     !(IJ|LK)
                  DK2EA(J,L,IMIC,IMO)=DK2EA(J,L,IMIC,IMO)-VALEX*DKIA     !(IJ|LK)
                  DK2EA(J,K,IMIC,IMO)=DK2EA(J,K,IMIC,IMO)-VALEX*DLIA     !(IJ|LK)
                  DK2EA(K,I,IMIC,IMO)=DK2EA(K,I,IMIC,IMO)-VALEX*DJLA     !(IJ|LK)
                  DK2EA(K,J,IMIC,IMO)=DK2EA(K,J,IMIC,IMO)-VALEX*DILA     !(IJ|LK)
                  DK2EA(L,I,IMIC,IMO)=DK2EA(L,I,IMIC,IMO)-VALEX*DJKA     !(IJ|LK)
                  DK2EA(L,J,IMIC,IMO)=DK2EA(L,J,IMIC,IMO)-VALEX*DIKA     !(IJ|LK)

                  DKJB=(BSJMB(K,J,IMIC,IMO)+BSJMB(J,K,IMIC,IMO))*HALF    !BSKMB(K,J,IMIC,IMO)
                  DLJB=(BSJMB(L,J,IMIC,IMO)+BSJMB(J,L,IMIC,IMO))*HALF    !BSKMB(L,J,IMIC,IMO)
                  DKIB=(BSJMB(K,I,IMIC,IMO)+BSJMB(I,K,IMIC,IMO))*HALF    !BSKMB(K,I,IMIC,IMO)
                  DLIB=(BSJMB(L,I,IMIC,IMO)+BSJMB(I,L,IMIC,IMO))*HALF    !BSKMB(L,I,IMIC,IMO)
                  DJLB=(BSJMB(J,L,IMIC,IMO)+BSJMB(L,J,IMIC,IMO))*HALF    !BSKMB(J,L,IMIC,IMO)
                  DILB=(BSJMB(I,L,IMIC,IMO)+BSJMB(L,I,IMIC,IMO))*HALF    !BSKMB(I,L,IMIC,IMO)
                  DJKB=(BSJMB(J,K,IMIC,IMO)+BSJMB(K,J,IMIC,IMO))*HALF    !BSKMB(J,K,IMIC,IMO)
                  DIKB=(BSJMB(I,K,IMIC,IMO)+BSJMB(K,I,IMIC,IMO))*HALF    !BSKMB(I,K,IMIC,IMO)

                  DK2EB(I,L,IMIC,IMO)=DK2EB(I,L,IMIC,IMO)-VALEX*DKJB     !(IJ|LK)
                  DK2EB(I,K,IMIC,IMO)=DK2EB(I,K,IMIC,IMO)-VALEX*DLJB     !(IJ|LK)
                  DK2EB(J,L,IMIC,IMO)=DK2EB(J,L,IMIC,IMO)-VALEX*DKIB     !(IJ|LK)
                  DK2EB(J,K,IMIC,IMO)=DK2EB(J,K,IMIC,IMO)-VALEX*DLIB     !(IJ|LK)
                  DK2EB(K,I,IMIC,IMO)=DK2EB(K,I,IMIC,IMO)-VALEX*DJLB     !(IJ|LK)
                  DK2EB(K,J,IMIC,IMO)=DK2EB(K,J,IMIC,IMO)-VALEX*DILB     !(IJ|LK)
                  DK2EB(L,I,IMIC,IMO)=DK2EB(L,I,IMIC,IMO)-VALEX*DJKB     !(IJ|LK)
                  DK2EB(L,J,IMIC,IMO)=DK2EB(L,J,IMIC,IMO)-VALEX*DIKB     !(IJ|LK)

                  DKL=BSJMA(K,L,IMIC,IMO)+BSJMB(K,L,IMIC,IMO)
                  DLK=BSJMA(L,K,IMIC,IMO)+BSJMB(L,K,IMIC,IMO)
                  DIJ=BSJMA(I,J,IMIC,IMO)+BSJMB(I,J,IMIC,IMO)
                  DJI=BSJMA(J,I,IMIC,IMO)+BSJMB(J,I,IMIC,IMO)

                 DJ2EA(I,J,IMIC,IMO)=DJ2EA(I,J,IMIC,IMO)+VALES*(DKL+DLK) !(IJ|LK)
                 DJ2EA(J,I,IMIC,IMO)=DJ2EA(J,I,IMIC,IMO)+VALES*(DKL+DLK) !(IJ|LK)
                 DJ2EA(K,L,IMIC,IMO)=DJ2EA(K,L,IMIC,IMO)+VALES*(DIJ+DJI) !(IJ|LK)
                 DJ2EA(L,K,IMIC,IMO)=DJ2EA(L,K,IMIC,IMO)+VALES*(DIJ+DJI) !(IJ|LK)

                 DJ2EB(I,J,IMIC,IMO)=DJ2EB(I,J,IMIC,IMO)+VALES*(DKL+DLK) !(IJ|LK)
                 DJ2EB(J,I,IMIC,IMO)=DJ2EB(J,I,IMIC,IMO)+VALES*(DKL+DLK) !(IJ|LK)
                 DJ2EB(K,L,IMIC,IMO)=DJ2EB(K,L,IMIC,IMO)+VALES*(DIJ+DJI) !(IJ|LK)
                 DJ2EB(L,K,IMIC,IMO)=DJ2EB(L,K,IMIC,IMO)+VALES*(DIJ+DJI) !(IJ|LK)
C
               END DO
            END DO
         END DO
C
      END IF
      END DO
C
      if(.not.LRINT) CALL SEQREW(IS)
      if(     LRINT) CALL SEQREW(LRFILE)
C
      RETURN
      END SUBROUTINE TSRADISK
C*MODULE MRSF  *DECK TDSRIM
C> @brief    Computing I_L for Ω calculation
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
C
      SUBROUTINE TDSRIM(DIWM,DEM,NMIC)
C
      USE constants, ONLY: half, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: NMIC
      REAL(KIND=dp), DIMENSION(NMIC) :: DEM, DIWM
C
      REAL(KIND=dp) :: TMP
C
      DIWM(1) = ONE
      DIWM(2) =-ONE
      TMP     = HALF*(DEM(1)-DEM(2))/(DEM(4)-DEM(3))
      DIWM(3) = TMP
      DIWM(4) =-TMP
      RETURN
      END SUBROUTINE TDSRIM

C*MODULE MRSF  *DECK TSROMAT
C> @brief    Computes n_{p,L}^{\sigma} - n_{q,L}^{\sigma} matrix
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
C
      SUBROUTINE TSROMAT(OMATA,OMATB,NCORE,NOCR,NOCS,LX,NMIC,RMDEEX)
C
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE constants, ONLY: zero, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
C
      INTEGER :: LX, NCORE, NMIC, NOCR, NOCS
      LOGICAL :: RMDEEX
      REAL(KIND=dp), DIMENSION(LX,LX,NMIC) :: OMATA, OMATB
C
      INTEGER :: I, IR, IRA, IRB, IS, ISA, ISB, J, JR, JS, L
      REAL(KIND=dp) :: NRA, NRB, NSA, NSB
C
      CALL VCLR(OMATA,1,LX*LX*NMIC)
      CALL VCLR(OMATB,1,LX*LX*NMIC)
C
      DO L=1,NMIC
         IRA=NMICRO(1,1,L)
         IF(IRA.EQ.1) NRA=ZERO
         IF(IRA.EQ.2) NRA=ONE
         ISA=NMICRO(2,1,L)
         IF(ISA.EQ.1) NSA=ZERO
         IF(ISA.EQ.2) NSA=ONE
         IRB=NMICRO(1,2,L)
         IF(IRB.EQ.1) NRB=ZERO
         IF(IRB.EQ.2) NRB=ONE
         ISB=NMICRO(2,2,L)
         IF(ISB.EQ.1) NSB=ZERO
         IF(ISB.EQ.2) NSB=ONE
C      CV BLOCK
C      VC BLOCK
         DO J=NOCS+1,LX
         DO I=1,NCORE
            OMATA(I,J,L)= ONE
            OMATB(I,J,L)= ONE
         ENDDO
         ENDDO
C      CR BLOCK
C      CS BLOCK
C      RC BLOCK
C      SC BLOCK
         JR=NOCR
         JS=NOCS
         DO I=1,NCORE
            OMATA(I,JR,L)= ONE-NRA
            OMATB(I,JR,L)= ONE-NRB
            OMATA(I,JS,L)= ONE-NSA
            OMATB(I,JS,L)= ONE-NSB
         ENDDO
C      RV BLOCK
C      SV BLOCK
C      VR BLOCK
C      VS BLOCK
         IR=NOCR
         IS=NOCS
         DO J=NOCS+1,LX
            OMATA(IR,J,L)= NRA
            OMATB(IR,J,L)= NRB
            OMATA(IS,J,L)= NSA
            OMATB(IS,J,L)= NSB
         ENDDO
C      RS BLOCK
C      SR BLOCK
         IR=NOCR
         JS=NOCS
         OMATA(IR,JS,L)= NRA-NSA
         OMATB(IR,JS,L)= NRB-NSB
clsh test sr transition
         IF(.NOT.RMDEEX) OMATA(JS,IR,L)=-NRA+NSA
         IF(.NOT.RMDEEX) OMATB(JS,IR,L)=-NRB+NSB
      ENDDO

      RETURN
      END SUBROUTINE TSROMAT

C*MODULE MRSF   *DECK TSRFAO2MO
C> @brief    Transforms Fock matrices from AO to MO
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
C>
      SUBROUTINE TSRFAO2MO(FAMMO,FBMMO,FAMAO,FBMAO,VA,SCR1,SCR2,L1,L2,  &
     &                     L3,LX,NMIC)
C
      USE constants, ONLY: zero, one
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: L1, L2, L3, LX, NMIC
      REAL(KIND=dp), DIMENSION(L2,NMIC) :: FAMAO, FBMAO
      REAL(KIND=dp), DIMENSION(LX*LX,NMIC) :: FAMMO, FBMMO
      REAL(KIND=dp), DIMENSION(L3) :: SCR1, SCR2
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA
C
      INTEGER :: L
C
      DO L=1,NMIC      ! microstate
         CALL EXPND(FAMAO(1,L),SCR1,L1,0)
         CALL DGEMM('N','N',L1,LX,L1,ONE,SCR1,L1,VA,L1,
     *              ZERO,SCR2,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR2,L1,
     *              ZERO,FAMMO(1,L),LX)
         CALL EXPND(FBMAO(1,L),SCR1,L1,0)
         CALL DGEMM('N','N',L1,LX,L1,ONE,SCR1,L1,VA,L1,
     *              ZERO,SCR2,L1)
         CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR2,L1,
     *              ZERO,FBMMO(1,L),LX)
      ENDDO

      RETURN
      END SUBROUTINE TSRFAO2MO
C*MODULE MRSF   *DECK TSRWSUM
C> @brief     Takes summation over -L-
C>
C> @author     Seunghoon Lee
C>
C> @date       2021, Oct
C>
C>
      SUBROUTINE TSRWSUM(WSUM,AMO,BMO,LMO,OMATA,OMATB,WGT,NCORE,NOCS,LX,&
     &                   NMIC,L7)
C
      USE constants, ONLY: two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: L7, LMO, LX, NCORE, NMIC, NOCS
      REAL(KIND=dp), DIMENSION(LMO,LMO,NMIC) :: AMO, BMO
      REAL(KIND=dp), DIMENSION(LX,LX,NMIC) :: OMATA, OMATB
      REAL(KIND=dp), DIMENSION(NMIC) :: WGT
      REAL(KIND=dp), DIMENSION(L7) :: WSUM
C
      INTEGER :: I, IJ, J, L
      REAL(KIND=dp) :: TMPA, TMPB, WMTMP
C
      DO L=1,NMIC       ! microstate
         WMTMP = WGT(L)
         IF(L.GE.3) WMTMP = TWO*WMTMP

         IJ=0
         DO J=NCORE+1,LX
            DO I=1,NOCS
               IJ = IJ+1
               TMPA = OMATA(I,J,L)*AMO(I,J,L)
               TMPB = OMATB(I,J,L)*BMO(I,J,L)
C
               WSUM(IJ)=WSUM(IJ)+WMTMP*(TMPA+TMPB)
            END DO
         END DO
      END DO
C
      RETURN
      END SUBROUTINE TSRWSUM
C*MODULE MRSF   *DECK TSRABVEC
C>
C> @brief     Computes A*b product, where b is an arbitrary vector
C>
C> @author  Seunghoon Lee
C>
      SUBROUTINE TSRABVEC(ITYP,AMO,BVEC,CM,WRK1,WRK2,OMATA,OMATB,       &
     &              BSMA,BSMB,VA,OMG,BVAO,FAM,FBM,BUF,IBUF,amotest,     &
     &              SCR1,SCR2,SCR3,SCR4,J2EA,J2EB,K2EA,K2EB,            &
     &              GRD,WGT,DCH,RHOI,TAUI,AOMAX,GMO,FXC,TRAI,COEF,      &
     &              EX,EC,EX0,EC0,VALGA,                                &
     &              IAO,VPRGA,VPRGB,                                    &
     &              L1,L2,L3,LX,L7,NOCR,NOCS,NCORE,                     &
     &              ISTART,IEND,NMAX,NV,MXVEC,NMIC)
C>
C> ITYP                  0 - SA Hessian, 1 - PPS Hessian, 2 - OSS Hessian
C> AMO                   output: A*b product in MO
C> BVEC                  input: B vector in MO
C> CM                    C_L
C> WRK1,WRK2             scratch
C> OMATA,OMATB           were defined before
C> BSMA,BSMB             scratch
C> VA                    MOs
C> OMG                   input: Ω matrix
C> BVAO                  scratch
C> FAM,FBM               input: Fock for all mircostates in AO (lower triangle)
C> BUF, IBUF             scratch for DFT
C> amotest               scratch for debugging
C> SCR1,SCR2,SCR3,SCR4   scratch
C> J2EA,J2EB,K2EA,K2EB   scratch for 2-e integrals
C> G1                    G1
C> WPPS,WOSS             obvious
C> L1,L2,L3,LX,L7        the same
C> NOCR,NOCS,NCORE       the same
C> ISTART,IEND           for Davidson
C> NMAX                  max size of trial vector
C> NV                    size of the trial vector
C> MXVEC                 for Davidson
C> NMIC                  the number of current microstate
C>
      USE comm_REKSCM, ONLY: NMICRO, MTTYP, WPPS, WOSS, G1, DNR,
     * DNS, DELTA, FR, FS
      USE mx_limits, ONLY: mxgrid
      USE constants, ONLY: zero, one, two
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      REAL(KIND=dp), DIMENSION(137) :: BSLRD
      REAL(KIND=dp) :: DFTGTHR, DFTTHR, RHOMIN, SW0, SWOFF
      INTEGER :: IGRDTYP, ILENG, INTTYP, JANS, MAXGRD, NDFTFG,          &
     &           NHEX, NINTMX, NPHI, NPHI0, NRAD, NRAD0, NTHE, NTHE0,   &
     &           NTUPL
      INTEGER, DIMENSION(MXGRID) :: NANGPT, NANGPT0
      LOGICAL :: PACK2E, SG1
      COMMON /DFGRID/ DFTTHR, DFTGTHR, SWOFF, SW0, BSLRD, NDFTFG, NRAD, &
     &                NTHE, NPHI, NRAD0, NTHE0, NPHI0, NANGPT, NANGPT0, &
     &                SG1, JANS
      COMMON /INFGRD/ RHOMIN, ILENG, MAXGRD
      COMMON /INTFIL/ NINTMX, NHEX, NTUPL, PACK2E, INTTYP, IGRDTYP
C
      INTEGER :: IEND, ISTART, ITYP, L1, L2, L3, L7, LX, MXVEC, NCORE,  &
     &           NMAX, NMIC, NOCR, NOCS, NV
      REAL(KIND=dp), DIMENSION(L7,MXVEC) :: AMO, BVEC
      REAL(KIND=dp), DIMENSION(l7,mxvec) :: AMOTEST
      REAL(KIND=dp), DIMENSION(L1) :: AOMAX
      REAL(KIND=dp), DIMENSION(L3,NMIC,NMAX) :: BSMA, BSMB, WRK1, WRK2
      REAL(KIND=dp), DIMENSION(NINTMX) :: BUF
      REAL(KIND=dp), DIMENSION(L3,NMAX) :: BVAO
      REAL(KIND=dp), DIMENSION(NMIC) :: CM
      REAL(KIND=dp), DIMENSION(32*ILENG) :: COEF
      REAL(KIND=dp), DIMENSION(4*L1*ILENG) :: DCH
      REAL(KIND=dp), DIMENSION(35*ILENG) :: EC
      REAL(KIND=dp), DIMENSION(ILENG) :: EC0, EX0
      REAL(KIND=dp), DIMENSION(18*ILENG) :: EX
      REAL(KIND=dp), DIMENSION(L2,NMIC) :: FAM, FBM
      REAL(KIND=dp), DIMENSION(L2,2) :: FXC
      REAL(KIND=dp), DIMENSION(L1*8) :: GMO
      REAL(KIND=dp), DIMENSION(MAXGRD*3) :: GRD
      INTEGER, DIMENSION(*) :: IAO
      INTEGER, DIMENSION(NINTMX) :: IBUF
      REAL(KIND=dp), DIMENSION(L3,NMIC,NMAX) :: J2EA, J2EB, K2EA, K2EB
      REAL(KIND=dp), DIMENSION(LX,LX,NMIC) :: OMATA, OMATB
      REAL(KIND=dp), DIMENSION(L7) :: OMG
      REAL(KIND=dp), DIMENSION(8*MAXGRD) :: RHOI
      REAL(KIND=dp), DIMENSION(L3) :: SCR1, SCR2, SCR3, SCR4
      REAL(KIND=dp), DIMENSION(2*MAXGRD) :: TAUI
      REAL(KIND=dp), DIMENSION(10*ILENG) :: TRAI
      REAL(KIND=dp), DIMENSION(L1,LX) :: VA
      REAL(KIND=dp), DIMENSION(L1,2) :: VALGA
      REAL(KIND=dp), DIMENSION(*) :: VPRGA, VPRGB
      REAL(KIND=dp), DIMENSION(MAXGRD) :: WGT
C
      INTEGER :: I, IJ, IV, IVEC, IW, J, JI, L, LAMO, LAMOTMP, LBVAO,   &
     &           LBVAOTMP, LBVEC, LBVECTMP, LIV, LWRK1, LWRK1TMP, LWRK2,&
     &           LWRK2TMP, NPTGRD
      REAL(KIND=dp) :: OMGTMP, RHO, SCALTMP, SCALTMP2, TMP
C
      LOGICAL, PARAMETER :: DBGAMAT = .false.                            !set to .true. if you need to see the Hessian matrix
C
C Option for types of A matrix
C ITYP   = 0.....A_{SA} *BVEC
C          1.....A_{PPS}*BVEC
C          2.....A_{OSS}*BVEC
C
C        ----- PREPARING DENSITY LIKE MATRIX bs FOR 2E PART IN AO -----
      IV = 0
      DO IVEC=ISTART,IEND  ! trial b vector
         IV = IV + 1
         CALL SFIATOGEN(BVEC,SCR1,LX,L7,NOCS,NCORE,IVEC,MXVEC,
     *                  .false.,.false.,.false.,.false.)                 !expand IVEC row of BVEC matrix in a square LX.LX matrix SCR1

         DO L=1,NMIC       ! microstate

            DO J=1,LX
            DO I=1,LX
               IJ=(J-1)*LX+I
               SCR2(IJ)=OMATA(I,J,L)*SCR1(IJ)
               SCR3(IJ)=OMATB(I,J,L)*SCR1(IJ)
            ENDDO
            ENDDO

            CALL DGEMM('N','N',L1,LX,LX,ONE,VA,L1,SCR2,LX,
     *                 ZERO,SCR4,L1)
            CALL DGEMM('N','T',L1,L1,LX,ONE,SCR4,L1,VA,L1,
     *                 ZERO,BSMA(1,L,IV),L1)
            CALL DGEMM('N','N',L1,LX,LX,ONE,VA,L1,SCR3,LX,
     *                 ZERO,SCR4,L1)
            CALL DGEMM('N','T',L1,L1,LX,ONE,SCR4,L1,VA,L1,
     *                 ZERO,BSMB(1,L,IV),L1)
         END DO
      END DO
C
C        ----- CALCULATING 2E PART: J, K -----
C     Ab^{I,alpha}_{L} = Jbs^{I}_{L} - c_HF * Kbs^{I,alpha}_{L} AT LWRK1
C     Ab^{I,beta }_{L} = Jbs^{I}_{L} - c_HF * Kbs^{I,beta }_{L} AT LWRK2
C
      CALL VCLR(WRK1,1,L3*NMIC*NMAX)         ! Ab^{I,alpha}_{L}
      CALL VCLR(WRK2,1,L3*NMIC*NMAX)         ! Ab^{I,beta }_{L}
      CALL VCLR(J2EA,1,L3*NMIC*NMAX)
      CALL VCLR(J2EB,1,L3*NMIC*NMAX)
      CALL VCLR(K2EA,1,L3*NMIC*NMAX)
      CALL VCLR(K2EB,1,L3*NMIC*NMAX)

      CALL TSR2E(BSMA,BSMB,WRK1,WRK2,J2EA,J2EB,K2EA,K2EB,BUF,IBUF,
     *           L1,NV,NMIC,NMAX)

C
C        ----- CALCULATING 2E PART: X, C -----
C     Ab^{I,alpha}_{L} = Ab^{I,alpha}_{L} + c_DFT * Mbs^{I,alpha}_{L} + Nbs^{I,alpha}_{L} AT LWRK1
C     Ab^{I,beta }_{L} = Ab^{I,beta }_{L} + c_DFT * Mbs^{I,beta }_{L} + Nbs^{I,beta }_{L} AT LWRK2
C
      IF(NDFTFG.EQ.1) THEN
         IV = 0
         DO IVEC=ISTART,IEND  ! trial b vector
            IV = IV + 1
C SYMMETRIZED DENSITY-LIKE MATRIX
            DO L=1,NMIC       ! microstate
               DO I=1,L1
               DO J=1,L1
                  IJ=(I-1)*L1+J
                  JI=(J-1)*L1+I
                  SCR1(IJ)=0.5d+0*(BSMA(IJ,L,IV)+BSMA(JI,L,IV))
                  SCR2(IJ)=0.5d+0*(BSMB(IJ,L,IV)+BSMB(JI,L,IV))
               ENDDO
               ENDDO
               CALL DCOPY(L3,SCR1,1,BSMA(1,L,IV),1)
               CALL DCOPY(L3,SCR2,1,BSMB(1,L,IV),1)
            END DO
         END DO
C
         DO L=1,4
            MTTYP   = L
C        ----- GRID SETTING FOR EACH MICROSTATE -----
            NPTGRD = MAXGRD
            CALL UTDDFTSET(GRD,WGT,DCH,VA,VA,RHOI,TAUI,AOMAX,GMO,        !prepares the density of the L-th microstate; L is communicated to UDENCNST via the MTTYP constant
     *                     ILENG,NPTGRD,L1)

            IV = 0
            DO IVEC=ISTART,IEND  ! TRIAL B VECTOR
               IV = IV + 1

               CALL VCLR(FXC,1,L2*2)
               CALL UTDFXCP2(FXC,RHO,GRD,WGT,DCH,BSMA(1,L,IV),
     *                    BSMB(1,L,IV),RHOI,TAUI,TRAI,COEF,EX,EC,EX0,
     *                    EC0,AOMAX,VPRGA,VPRGB,IAO,ILENG,NPTGRD,
     *                    L1,L2,2,.FALSE.,0)
C
               CALL EXPND(FXC(1,1),SCR1,L1,0)
               CALL EXPND(FXC(1,2),SCR2,L1,0)

               CALL VADD(WRK1(1,L,IV),1,SCR1,1,WRK1(1,L,IV),1,L3)
               CALL VADD(WRK2(1,L,IV),1,SCR2,1,WRK2(1,L,IV),1,L3)

            ENDDO
         ENDDO
      ENDIF
C
C        ----- END OF 2E PART ----
C
      IV = 0
      DO IVEC=ISTART,IEND
         IV = IV + 1
         CALL SFIATOGEN(BVEC,SCR1,LX,L7,NOCS,NCORE,IVEC,MXVEC,
     *                  .false.,.false.,.false.,.false.)                 !expand IVEC row of BVEC matrix in a square LX.LX matrix SCR1
         CALL ASYMTRZE(SCR1,LX,LX,TMP)
         CALL DGEMM('N','T',LX,L1,LX,ONE,SCR1,LX,VA,L1,
     *              ZERO,BVAO(1,IV),LX)
      ENDDO
C
C        ----- CALCULATING 1E PART: Fb - bF -----
C     Ab^{I,alpha}_{L} = Ab^{I,alpha}_{L} + Fb^{I,alpha}_{L} - bF^{I,alpha}_{L} AT LWRK1
C     Ab^{I,beta }_{L} = Ab^{I,beta }_{L} + Fb^{I,beta }_{L} - bF^{I,beta }_{L} AT LWRK2
C
      DO L=1,NMIC      ! microstate
         CALL EXPND(FAM(1,L),SCR3,L1,0)
         CALL DGEMM('N','N',L1,LX,L1,ONE,SCR3,L1,VA,L1,ZERO,SCR1,L1)
         CALL EXPND(FBM(1,L),SCR3,L1,0)
         CALL DGEMM('N','N',L1,LX,L1,ONE,SCR3,L1,VA,L1,ZERO,SCR2,L1)
C
         IV = 0
         DO IVEC=ISTART,IEND  ! trial b vector
            IV = IV + 1
            LBVAOTMP = LBVAO + (IV-1)*L3
            LIV = (IV-1)*NMIC + L
            LWRK1TMP = LWRK1 + (LIV-1)*L3
            LWRK2TMP = LWRK2 + (LIV-1)*L3
C
            CALL DGEMM('N','N',L1,LX,L1,ONE,WRK1(1,L,IV),L1,VA,L1,
     *                 ZERO,SCR3,L1)
            CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR3,L1,
     *                 ZERO,WRK1(1,L,IV),L1)
            CALL DGEMM('N','N',L1,LX,L1,ONE,WRK2(1,L,IV),L1,VA,L1,
     *                 ZERO,SCR3,L1)
            CALL DGEMM('T','N',LX,LX,L1,ONE,VA,L1,SCR3,L1,
     *                 ZERO,WRK2(1,L,IV),L1)
C           Alpha
            CALL DGEMM('N','N',LX,LX,L1,ONE,BVAO(1,IV),LX,SCR1,L1,
     *                 ZERO,SCR3,LX)
            CALL SYMTRZE(SCR3,LX,LX)
            CALL DSCAL(LX*LX,TWO,SCR3,1)
            CALL VADD(WRK1(1,L,IV),1,SCR3,1,WRK1(1,L,IV),1,L3)
C           Beta
            CALL DGEMM('N','N',LX,L1,L1,ONE,BVAO(1,IV),L1,SCR2,L1,
     *                 ZERO,SCR3,L1)
            CALL SYMTRZE(SCR3,LX,LX)
            CALL DSCAL(LX*LX,TWO,SCR3,1)
            CALL VADD(WRK2(1,L,IV),1,SCR3,1,WRK2(1,L,IV),1,L3)
         ENDDO
      ENDDO
C        ----- PREPARING DENSITY LIKE MATRIX bn FOR OCCUPATION PART IN AO -----
      IF(ITYP.EQ.0) SCALTMP=-G1*WPPS
      IF(ITYP.EQ.1) SCALTMP=-G1
      IF(ITYP.EQ.2) SCALTMP=0
      IV = 0
      DO IVEC=ISTART,IEND  ! trial b vector
         IV = IV + 1
         LAMOTMP = LAMO + (IVEC-1)*L7
         LBVECTMP= LBVEC+ (IVEC-1)*L7
C
         OMGTMP=0.0D+00
         IJ = 0
         DO J=NCORE+1,LX
            DO I=1,NOCS
               IJ = IJ+1
               TMP= OMG(IJ)*BVEC(IJ,IVEC)
               OMGTMP = OMGTMP + TMP
            END DO
         END DO
C
         SCALTMP2=SCALTMP*OMGTMP
         IJ = 0
         DO J=NCORE+1,LX
            DO I=1,NOCS
               IJ = IJ + 1
               AMO(IJ,IVEC)=AMO(IJ,IVEC)+SCALTMP2*OMG(IJ)
            END DO
         END DO
C
         CALL TSRWSUM(AMO(1,IVEC),WRK1(1,1,IV),WRK2(1,1,IV),L1,OMATA,
     *                OMATB,CM,NCORE,NOCS,LX,NMIC,L7)
      END DO
C        ----- END OF OCCUPATION PART ----
      RETURN
      END SUBROUTINE TSRABVEC
c
C*MODULE MRSF  *DECK JKDMRSPNFLP
C>
C>    @brief    clone of JKDSPNFLP
C>              for singlet and triplet MRSF
C>
C>    @author   Seunghoon Lee
C>
      subroutine JKDMRSPNFLP(DA,DB,PA,PB,V,VSP,L3,L2)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IDAF, IJK, IP, IPK, IR, IW, NAV
      INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
C
      INTEGER :: L2, L3
      REAL(KIND=dp), DIMENSION(*) :: DA, DB, PA, PB, V, VSP
C
      REAL(KIND=dp) :: DUMA, DUMB
      INTEGER :: I
C
      CALL DAREAD(IDAF,IODA, DA,  L2,418,0)
      CALL DAREAD(IDAF,IODA, DB,  L2,428,0)
      CALL DAREAD(IDAF,IODA, PA,  L2,417,0)
      CALL DAREAD(IDAF,IODA, PB,  L2,427,0)
      CALL DAREAD(IDAF,IODA,  V,  L3,474,0)
      CALL DAREAD(IDAF,IODA,VSP,6*L3,479,0)
C
C FOR THE DENSITY MATRIX OF REFERENCE STATE,
C     DA = TOTAL (ALPHA+BETA) DENSITY, DB = SPIN (ALPHA-BETA) DENSITY
C FOR THE DIFFERENCE DENSITY MATRIX,
C     PA = TOTAL (ALPHA+BETA) DENSITY, PB = SPIN (ALPHA-BETA) DENSITY
C
      DO I=1,L2
         DUMA = DA(I)
         DUMB = DB(I)
         DA(I) = DUMA + DUMB
         DB(I) = DUMA - DUMB
C
         DUMA = PA(I)
         DUMB = PB(I)
         PA(I) = DUMA + DUMB
         PB(I) = DUMA - DUMB
      END DO
C
      RETURN
      END SUBROUTINE JKDMRSPNFLP
C
C*MODULE MRSF  *DECK JKDMRSPNFLPQ
C>
C>    @brief    clone of JKDSPNFLP
C>              for quintet MRSF
C>
C>    @author  Seunghoon Lee
C>
      SUBROUTINE JKDMRSPNFLPQ(DA,DB,PA,PB,V,L3,L2)
C
      USE prec, ONLY: dp
      IMPLICIT NONE
C
      INTEGER :: IDAF, IJK, IP, IPK, IR, IW, NAV
      INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
C
      INTEGER :: L2, L3
      REAL(KIND=dp), DIMENSION(*) :: DA, DB, PA, PB, V
C
      REAL(KIND=dp) :: DUMA, DUMB
      INTEGER :: I
C
      CALL DAREAD(IDAF,IODA,DA,L2,418,0)
      CALL DAREAD(IDAF,IODA,DB,L2,428,0)
      CALL DAREAD(IDAF,IODA,PA,L2,417,0)
      CALL DAREAD(IDAF,IODA,PB,L2,427,0)
      CALL DAREAD(IDAF,IODA, V,L3,474,0)
C
C FOR THE DENSITY MATRIX OF REFERENCE STATE,
C     DA = TOTAL (ALPHA+BETA) DENSITY, DB = SPIN (ALPHA-BETA) DENSITY
C FOR THE DIFFERENCE DENSITY MATRIX,
C     PA = TOTAL (ALPHA+BETA) DENSITY, PB = SPIN (ALPHA-BETA) DENSITY
C
      DO I=1,L2
         DUMA = DA(I)
         DUMB = DB(I)
         DA(I) = DUMA + DUMB
         DB(I) = DUMA - DUMB
C
         DUMA = PA(I)
         DUMB = PB(I)
         PA(I) = DUMA + DUMB
         PB(I) = DUMA - DUMB
      END DO
C
      RETURN
      END SUBROUTINE JKDMRSPNFLPQ
C
C*MODULE MRSF   *DECK DABMRSFDFT
C>
C>    @brief   clone of DABSFDFT
C>             for singlet and triplet MRSF
C>
C>    @author  Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Terms co2v, co1v, cco1 and cco2 was replaced
C>             by o21v, co12. mrsfcbc is initial subroutine of this.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE DABMRSFDFT(II,JJ,KK,LL,DA,DB,PA,PB,V,
     &                      BO2V,BO1V,BCO1,BCO2,O21V,CO12,
     &                      DAB,DABMAX,L1,L2,
     &                      Q4,POPLE,SPCP,MRSFS,MRSFT)
C
      USE mx_limits, ONLY: mxsh, mxgtot, mxao
      USE lrcdft, ONLY: LCFLAG
      USE camdft, ONLY: CAMFLAG
      USE constants, ONLY: zero, quarter, half, one, two, four
      USE prec, ONLY: dp
      IMPLICIT NONE
C
C DA: TOTAL (ALPHA+BETA) DENSITY FOR REFERENCE STATE
C DB: SPIN  (ALPHA-BETA) DENSITY FOR REFERENCE STATE
C PA: TOTAL (ALPHA+BETA) SF-CIS DIFFERENCE DENSITY
C PB: SPIN  (ALOHA-BETA) SF-CIS DIFFERENCE DENSITY
C  V: SF-CIS TRANSITION DENSITY
C
      REAL(KIND=dp), DIMENSION(MXGTOT) :: CCG, CD, CF, CH, CI, CP, CS,  &
     &                                    EX
      LOGICAL :: DBUG, IIEQJJ, IJEQKL, IJGTKL, IJLTKL, KKEQLL, LRINT,   &
     &           OUT, SOME
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      REAL(KIND=dp) :: EXENA, EXENB, EXENC
      INTEGER, DIMENSION(MXAO) :: IA
      INTEGER :: IDAF, IDFT34, IJK, IP, IPK, IR, IW, LA, LB, LC, LD,    &
     &           NAUXFUN, NAUXSHL, NAV, NSHELL
      INTEGER, DIMENSION(4,35) :: IGXYZ, JGXYZ, KGXYZ, LGXYZ
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG, KSTART, &
     &                            KTYPE
      REAL(KIND=dp), DIMENSION(84) :: PNRM
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /IJPAIR/ IA
      COMMON /INDD80/ LA, LB, LC, LD
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
      COMMON /MP2PRT/ SOME, OUT, DBUG
      COMMON /NLRCF / LRINT
      COMMON /NSHEL / EX, CS, CP, CD, CF, CCG, CH, CI, KSTART, KATOM,   &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
      COMMON /SHLEQU/ IIEQJJ, KKEQLL, IJEQKL, IJGTKL, IJLTKL
      COMMON /SHLLMN/ IGXYZ, JGXYZ, KGXYZ, LGXYZ
      COMMON /SHLNRM/ PNRM
C
      REAL(KIND=dp) :: DABMAX, Q4
      INTEGER :: II, JJ, KK, L1, L2, LL
      LOGICAL :: MRSFS, MRSFT, POPLE
      REAL(KIND=dp), DIMENSION(L1,*) :: BCO1, BCO2, BO1V, BO2V,
     &                                  O21V, CO12, V
      REAL(KIND=dp), DIMENSION(L2) :: DA, DB, PA, PB
      REAL(KIND=dp), DIMENSION(*) :: DAB
      REAL(KIND=dp), DIMENSION(3) :: SPCP
C
      REAL(KIND=dp) :: DB1, DB2, DC1, DC2, DC3, DC4, DD1, DD2, DD3,     &
     &                 DD4, DF1, DQ1, DT2, P1I, P2J, P3K, P4L, QFS,     &
     &                 QFSPCP1, QFSPCP2, QFSPCP3, SGNK
      INTEGER :: I, IJKL, J, JMAX, K, KKMAX, L, LANU, LASI, LMAX,       &
     &           LOCI, LOCJ, LOCK, LOCL, MAXI, MAXJ, MAXK, MAXL,        &
     &           MINI, MINJ, MINK, MINL, MULA, MUNU, MUSI, NLA,         &
     &           NMU, NNU, NSI, NUSI
C
CJMS  LABELLED COMMON GSPG80 DEFINED FOR COMPUTATIONAL EFFICIENCY.
CJMS  FOR SP BASES ONLY, IT CONTAINS THE E ARRAY WHICH IS THE DAB
CJMS  ARRAY WITH INDICES IN REVERSE ORDER: E(I,J,K,L)= DAB(L,K,J,I)
CJMS  AND IS USED IN SUB JKDG80 (MODULE GRD2B). IT ORIGINATES IN:
CJMS
CJMS     1. SUBS DABCLU, DABDFT, DABGVB, DABMC, DABMP2 AND DABUMP
CJMS        (MODULE GRD2A) AND SUB DABPAU (MODULE EFPAUL) WHICH ARE
CJMS        ALL CALLED BY SUB JKDER (MODULE GRD2A)
CJMS
CJMS     2. SUB DABCLU (MODULE GRD2A) WHICH IS CALLED BY SUB EFDEN OF
CJMS        MODULE EFGRD2
CJMS
CJMS     3. SUB PAR2PDM (MODULE MP2DDI) WHICH IS CALLED BY SUB PJKDMP2
CJMS        OF MODULE MP2DDI
C
      REAL(KIND=dp), DIMENSION(4,4,4,4) :: E
      COMMON /GSPG80/ E
C
C
      QFS = DFTTYP(3)*QUARTER
      QFSPCP1 = QFS
      QFSPCP2 = QFS
      QFSPCP3 = QFS
      IF (SPCP(1)/=-ONE) QFSPCP1 = SPCP(1)*QUARTER
      IF (SPCP(2)/=-ONE) QFSPCP2 = SPCP(2)*QUARTER
      IF (SPCP(3)/=-ONE) QFSPCP3 = SPCP(3)*QUARTER
C
      SGNK = 1.0
      IF (MRSFT) SGNK = -1.0
C
C     ----- FORM TWO-PARTICLE DENSITY MATRIX FOR CIS GRADIENT -----
C
      DABMAX = ZERO
      MINI = KMIN(II)
      MINJ = KMIN(JJ)
      MINK = KMIN(KK)
      MINL = KMIN(LL)
      LOCI = KLOC(II) - MINI
      LOCJ = KLOC(JJ) - MINJ
      LOCK = KLOC(KK) - MINK
      LOCL = KLOC(LL) - MINL
C
      IF (POPLE) THEN
         DO L = 1, LD
            NNU = LOCL + L
            DO K = 1, LC
               NMU = LOCK + K
               MUNU = IA(MAX0(NMU,NNU)) + MIN0(NMU,NNU)
               DO J = 1, LB
                  NSI = LOCJ + J
                  DO I = 1, LA
                     NLA = LOCI + I
                     LASI = IA(MAX0(NLA,NSI)) + MIN0(NLA,NSI)
C
C                    WRITE(IW,*) 'NLA,NSI,NMU,NNU,LASI,MUNU=',
C    *                            NLA,NSI,NMU,NNU,LASI,MUNU
C
                     MUSI = IA(MAX0(NMU,NSI)) + MIN0(NMU,NSI)
                     LANU = IA(MAX0(NLA,NNU)) + MIN0(NLA,NNU)
                     MULA = IA(MAX0(NMU,NLA)) + MIN0(NMU,NLA)
                     NUSI = IA(MAX0(NNU,NSI)) + MIN0(NNU,NSI)
C
                     DF1 =  (DA(MUNU) + PA(MUNU))*DA(LASI)
     *                    +  DA(MUNU)            *PA(LASI)
C
                     IF(QFS .NE.ZERO) THEN
                        DQ1 =  (DA(MUSI)+PA(MUSI))*DA(LANU)
     *                       +  DA(MUSI)          *PA(LANU)
     *                       + (DA(MULA)+PA(MULA))*DA(NUSI)
     *                       +  DA(MULA)          *PA(NUSI)
     *                       + (DB(MUSI)+PB(MUSI))*DB(LANU)
     *                       +  DB(MUSI)          *PB(LANU)
     *                       + (DB(MULA)+PB(MULA))*DB(NUSI)
     *                       +  DB(MULA)          *PB(NUSI)
C
                        DT2 =  V(NMU,NLA)*V(NNU,NSI)  ! FOR (MU NU|SI LA)
     *                       + V(NLA,NMU)*V(NSI,NNU)  ! FOR (LA SI|NU MU)
     *                       + V(NMU,NSI)*V(NNU,NLA)  ! FOR (MU NU|SI LA)
     *                       + V(NSI,NMU)*V(NLA,NNU)  ! FOR (LA SI|NU MU)
C
                        IF(LCFLAG.or.(CAMFLAG.and.LRINT)) THEN
                           IF(LRINT) DF1= -QFS*(DQ1 + DT2 + DT2 )
                        ELSE
                          DF1 = DF1 - QFS*(DQ1 + TWO*DT2)
                        ENDIF
                     ENDIF
C
                     IF(.NOT.lrint) THEN
                       IF(QFSPCP1 .NE.ZERO) THEN
                          DB1 =  CO12(NMU,NLA)*CO12(NSI,NNU)
     *                         + CO12(NMU,NSI)*CO12(NLA,NNU)
     *                         + CO12(NNU,NLA)*CO12(NSI,NMU)
     *                         + CO12(NNU,NSI)*CO12(NLA,NMU)
     *                         + CO12(NSI,NNU)*CO12(NMU,NLA)
     *                         + CO12(NLA,NNU)*CO12(NMU,NSI)
     *                         + CO12(NSI,NMU)*CO12(NNU,NLA)
     *                         + CO12(NLA,NMU)*CO12(NNU,NSI)
C
                          DF1 = DF1 + SGNK*QFSPCP1*DB1
                       ENDIF
C
                       IF(QFSPCP2 .NE.ZERO) THEN
                          DB2 =  O21V(NMU,NLA)*O21V(NSI,NNU)
     *                         + O21V(NMU,NSI)*O21V(NLA,NNU)
     *                         + O21V(NNU,NLA)*O21V(NSI,NMU)
     *                         + O21V(NNU,NSI)*O21V(NLA,NMU)
     *                         + O21V(NSI,NNU)*O21V(NMU,NLA)
     *                         + O21V(NLA,NNU)*O21V(NMU,NSI)
     *                         + O21V(NSI,NMU)*O21V(NNU,NLA)
     *                         + O21V(NLA,NMU)*O21V(NNU,NSI)
C
                          DF1 = DF1 - SGNK*QFSPCP2*DB2
                       ENDIF
C
                       IF(QFSPCP3 .NE.ZERO) THEN
                          DC1 =  BCO1(NMU,NLA)*BO2V(NNU,NSI)
     *                         + BCO1(NMU,NSI)*BO2V(NNU,NLA)
     *                         + BCO1(NNU,NLA)*BO2V(NMU,NSI)
     *                         + BCO1(NNU,NSI)*BO2V(NMU,NLA)
     *                         + BCO1(NSI,NNU)*BO2V(NLA,NMU)
     *                         + BCO1(NLA,NNU)*BO2V(NSI,NMU)
     *                         + BCO1(NSI,NMU)*BO2V(NLA,NNU)
     *                         + BCO1(NLA,NMU)*BO2V(NSI,NNU)
C
                          DC2 =  BCO2(NMU,NLA)*BO1V(NNU,NSI)
     *                         + BCO2(NMU,NSI)*BO1V(NNU,NLA)
     *                         + BCO2(NNU,NLA)*BO1V(NMU,NSI)
     *                         + BCO2(NNU,NSI)*BO1V(NMU,NLA)
     *                         + BCO2(NSI,NNU)*BO1V(NLA,NMU)
     *                         + BCO2(NLA,NNU)*BO1V(NSI,NMU)
     *                         + BCO2(NSI,NMU)*BO1V(NLA,NNU)
     *                         + BCO2(NLA,NMU)*BO1V(NSI,NNU)
C
                          DC3 =  BO2V(NMU,NLA)*BCO1(NNU,NSI)
     *                         + BO2V(NMU,NSI)*BCO1(NNU,NLA)
     *                         + BO2V(NNU,NLA)*BCO1(NMU,NSI)
     *                         + BO2V(NNU,NSI)*BCO1(NMU,NLA)
     *                         + BO2V(NSI,NNU)*BCO1(NLA,NMU)
     *                         + BO2V(NLA,NNU)*BCO1(NSI,NMU)
     *                         + BO2V(NSI,NMU)*BCO1(NLA,NNU)
     *                         + BO2V(NLA,NMU)*BCO1(NSI,NNU)
C
                          DC4 =  BO1V(NMU,NLA)*BCO2(NNU,NSI)
     *                         + BO1V(NMU,NSI)*BCO2(NNU,NLA)
     *                         + BO1V(NNU,NLA)*BCO2(NMU,NSI)
     *                         + BO1V(NNU,NSI)*BCO2(NMU,NLA)
     *                         + BO1V(NSI,NNU)*BCO2(NLA,NMU)
     *                         + BO1V(NLA,NNU)*BCO2(NSI,NMU)
     *                         + BO1V(NSI,NMU)*BCO2(NLA,NNU)
     *                         + BO1V(NLA,NMU)*BCO2(NSI,NNU)
C
                          DD1 =  BCO1(NMU,NNU)*BO2V(NSI,NLA)
     *                         + BCO1(NMU,NNU)*BO2V(NLA,NSI)
     *                         + BCO1(NNU,NMU)*BO2V(NSI,NLA)
     *                         + BCO1(NNU,NMU)*BO2V(NLA,NSI)
     *                         + BCO1(NSI,NLA)*BO2V(NMU,NNU)
     *                         + BCO1(NLA,NSI)*BO2V(NMU,NNU)
     *                         + BCO1(NSI,NLA)*BO2V(NNU,NMU)
     *                         + BCO1(NLA,NSI)*BO2V(NNU,NMU)
C
                          DD2 =  BCO2(NMU,NNU)*BO1V(NSI,NLA)
     *                         + BCO2(NMU,NNU)*BO1V(NLA,NSI)
     *                         + BCO2(NNU,NMU)*BO1V(NSI,NLA)
     *                         + BCO2(NNU,NMU)*BO1V(NLA,NSI)
     *                         + BCO2(NSI,NLA)*BO1V(NMU,NNU)
     *                         + BCO2(NLA,NSI)*BO1V(NMU,NNU)
     *                         + BCO2(NSI,NLA)*BO1V(NNU,NMU)
     *                         + BCO2(NLA,NSI)*BO1V(NNU,NMU)
C
                          DD3 =  BO2V(NMU,NNU)*BCO1(NSI,NLA)
     *                         + BO2V(NMU,NNU)*BCO1(NLA,NSI)
     *                         + BO2V(NNU,NMU)*BCO1(NSI,NLA)
     *                         + BO2V(NNU,NMU)*BCO1(NLA,NSI)
     *                         + BO2V(NSI,NLA)*BCO1(NMU,NNU)
     *                         + BO2V(NLA,NSI)*BCO1(NMU,NNU)
     *                         + BO2V(NSI,NLA)*BCO1(NNU,NMU)
     *                         + BO2V(NLA,NSI)*BCO1(NNU,NMU)
C
                          DD4 =  BO1V(NMU,NNU)*BCO2(NSI,NLA)
     *                         + BO1V(NMU,NNU)*BCO2(NLA,NSI)
     *                         + BO1V(NNU,NMU)*BCO2(NSI,NLA)
     *                         + BO1V(NNU,NMU)*BCO2(NLA,NSI)
     *                         + BO1V(NSI,NLA)*BCO2(NMU,NNU)
     *                         + BO1V(NLA,NSI)*BCO2(NMU,NNU)
     *                         + BO1V(NSI,NLA)*BCO2(NNU,NMU)
     *                         + BO1V(NLA,NSI)*BCO2(NNU,NMU)

                          DF1  = DF1 + SGNK*QFSPCP3*(-DC1-DC2-DC3-DC4
     *                                              +DD1+DD2+DD3+DD4)
                       ENDIF
                     ENDIF

                     DF1 = DF1*Q4
                     IF(DABMAX.LT. ABS(DF1)) DABMAX = ABS(DF1)
                     E(I,J,K,L)= DF1
                     IF(OUT) WRITE(IW,9010) II,JJ,KK,LL,I,J,K,L,DF1
                  END DO
               END DO
            END DO
         END DO
      ELSE
C
C D AND HIGHER FUNCTIONS OR HONDO ONLY RUNS
C
         MAXI = KMAX(II)
         MAXJ = KMAX(JJ)
         MAXK = KMAX(KK)
         MAXL = KMAX(LL)
         DO I = MINI, MAXI
            P1I = PNRM(I)
            NLA = LOCI + I
            JMAX = MAXJ
            IF (IIEQJJ) JMAX = I
            DO J = MINJ, JMAX
               P2J = P1I*PNRM(J)
               NSI = LOCJ + J
               LASI = IA(MAX0(NLA,NSI)) + MIN0(NLA,NSI)
               KKMAX = MAXK
               IF (IJEQKL) KKMAX = I
               DO K = MINK, KKMAX
                  P3K = P2J*PNRM(K)
                  NMU = LOCK + K
                  LMAX = MAXL
                  IF (KKEQLL) LMAX = K
                  IF (IJEQKL .AND. K.EQ.I) LMAX = J
                  DO L = MINL, LMAX
                     P4L = P3K*PNRM(L)
                     NNU = LOCL + L
                     MUNU = IA(MAX0(NMU,NNU)) + MIN0(NMU,NNU)
C
C                     WRITE(IW,*) 'NLA,NSI,NMU,NNU,LASI,MUNU=',
C     *                            NLA,NSI,NMU,NNU,LASI,MUNU
C
                     MUSI = IA(MAX0(NMU,NSI)) + MIN0(NMU,NSI)
                     LANU = IA(MAX0(NLA,NNU)) + MIN0(NLA,NNU)
                     MULA = IA(MAX0(NMU,NLA)) + MIN0(NMU,NLA)
                     NUSI = IA(MAX0(NNU,NSI)) + MIN0(NNU,NSI)
C
C                    WRITE(IW,*) 'MUSI,LANU,MULA,NUSI=',
C     *                              MUSI,LANU,MULA,NUSI
C
                     DF1 =  (DA(MUNU) + PA(MUNU))*DA(LASI)
     *                    +  DA(MUNU)            *PA(LASI)
C
                     IF(QFS .NE.ZERO) THEN
                        DQ1 =  (DA(MUSI)+PA(MUSI))*DA(LANU)
     *                       +  DA(MUSI)          *PA(LANU)
     *                       + (DA(MULA)+PA(MULA))*DA(NUSI)
     *                       +  DA(MULA)          *PA(NUSI)
     *                       + (DB(MUSI)+PB(MUSI))*DB(LANU)
     *                       +  DB(MUSI)          *PB(LANU)
     *                       + (DB(MULA)+PB(MULA))*DB(NUSI)
     *                       +  DB(MULA)          *PB(NUSI)
C
                        DT2 =  V(NMU,NLA)*V(NNU,NSI)  ! FOR (MU NU|SI LA)
     *                       + V(NLA,NMU)*V(NSI,NNU)  ! FOR (LA SI|NU MU)
     *                       + V(NMU,NSI)*V(NNU,NLA)  ! FOR (MU NU|SI LA)
     *                       + V(NSI,NMU)*V(NLA,NNU)  ! FOR (LA SI|NU MU)
C
                        if(lcflag.OR.(camflag.AND.LRINT)) then
                          if(     lrint) df1= -qfs*(dq1 + dt2 + dt2 )
                        else
                          DF1 = DF1 - QFS*(DQ1 + TWO*DT2)
                        endif
                     ENDIF
C
                     IF(.NOT.lrint) THEN
                       IF(QFSPCP1 .NE.ZERO) THEN
                          DB1 =  CO12(NMU,NLA)*CO12(NSI,NNU)
     *                         + CO12(NMU,NSI)*CO12(NLA,NNU)
     *                         + CO12(NNU,NLA)*CO12(NSI,NMU)
     *                         + CO12(NNU,NSI)*CO12(NLA,NMU)
     *                         + CO12(NSI,NNU)*CO12(NMU,NLA)
     *                         + CO12(NLA,NNU)*CO12(NMU,NSI)
     *                         + CO12(NSI,NMU)*CO12(NNU,NLA)
     *                         + CO12(NLA,NMU)*CO12(NNU,NSI)

                          DF1 = DF1 + SGNK*QFSPCP1*DB1
                       ENDIF
C
                       IF(QFSPCP2 .NE.ZERO) THEN
                          DB2 =  O21V(NMU,NLA)*O21V(NSI,NNU)
     *                         + O21V(NMU,NSI)*O21V(NLA,NNU)
     *                         + O21V(NNU,NLA)*O21V(NSI,NMU)
     *                         + O21V(NNU,NSI)*O21V(NLA,NMU)
     *                         + O21V(NSI,NNU)*O21V(NMU,NLA)
     *                         + O21V(NLA,NNU)*O21V(NMU,NSI)
     *                         + O21V(NSI,NMU)*O21V(NNU,NLA)
     *                         + O21V(NLA,NMU)*O21V(NNU,NSI)

                          DF1 = DF1 - SGNK*QFSPCP2*DB2
                       ENDIF
C
                       IF(QFSPCP3 .NE.ZERO) THEN
                          DC1 =  BCO1(NMU,NLA)*BO2V(NNU,NSI)
     *                         + BCO1(NMU,NSI)*BO2V(NNU,NLA)
     *                         + BCO1(NNU,NLA)*BO2V(NMU,NSI)
     *                         + BCO1(NNU,NSI)*BO2V(NMU,NLA)
     *                         + BCO1(NSI,NNU)*BO2V(NLA,NMU)
     *                         + BCO1(NLA,NNU)*BO2V(NSI,NMU)
     *                         + BCO1(NSI,NMU)*BO2V(NLA,NNU)
     *                         + BCO1(NLA,NMU)*BO2V(NSI,NNU)
C
                          DC2 =  BCO2(NMU,NLA)*BO1V(NNU,NSI)
     *                         + BCO2(NMU,NSI)*BO1V(NNU,NLA)
     *                         + BCO2(NNU,NLA)*BO1V(NMU,NSI)
     *                         + BCO2(NNU,NSI)*BO1V(NMU,NLA)
     *                         + BCO2(NSI,NNU)*BO1V(NLA,NMU)
     *                         + BCO2(NLA,NNU)*BO1V(NSI,NMU)
     *                         + BCO2(NSI,NMU)*BO1V(NLA,NNU)
     *                         + BCO2(NLA,NMU)*BO1V(NSI,NNU)
C
                          DC3 =  BO2V(NMU,NLA)*BCO1(NNU,NSI)
     *                         + BO2V(NMU,NSI)*BCO1(NNU,NLA)
     *                         + BO2V(NNU,NLA)*BCO1(NMU,NSI)
     *                         + BO2V(NNU,NSI)*BCO1(NMU,NLA)
     *                         + BO2V(NSI,NNU)*BCO1(NLA,NMU)
     *                         + BO2V(NLA,NNU)*BCO1(NSI,NMU)
     *                         + BO2V(NSI,NMU)*BCO1(NLA,NNU)
     *                         + BO2V(NLA,NMU)*BCO1(NSI,NNU)
C
                          DC4 =  BO1V(NMU,NLA)*BCO2(NNU,NSI)
     *                         + BO1V(NMU,NSI)*BCO2(NNU,NLA)
     *                         + BO1V(NNU,NLA)*BCO2(NMU,NSI)
     *                         + BO1V(NNU,NSI)*BCO2(NMU,NLA)
     *                         + BO1V(NSI,NNU)*BCO2(NLA,NMU)
     *                         + BO1V(NLA,NNU)*BCO2(NSI,NMU)
     *                         + BO1V(NSI,NMU)*BCO2(NLA,NNU)
     *                         + BO1V(NLA,NMU)*BCO2(NSI,NNU)
C
                          DD1 =  BCO1(NMU,NNU)*BO2V(NSI,NLA)
     *                         + BCO1(NMU,NNU)*BO2V(NLA,NSI)
     *                         + BCO1(NNU,NMU)*BO2V(NSI,NLA)
     *                         + BCO1(NNU,NMU)*BO2V(NLA,NSI)
     *                         + BCO1(NSI,NLA)*BO2V(NMU,NNU)
     *                         + BCO1(NLA,NSI)*BO2V(NMU,NNU)
     *                         + BCO1(NSI,NLA)*BO2V(NNU,NMU)
     *                         + BCO1(NLA,NSI)*BO2V(NNU,NMU)
C
                          DD2 =  BCO2(NMU,NNU)*BO1V(NSI,NLA)
     *                         + BCO2(NMU,NNU)*BO1V(NLA,NSI)
     *                         + BCO2(NNU,NMU)*BO1V(NSI,NLA)
     *                         + BCO2(NNU,NMU)*BO1V(NLA,NSI)
     *                         + BCO2(NSI,NLA)*BO1V(NMU,NNU)
     *                         + BCO2(NLA,NSI)*BO1V(NMU,NNU)
     *                         + BCO2(NSI,NLA)*BO1V(NNU,NMU)
     *                         + BCO2(NLA,NSI)*BO1V(NNU,NMU)
C
                          DD3 =  BO2V(NMU,NNU)*BCO1(NSI,NLA)
     *                         + BO2V(NMU,NNU)*BCO1(NLA,NSI)
     *                         + BO2V(NNU,NMU)*BCO1(NSI,NLA)
     *                         + BO2V(NNU,NMU)*BCO1(NLA,NSI)
     *                         + BO2V(NSI,NLA)*BCO1(NMU,NNU)
     *                         + BO2V(NLA,NSI)*BCO1(NMU,NNU)
     *                         + BO2V(NSI,NLA)*BCO1(NNU,NMU)
     *                         + BO2V(NLA,NSI)*BCO1(NNU,NMU)
C
                          DD4 =  BO1V(NMU,NNU)*BCO2(NSI,NLA)
     *                         + BO1V(NMU,NNU)*BCO2(NLA,NSI)
     *                         + BO1V(NNU,NMU)*BCO2(NSI,NLA)
     *                         + BO1V(NNU,NMU)*BCO2(NLA,NSI)
     *                         + BO1V(NSI,NLA)*BCO2(NMU,NNU)
     *                         + BO1V(NLA,NSI)*BCO2(NMU,NNU)
     *                         + BO1V(NSI,NLA)*BCO2(NNU,NMU)
     *                         + BO1V(NLA,NSI)*BCO2(NNU,NMU)

                          DF1 = DF1 + SGNK*QFSPCP3*(-DC1-DC2-DC3-DC4
     *                                              +DD1+DD2+DD3+DD4)
                       ENDIF
                     ENDIF
C
                     DF1= DF1*FOUR
                     IF(NMU .EQ.NNU ) DF1= DF1*HALF
                     IF(NLA .EQ.NSI ) DF1= DF1*HALF
                     IF(MUNU.EQ.LASI) DF1= DF1*HALF
C
C                     WRITE(IW,*) '** DFAC=',DF1
C
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
                     IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                     DAB(IJKL)= DF1*P4L
                     IF(OUT) WRITE(IW,9020) II,JJ,KK,LL,I,J,K,L,IJKL,DF1
                  END DO
               END DO
            END DO
         END DO
      END IF
      RETURN
 9010 FORMAT(' -DABMRSFDFT,POPLE- ',4I4,4I3,D20.12)
 9020 FORMAT(' -DABMRSFDFT,HONDO- ',4I4,4I3,I5,D20.12)
      END SUBROUTINE DABMRSFDFT

C*MODULE MRSF   *DECK DABMRSFQDFT
C>
C>    @brief    clone of DABSFDFT
C>              for quintet MRSF
C>
C>    @author   Seunghoon Lee
C>
      SUBROUTINE DABMRSFQDFT(II,JJ,KK,LL,DA,DB,PA,PB,V,DAB,DABMAX,L1,L2,&
     &                       Q4,POPLE)
C
      USE mx_limits, ONLY: mxsh, mxgtot, mxao
      USE prec, ONLY: dp
      USE constants, ONLY: ZERO, quarter, half, four
      IMPLICIT NONE
C
C DA: TOTAL (ALPHA+BETA) DENSITY FOR REFERENCE STATE
C DB: SPIN  (ALPHA-BETA) DENSITY FOR REFERENCE STATE
C PA: TOTAL (ALPHA+BETA) SF-CIS DIFFERENCE DENSITY
C PB: SPIN  (ALOHA-BETA) SF-CIS DIFFERENCE DENSITY
C  V: SF-CIS TRANSITION DENSITY
C
      REAL(KIND=dp), DIMENSION(MXGTOT) :: CCG, CD, CF, CH, CI, CP, CS,  &
     &                                    EX
      LOGICAL :: DBUG, IIEQJJ, IJEQKL, IJGTKL, IJLTKL, KKEQLL, OUT, SOME
      REAL(KIND=dp), DIMENSION(20) :: DFTTYP
      REAL(KIND=dp) :: EXENA, EXENB, EXENC
      INTEGER, DIMENSION(MXAO) :: IA
      INTEGER :: IDAF, IDFT34, IJK, IP, IPK, IR, IW, LA, LB, LC, LD,    &
     &           NAUXFUN, NAUXSHL, NAV, NSHELL
      INTEGER, DIMENSION(4,35) :: IGXYZ, JGXYZ, KGXYZ, LGXYZ
      INTEGER, DIMENSION(950) :: IODA
      INTEGER, DIMENSION(MXSH) :: KATOM, KLOC, KMAX, KMIN, KNG, KSTART, &
     &                            KTYPE
      REAL(KIND=dp), DIMENSION(84) :: PNRM
      COMMON /DFTPAR/ DFTTYP, EXENA, EXENB, EXENC, IDFT34, NAUXFUN,     &
     &                NAUXSHL
      COMMON /IJPAIR/ IA
      COMMON /INDD80/ LA, LB, LC, LD
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA
      COMMON /MP2PRT/ SOME, OUT, DBUG
      COMMON /NSHEL / EX, CS, CP, CD, CF, CCG, CH, CI, KSTART, KATOM,   &
     &                KTYPE, KNG, KLOC, KMIN, KMAX, NSHELL
      COMMON /SHLEQU/ IIEQJJ, KKEQLL, IJEQKL, IJGTKL, IJLTKL
      COMMON /SHLLMN/ IGXYZ, JGXYZ, KGXYZ, LGXYZ
      COMMON /SHLNRM/ PNRM
C
      REAL(KIND=dp) :: DABMAX, Q4
      INTEGER :: II, JJ, KK, L1, L2, LL
      LOGICAL :: POPLE
      REAL(KIND=dp), DIMENSION(L2) :: DA, DB, PA, PB
      REAL(KIND=dp), DIMENSION(*) :: DAB
      REAL(KIND=dp), DIMENSION(L1,*) :: V
C
      REAL(KIND=dp) :: DF1, DQ1, DT2, P1I, P2J, P3K, P4L, QFS
      INTEGER :: I, IJKL, J, JMAX, K, KKMAX, L, LANU, LASI, LMAX, LOCI, &
     &           LOCJ, LOCK, LOCL, MAXI, MAXJ, MAXK, MAXL, MINI, MINJ,  &
     &           MINK, MINL, MULA, MUNU, MUSI, NLA, NMU, NNU, NSI, NUSI
C
CJMS  LABELLED COMMON GSPG80 DEFINED FOR COMPUTATIONAL EFFICIENCY.
CJMS  FOR SP BASES ONLY, IT CONTAINS THE E ARRAY WHICH IS THE DAB
CJMS  ARRAY WITH INDICES IN REVERSE ORDER: E(I,J,K,L)= DAB(L,K,J,I)
CJMS  AND IS USED IN SUB JKDG80 (MODULE GRD2B). IT ORIGINATES IN:
CJMS
CJMS     1. SUBS DABCLU, DABDFT, DABGVB, DABMC, DABMP2 AND DABUMP
CJMS        (MODULE GRD2A) AND SUB DABPAU (MODULE EFPAUL) WHICH ARE
CJMS        ALL CALLED BY SUB JKDER (MODULE GRD2A)
CJMS
CJMS     2. SUB DABCLU (MODULE GRD2A) WHICH IS CALLED BY SUB EFDEN OF
CJMS        MODULE EFGRD2
CJMS
CJMS     3. SUB PAR2PDM (MODULE MP2DDI) WHICH IS CALLED BY SUB PJKDMP2
CJMS        OF MODULE MP2DDI
C
      REAL(KIND=dp), DIMENSION(4,4,4,4) :: E
      COMMON /GSPG80/ E
C
      QFS = DFTTYP(3)*QUARTER
C
C     ----- FORM TWO-PARTICLE DENSITY MATRIX FOR CIS GRADIENT -----
C
      DABMAX = ZERO
      MINI = KMIN(II)
      MINJ = KMIN(JJ)
      MINK = KMIN(KK)
      MINL = KMIN(LL)
      LOCI = KLOC(II) - MINI
      LOCJ = KLOC(JJ) - MINJ
      LOCK = KLOC(KK) - MINK
      LOCL = KLOC(LL) - MINL
C
      IF (POPLE) THEN
         DO L = 1, LD
            NNU = LOCL + L
            DO K = 1, LC
               NMU = LOCK + K
               MUNU = IA(MAX0(NMU,NNU)) + MIN0(NMU,NNU)
               DO J = 1, LB
                  NSI = LOCJ + J
                  DO I = 1, LA
                     NLA = LOCI + I
                     LASI = IA(MAX0(NLA,NSI)) + MIN0(NLA,NSI)
C
C                    WRITE(IW,*) 'NLA,NSI,NMU,NNU,LASI,MUNU=',
C    *                            NLA,NSI,NMU,NNU,LASI,MUNU
C
                     MUSI = IA(MAX0(NMU,NSI)) + MIN0(NMU,NSI)
                     LANU = IA(MAX0(NLA,NNU)) + MIN0(NLA,NNU)
                     MULA = IA(MAX0(NMU,NLA)) + MIN0(NMU,NLA)
                     NUSI = IA(MAX0(NNU,NSI)) + MIN0(NNU,NSI)
C
                     DF1 =  (DA(MUNU) + PA(MUNU))*DA(LASI)
     *                    +  DA(MUNU)            *PA(LASI)
C
                     IF(QFS.NE.ZERO) THEN
                        DQ1 =  (DA(MUSI)+PA(MUSI))*DA(LANU)
     *                       +  DA(MUSI)          *PA(LANU)
     *                       + (DA(MULA)+PA(MULA))*DA(NUSI)
     *                       +  DA(MULA)          *PA(NUSI)
     *                       + (DB(MUSI)+PB(MUSI))*DB(LANU)
     *                       +  DB(MUSI)          *PB(LANU)
     *                       + (DB(MULA)+PB(MULA))*DB(NUSI)
     *                       +  DB(MULA)          *PB(NUSI)
C
                        DT2 =  V(NMU,NLA)*V(NNU,NSI)  ! for (mu nu|si la)
     *                       + V(NLA,NMU)*V(NSI,NNU)  ! for (la si|nu mu)
     *                       + V(NMU,NSI)*V(NNU,NLA)  ! for (mu nu|la si)
     *                       + V(NSI,NMU)*V(NLA,NNU)  ! for (si la|nu mu)
C
                        DF1 = DF1 - QFS*(DQ1 + DT2 + DT2 )
                     END IF
C
                     DF1 = DF1*Q4
                     IF(DABMAX.LT. ABS(DF1)) DABMAX = ABS(DF1)
                     E(I,J,K,L)= DF1
                     IF(OUT) WRITE(IW,9010) II,JJ,KK,LL,I,J,K,L,DF1
                  END DO
               END DO
            END DO
         END DO
      ELSE
C
C D AND HIGHER FUNCTIONS OR HONDO ONLY RUNS
C
         MAXI = KMAX(II)
         MAXJ = KMAX(JJ)
         MAXK = KMAX(KK)
         MAXL = KMAX(LL)
         DO I = MINI, MAXI
            P1I = PNRM(I)
            NLA = LOCI + I
            JMAX = MAXJ
            IF (IIEQJJ) JMAX = I
            DO J = MINJ, JMAX
               P2J = P1I*PNRM(J)
               NSI = LOCJ + J
               LASI = IA(MAX0(NLA,NSI)) + MIN0(NLA,NSI)
               KKMAX = MAXK
               IF (IJEQKL) KKMAX = I
               DO K = MINK, KKMAX
                  P3K = P2J*PNRM(K)
                  NMU = LOCK + K
                  LMAX = MAXL
                  IF (KKEQLL) LMAX = K
                  IF (IJEQKL .AND. K==I) LMAX = J
                  DO L = MINL, LMAX
                     P4L = P3K*PNRM(L)
                     NNU = LOCL + L
                     MUNU = IA(MAX0(NMU,NNU)) + MIN0(NMU,NNU)
C
C                     WRITE(IW,*) 'NLA,NSI,NMU,NNU,LASI,MUNU=',
C     *                            NLA,NSI,NMU,NNU,LASI,MUNU
C
                     MUSI = IA(MAX0(NMU,NSI)) + MIN0(NMU,NSI)
                     LANU = IA(MAX0(NLA,NNU)) + MIN0(NLA,NNU)
                     MULA = IA(MAX0(NMU,NLA)) + MIN0(NMU,NLA)
                     NUSI = IA(MAX0(NNU,NSI)) + MIN0(NNU,NSI)
C
C                    WRITE(IW,*) 'MUSI,LANU,MULA,NUSI=',
C     *                              MUSI,LANU,MULA,NUSI
C
                     DF1 =  (DA(MUNU) + PA(MUNU))*DA(LASI)
     *                    +  DA(MUNU)            *PA(LASI)
C
                     IF(QFS.NE.ZERO) THEN
C
                        DQ1 =  (DA(MUSI)+PA(MUSI))*DA(LANU)
     *                       +  DA(MUSI)          *PA(LANU)
     *                       + (DA(MULA)+PA(MULA))*DA(NUSI)
     *                       +  DA(MULA)          *PA(NUSI)
     *                       + (DB(MUSI)+PB(MUSI))*DB(LANU)
     *                       +  DB(MUSI)          *PB(LANU)
     *                       + (DB(MULA)+PB(MULA))*DB(NUSI)
     *                       +  DB(MULA)          *PB(NUSI)
C
                        DT2 =  V(NMU,NLA)*V(NNU,NSI)
     *                       + V(NLA,NMU)*V(NSI,NNU)
     *                       + V(NMU,NSI)*V(NNU,NLA)
     *                       + V(NSI,NMU)*V(NLA,NNU)
C
                        DF1 = DF1 - QFS*(DQ1 + DT2 + DT2)
                     END IF
                     DF1= DF1*FOUR
                     IF(NMU .EQ.NNU ) DF1= DF1*HALF
                     IF(NLA .EQ.NSI ) DF1= DF1*HALF
                     IF(MUNU.EQ.LASI) DF1= DF1*HALF
C
C                     WRITE(IW,*) '** DFAC=',DF1
C
                     IF(DABMAX.LT. ABS(DF1)) DABMAX= ABS(DF1)
                     IJKL=IGXYZ(1,I)+JGXYZ(1,J)+KGXYZ(1,K)+LGXYZ(1,L)
                     DAB(IJKL)= DF1*P4L
                     IF(OUT) WRITE(IW,9020) II,JJ,KK,LL,I,J,K,L,IJKL,DF1
                  END DO
               END DO
            END DO
         END DO
      END IF
      RETURN
 9010 FORMAT(' -DABSFDFT,POPLE- ',4I4,4I3,D20.12)
 9020 FORMAT(' -DABSFDFT,HONDO- ',4I4,4I3,I5,D20.12)
      END SUBROUTINE DABMRSFQDFT
C
C*MODULE MRSF  *DECK MRSFQRORHS
C>
C>    @brief    clone of SFRORHS
C>              for quintet MRSF
C>
C>    @author   Seunghoon Lee
C>
      SUBROUTINE MRSFQRORHS(R,HPTA,HPTB,XHXA,XHXB,TA,TB,FA,FB,SCR,NOCA, &
     &                      NOCB,NVIRA,NVIRB,LX,LZDIM)
C
      USE prec, ONLY: dp
      USE constants, ONLY: one, two
      IMPLICIT NONE
C
      INTEGER :: LX, LZDIM, NOCA, NOCB, NVIRA, NVIRB
      REAL(KIND=dp), DIMENSION(LX,LX) :: FA, FB, SCR, XHXA
      REAL(KIND=dp), DIMENSION(NOCA,NVIRA) :: HPTA
      REAL(KIND=dp), DIMENSION(NOCB,NVIRB) :: HPTB
      REAL(KIND=dp), DIMENSION(*) :: R
      REAL(KIND=dp), DIMENSION(NVIRA,NVIRA) :: TA
      REAL(KIND=dp), DIMENSION(NOCB,NOCB) :: TB
      REAL(KIND=dp), DIMENSION(LX,NOCB) :: XHXB
C
      INTEGER :: I, IJ, J, MA, MI, MX, NCONF
C
C     ----- ALPHA -----
C
C XHXA+= 2*FA(P+,A+)*TA(A+,B+)
      CALL VCLR(SCR,1,LX*LX)
      DO J=NOCA+1,LX
         DO I=NOCA+1,LX
            SCR(I,J) = TA(I-NOCA,J-NOCA)
         END DO
      END DO
      CALL DGEMM('N','N',LX,LX,LX,TWO,FA,LX,SCR,LX,ONE,XHXA,LX)
C
C     ----- BETA -----
C
C XHXB+= 2*FB(P-,I-)*TB(I-,J-)
      CALL DGEMM('N','N',LX,NOCB,NOCB,TWO,FB,LX,TB,NOCB,ONE,XHXB,LX)
C
c      CALL VCLR(R,1,LZDIM)
C
C     ----- DOC-SOCC -----
C
      IJ = 0
      DO MX=NOCB+1,NOCA
         DO MI=1,NOCB
            IJ = IJ + 1
            R(IJ) = HPTB(MI,MX-NOCB) + XHXB(MX,MI)
         END DO
      END DO
C
C     ----- DOC-VIRT -----
C
      DO MA=NOCA+1,LX
         DO MI=1,NOCB
            IJ = IJ + 1
            R(IJ) = HPTA(MI,MA-NOCA) + HPTB(MI,MA-NOCB)
     *             + XHXB(MA,MI) - XHXA(MI,MA)
         END DO
      END DO
C
C     ----- SOC-VIRT -----
C
      DO MA=NOCA+1,LX
         DO MX=NOCB+1,NOCA
            IJ = IJ + 1
            R(IJ) = HPTA(MX,MA-NOCA) - XHXA(MX,MA)
         END DO
      END DO
C
C     ----- MULTIPLIED BY -1 I.E., RHS OF Z-VECTOR EQ. -----
C
      NCONF = IJ
      CALL DSCAL(NCONF,-1.0D+00,R,1)
C
      RETURN
      END SUBROUTINE MRSFQRORHS
c
C*MODULE MRSF   *DECK MRSFQROPCAL
C>
C>    @brief    clone of SFROPCAL
C>              for quintet MRSF
C>
C>    @author   Seunghoon Lee
C>
      SUBROUTINE MRSFQROPCAL(PA,PB,ZA,ZB,TA,TB,Z,LX,NOCA,NOCB,NVIRA,    &
     &                       NVIRB)
C
      USE prec, ONLY: dp
      USE constants, ONLY: half
      IMPLICIT NONE
C
      INTEGER :: LX, NOCA, NOCB, NVIRA, NVIRB
      REAL(KIND=dp), DIMENSION(LX,LX) :: PA, PB, ZA, ZB
      REAL(KIND=dp), DIMENSION(NVIRA,NVIRA) :: TA
      REAL(KIND=dp), DIMENSION(NOCB,NOCB) :: TB
      REAL(KIND=dp), DIMENSION(*) :: Z
C
      INTEGER :: I, IJ, J, MA, MI, MX
C
C     ----- COPY T -----
C
C ALPHA
      CALL VCLR(PA,1,LX*LX)
      DO J=NOCA+1,LX
         DO I=NOCA+1,LX
            PA(I,J) = TA(I-NOCA,J-NOCA)
         END DO
      END DO
C BETA
      CALL VCLR(PB,1,LX*LX)
      DO J=1,NOCB
         DO I=1,NOCB
            PB(I,J) = TB(I,J)
         END DO
      END DO
C
      CALL VCLR(ZA,1,LX*LX)
      CALL VCLR(ZB,1,LX*LX)
C
C     ----- ADD Z CONTRIBUTION -----
C
C DOC-SOCC
C
      IJ = 0
      DO MX=NOCB+1,NOCA
         DO MI=1,NOCB
            IJ = IJ + 1
            PB(MI,MX) = PB(MI,MX) + HALF*Z(IJ)
            ZB(MI,MX) = HALF*Z(IJ)
         END DO
      END DO
C
C DOC-VIRT
C
      DO MA=NOCA+1,LX
         DO MI=1,NOCB
            IJ = IJ + 1
            PA(MI,MA) = PA(MI,MA) + HALF*Z(IJ)
            PB(MI,MA) = PB(MI,MA) + HALF*Z(IJ)
            ZA(MI,MA) = HALF*Z(IJ)
            ZB(MI,MA) = HALF*Z(IJ)
         END DO
      END DO
C
C SOCC-VIRT
C
      DO MA=NOCA+1,LX
         DO MX=NOCB+1,NOCA
            IJ = IJ + 1
            PA(MX,MA) = PA(MX,MA) + HALF*Z(IJ)
            ZA(MX,MA) = HALF*Z(IJ)
         END DO
      END DO
C
      RETURN
      END SUBROUTINE MRSFQROPCAL
c
c*MODULE MRSF   *deck mrsfrowcal
C>
C>    @brief    clone of SFROWCAL
C>              for singlet and triplet MRSF
C>
C>    @details  Eq. (III.18) of JCP 150 184111 (2019)
C>
C>    @author   Seunghoon Lee
C>
      subroutine mrsfrowcal(w,ee,eorb,fa,fb,x12,z,wrk1,wrk2,wrk3,xhxa,  &
     &                      xhxb,hppija,hppijb,lx,noca,nocb,nvirb)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, half
      implicit none
c
      real(kind=dp) :: ee
      integer :: lx, noca, nocb, nvirb
      real(kind=dp), dimension(*) :: eorb, z
      real(kind=dp), dimension(lx,lx) :: fa, fb, w, wrk1, wrk2, xhxb
      real(kind=dp), dimension(noca,noca) :: hppija
      real(kind=dp), dimension(nocb,nocb) :: hppijb
      real(kind=dp), dimension(noca,*) :: wrk3
      real(kind=dp), dimension(noca,nvirb) :: x12
      real(kind=dp), dimension(lx,noca) :: xhxa
c
      real(kind=dp) :: dum
      integer :: i, ij, ma, mb, mc, mi, mj, mk, mw, mx, my
c
c     ----- copy z -----
c
      call vclr(wrk1,1,lx*lx)
c
      ij = 0
      do mx=nocb+1,noca
         do mi=1,nocb
            ij = ij + 1
            wrk1(mi,mx) = z(ij)
         end do
      end do
c
      do ma=noca+1,lx
         do mi=1,nocb
            ij = ij + 1
            wrk1(mi,ma) = z(ij)
         end do
      end do
c
      do ma=noca+1,lx
         do mx=nocb+1,noca
            ij = ij + 1
            wrk1(mx,ma) = z(ij)
         end do
      end do
c
c     ----- w_ix -----
c
      do mx=nocb+1,noca
         do mi=1,nocb
            dum = zero
            do mk=1,nocb
               dum = dum - fa(mk,mi)*wrk1(mk,mx)
c               dum = dum + fa(mk,mi)*wrk1(mk,mx)
            end do
            do mc=noca+1,lx
               dum = dum + fa(mc,mi)*wrk1(mx,mc)
            end do
            w(mi,mx) = eorb(mi)*wrk1(mi,mx) + half*dum
     *                 + xhxa(mi,mx) + xhxb(mi,mx)
     *                 + hppija(mi,mx)
         end do
      end do
c
c     ----- w_ia -----
c
      do ma=noca+1,lx
         do mi=1,nocb
            dum = zero
            do mw=nocb+1,noca
               dum = dum + fa(mw,mi)*wrk1(mw,ma)
            end do
            w(mi,ma) = eorb(mi)*wrk1(mi,ma) + half*dum
     *                 + xhxb(mi,ma)
         end do
      end do
c
c     ----- w_xa -----
c
      do ma=noca+1,lx
         do mx=nocb+1,noca
            dum = zero
            do mk=1,nocb
               dum = dum + fa(mk,mx)*wrk1(mk,ma)
            end do
            do mw=nocb+1,noca
               dum = dum - fb(mw,mx)*wrk1(mw,ma)
            end do
            w(mx,ma) = eorb(mx)*wrk1(mx,ma) + half*dum
     *                 + xhxb(mx,ma)
         end do
      end do
C
c     ----- w_ij -----
c
      do mi=1,nocb
      do mj=1,mi
         w(mi,mj) = hppija(mi,mj) + hppijb(mi,mj) + xhxa(mj,mi)
      end do
      end do
c
c     ----- w_xy -----
c
      do mx=nocb+1,noca
      do my=nocb+1,mx
         w(mx,my) = xhxa(my,mx) + xhxb(my,mx)
     *            + hppija(mx,my)
      end do
      end do
c
c     ----- w_ab -----
c
      do ma=noca+1,lx
         do mb=noca+1,ma
            w(ma,mb) = xhxb(mb,ma)
         end do
      end do
c
c     ----- scale diagonal elements -----
c
      do i=1,lx
         w(i,i) = half*w(i,i)
      end do
c
      call dscal(lx*lx,-1.0d+00,w,1)
c
      return
      end subroutine mrsfrowcal
c
c*MODULE MRSF   *deck mrsfqrowcal
C>
C>    @brief    clone of SFROWCAL
C>              for quintet MRSF
C>
C>    @author   Seunghoon Lee
C>
      subroutine mrsfqrowcal(w,eorb,fa,fb,z,wrk1,wrk2,wrk3,xhxa,xhxb,   &
     &                       hppija,hppijb,hpzija,hpzijb,lx,noca,nocb,  &
     &                       nvira,nvirb)
C
      USE prec, ONLY: dp
      USE constants, ONLY: zero, half
      implicit none
c
      integer :: lx, noca, nocb, nvira, nvirb
      real(kind=dp), dimension(*) :: eorb, z
      real(kind=dp), dimension(lx,lx) :: fa, fb, w, wrk1, wrk2, xhxa
      real(kind=dp), dimension(noca,noca) :: hppija, hpzija
      real(kind=dp), dimension(nocb,nocb) :: hppijb, hpzijb
      real(kind=dp), dimension(nocb,*) :: wrk3
      real(kind=dp), dimension(lx,nocb) :: xhxb
c
      real(kind=dp) :: dum
      integer :: i, ij, ma, mb, mc, mi, mj, mk, mw, mx, my
c
c     ----- copy z -----
c
      call vclr(wrk1,1,lx*lx)
c
      ij = 0
      do mx=nocb+1,noca
         do mi=1,nocb
            ij = ij + 1
            wrk1(mi,mx) = z(ij)
         end do
      end do
c
      do ma=noca+1,lx
         do mi=1,nocb
            ij = ij + 1
            wrk1(mi,ma) = z(ij)
         end do
      end do
c
      do ma=noca+1,lx
         do mx=nocb+1,noca
            ij = ij + 1
            wrk1(mx,ma) = z(ij)
         end do
      end do
c
c     ----- w_ix -----
c
      do mx=nocb+1,noca
         do mi=1,nocb
            dum = zero
            do mk=1,nocb
               dum = dum - fa(mk,mi)*wrk1(mk,mx)
c               dum = dum + fa(mk,mi)*wrk1(mk,mx)
            end do
            do mc=noca+1,lx
               dum = dum + fa(mc,mi)*wrk1(mx,mc)
            end do
            w(mi,mx) = eorb(mi)*wrk1(mi,mx) + half*dum
     *                 + hppija(mi,mx)
         end do
      end do
c
c     ----- w_ia -----
c
      do ma=noca+1,lx
         do mi=1,nocb
            dum = zero
            do mw=nocb+1,noca
               dum = dum + fa(mw,mi)*wrk1(mw,ma)
            end do
            w(mi,ma) = eorb(mi)*wrk1(mi,ma) + half*dum
     *                 + xhxa(mi,ma)
         end do
      end do
c
c     ----- w_xa -----
c
      do ma=noca+1,lx
         do mx=nocb+1,noca
            dum = zero
            do mk=1,nocb
               dum = dum + fa(mk,mx)*wrk1(mk,ma)
            end do
            do mw=nocb+1,noca
               dum = dum - fb(mw,mx)*wrk1(mw,ma)
            end do
            w(mx,ma) = eorb(mx)*wrk1(mx,ma) + half*dum
     *                 + xhxa(mx,ma)
         end do
      end do
C
c     ----- w_ij -----
c
      do mi=1,nocb
      do mj=1,mi
         w(mi,mj) = hppija(mi,mj) + hppijb(mi,mj) + xhxb(mj,mi)
      end do
      end do
c
c     ----- w_xy -----
c
      do mx=nocb+1,noca
      do my=nocb+1,mx
         w(mx,my) = hppija(mx,my)
      end do
      end do
c
c     ----- w_ab -----
c
      do ma=noca+1,lx
         do mb=noca+1,ma
            w(ma,mb) = xhxa(mb,ma)
         end do
      end do
c
c     ----- scale diagonal elements -----
c
      do i=1,lx
         w(i,i) = half*w(i,i)
      end do
c
      call dscal(lx*lx,-1.0d+00,w,1)
c
      return
      end subroutine mrsfqrowcal
c
!*MODULE MRSF  *deck mrsfxvec
!>
!>    @brief    dimensional transformed X amplitudes
!>              for singlet and triplet MRSF
!>
!>    @details  Eqs. (II.4)--(II.7) of JCP 150 184111 (2019)
!>
!>    @author   Seunghoon Lee, Konstantin Komarov
!>
      subroutine mrsfxvec(xv,xv12,noca,nocb,lx,mrsfs,mrsft)

        use prec, only: dp

        implicit none

        real(kind=dp), dimension(*) :: xv
        real(kind=dp), dimension(noca,*) :: xv12
        integer :: noca, nocb, lx
        logical :: mrsfs, mrsft

        integer :: i, ij, ijd, ijg, ijlr1, ijlr2, j
        real(kind=dp), parameter :: sqrt2 = 1/sqrt(2.0_dp)
        real(kind=dp), allocatable, dimension(:,:) :: tmp

        ijlr1 = (noca-1-nocb-1)*noca+noca-1
        ijg   = (noca-1-nocb-1)*noca+noca
        ijd   = (noca  -nocb-1)*noca+noca-1
        ijlr2 = (noca  -nocb-1)*noca+noca

        allocate(tmp(noca,lx-nocb), source=0.0_dp)

        if (mrsfs) then
          do i = 1, noca
            do j = nocb+1, lx
              ij = (j-nocb-1)*noca+i
              if (ij==ijlr1) then
                tmp(i,j-nocb) = xv(ijlr1)*sqrt2
                cycle
              else if(ij==ijlr2) then
                tmp(i,j-nocb) = -xv(ijlr1)*sqrt2
                cycle
              end if
              if (ij>ijlr2) ij = ij-1
              tmp(i,j-nocb) = xv(ij)
            end do
          end do
        else if (mrsft) then
          do i = 1, noca
            do j = nocb+1, lx
              ij = (j-nocb-1)*noca+i
              if (ij==ijlr1) then
                tmp(i,j-nocb) = xv(ijlr1)*sqrt2
                cycle
              else if (ij==ijlr2) then
                tmp(i,j-nocb) = xv(ijlr1)*sqrt2
                cycle
              else if (ij==ijg) then
                tmp(i,j-nocb) = 0.0_dp
                cycle
              else if (ij==ijd) then
                tmp(i,j-nocb) = 0.0_dp
                cycle
              end if
              if (ij>ijg .and. ij<ijd) then
                ij = ij-1
              else if (ij>ijd .and. ij<ijlr2) then
                ij = ij-2
              else if (ij>ijlr2) then
                ij = ij-3
              end if
              tmp(i,j-nocb) = xv(ij)
            end do
          end do
        end if

        xv12(1:noca,1:lx-nocb) = tmp(1:noca,1:lx-nocb)

        return

      end subroutine mrsfxvec
c
c*MODULE MRSF  *deck mrsfsp
C>
C>    @brief    Spin-pairing parts
C>              of singlet and triplet MRSF Lagrangian
C>
C>    @details  Eqs. (III.14), (III.15) of JCP 150 184111 (2019)
C>
C>    @author   Seunghoon Lee
C>
C>    @date    Nov, 2021 Initial release
C>
C>    @author  Konstantin Komarov
c>
C>    @details Terms co2v, co1v, cco1 and cco2 was replaced
C>             by o21v, co12. mrsfcbc is initial subroutine of this.
C>
C>    @date    Mar, 2022 Performance improvements
C>
      SUBROUTINE mrsfsp(xhxa,xhxb,ca,cb,xv,
     &                  ao21v,aco12,
     &                  ado2v,ado1v,adco1,adco2,
     &                  scr,scr2,l1,lx,noca,nocb,l3)

      USE prec, ONLY: dp
      USE constants, ONLY: zero, one, two
      implicit none

      integer :: l1, l3, lx, noca, nocb
      real(kind=dp), dimension(l1,*) :: adco1, adco2, ado1v, ado2v,
     &                                  aco12, ao21v
      real(kind=dp), dimension(l1,*) :: ca, cb
      real(kind=dp), dimension(lx,*) :: scr, scr2
      real(kind=dp), dimension(l1,l1) :: tmp
      real(kind=dp), dimension(lx,noca) :: xhxa
      real(kind=dp), dimension(lx,*) :: xhxb, xv

      integer :: i, i1, i2, j, j1, j2

c spin-pairing coupling contributions of xhxa

C o1v
      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,AO21V,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX) ! 1 for nocb+1
                                                                 ! 2 for noca
      i2=nocb+1
      do j=noca+1,lx
         do i1=1,lx
            xhxa(i1,i2+1)=xhxa(i1,i2+1)+scr(i1,j)*xv(i2,j)     ! 1
            xhxa(i1,i2  )=xhxa(i1,i2  )-scr(i1,j)*xv(i2+1,j)   ! 2
         enddo
      enddo

c co1
      CALL DCOPY(L3,ACO12,1,TMP,1)
      CALL DSCAL(L3,-ONE,TMP,1)
      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,TMP,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,lX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX) ! NOCA for 1
                                                                 ! NOCB+1 for 2
      j=nocb+1
      do i2=1,nocb
         do i1=1,lx
            xhxa(i1,i2)=xhxa(i1,i2)+scr(i1,j+1)*xv(i2,j)     ! 1
            xhxa(i1,i2)=xhxa(i1,i2)-scr(i1,j)  *xv(i2,j+1)   ! 2
         enddo
      enddo

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ADCO2,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      i2=nocb+1
      do j=noca+1,lx
         do i1=1,lx
            xhxa(i1,i2)=xhxa(i1,i2)+scr(i1,j)*xv(i2,j)
         enddo
      enddo

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ADCO1,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      i2=noca
      do j=noca+1,lx
         do i1=1,lx
            xhxa(i1,i2)=xhxa(i1,i2)+scr(i1,j)*xv(i2,j)
         enddo
      enddo

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ADO2V,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,1,L1,TWO,SCR2,LX,CB(1,NOCB+1),L1,
     *           ZERO,SCR,LX)

      j=nocb+1
      do i2=1,nocb
         do i1=1,lx
            xhxa(i1,i2)=xhxa(i1,i2)+scr(i1,1)*xv(i2,j)
         enddo
      enddo

C co2
      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ADO1V,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,1,L1,TWO,SCR2,LX,CB(1,NOCA),L1,
     *           ZERO,SCR,LX)

      j=noca
      do i2=1,nocb
         do i1=1,lx
            xhxa(i1,i2)=xhxa(i1,i2)+scr(i1,1)*xv(i2,j)
         enddo
      enddo

C SPIN-PAIRING COUPLING CONTRIBUTIONS OF XHXB

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,AO21V,L1,ZERO,SCR2,lX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      I=NOCB+1
      DO J1=1,LX
         DO J2=NOCA+1,LX
            XHXB(J1,J2)=XHXB(J1,J2)+SCR(I+1,J1)*XV(I,J2)
         ENDDO
      ENDDO

      CALL DCOPY(L3,AO21V,1,TMP,1)
      CALL DSCAL(L3,-ONE,TMP,1)
      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,TMP,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      I=NOCA
      DO J1=1,LX
         DO J2=NOCA+1,LX
            XHXB(J1,J2)=XHXB(J1,J2)+SCR(I-1,J1)*XV(I,J2)
         ENDDO
      ENDDO

C CO1
      CALL DCOPY(L3,ACO12,1,TMP,1)
      CALL DSCAL(L3,-ONE,TMP,1)
      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,TMP,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      j2=nocb+1
      do j1=1,lx
         do i=1,nocb
            xhxb(j1,j2+1)=xhxb(j1,j2+1)+scr(i,j1)*xv(i,j2)
         enddo
      enddo

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ACO12,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      J2=NOCA
      DO J1=1,LX
         DO I=1,NOCB
            XHXB(J1,J2-1)=XHXB(J1,J2-1)+SCR(I,J1)*XV(I,J2)
         ENDDO
      ENDDO

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ADO2V,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      j2=nocb+1
      do j1=1,lx
         do i=1,nocb
            xhxb(j1,j2)=xhxb(j1,j2)+scr(i,j1)*xv(i,j2)
         enddo
      enddo

      CALL DGEMM('T','N',LX,L1,L1,ONE,CA,L1,ADO1V,L1,ZERO,SCR2,LX)
      CALL DGEMM('N','N',LX,LX,L1,TWO,SCR2,LX,CB,L1,ZERO,SCR,LX)

      J2=NOCA
      DO J1=1,LX
         DO I=1,NOCB
            XHXB(J1,J2)=XHXB(J1,J2)+SCR(I,J1)*XV(I,J2)
         ENDDO
      ENDDO

C O1V
      CALL DGEMM('T','N',1,L1,L1,ONE,CA(1,NOCB+1),L1,ADCO2,L1,
     *           ZERO,SCR2,1)
      CALL DGEMM('N','N',1,LX,L1,TWO,SCR2,1,CB,L1,ZERO,SCR,1)

      I=NOCB+1
      DO J1=1,LX
         DO J2=NOCA+1,LX
            XHXB(J1,J2)=XHXB(J1,J2)+SCR(J1,1)*XV(I,J2)
         ENDDO
      ENDDO

      CALL DGEMM('T','N',1,L1,L1,ONE,CA(1,NOCA),L1,ADCO1,L1,
     *           ZERO,SCR2,1)
      CALL DGEMM('N','N',1,LX,L1,TWO,SCR2,1,CB,L1,ZERO,SCR,1)

      I=NOCA
      DO J1=1,LX
         DO J2=NOCA+1,LX
            XHXB(J1,J2)=XHXB(J1,J2)+SCR(J1,1)*XV(I,J2)
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE MRSFSP
C*MODULE SFDFT   *DECK UNXMDRST
      SUBROUTINE UNXMDRST

      USE comm_NONAD, ONLY: NAMD,NDSWCH,NDRST,NDTLF,COLD,THRSHE
      use mx_limits, only: mxatm, mxrt, mxao
      USE prec, only: dp

      IMPLICIT NONE

      REAL(KIND=dp), PARAMETER :: UNITS = 0.52917724924D+00
      REAL(KIND=dp), PARAMETER :: AU2SEC = 2.41889D-17

        LOGICAL :: DSKWRK, GOPARR, MASWRK
        INTEGER :: ME, MASTER, NPROC, IBTYP, IPTIM
      COMMON /PAR   / ME, MASTER, NPROC, IBTYP, IPTIM, GOPARR, DSKWRK,  &
     &                MASWRK

        INTEGER :: IDAF, IJK, IPK, IP, IR, IW, NAV
        INTEGER, DIMENSION(950) :: IODA
      COMMON /IOFILE/ IR, IW, IP, IJK, IPK, IDAF, NAV, IODA

        REAL(KIND=dp), DIMENSION(MXATM) :: ANAM, BNAM
        REAL(KIND=dp), DIMENSION(MXAO) :: BFLAB
        REAL(KIND=dp), DIMENSION(10) :: TITLE
      COMMON /RUNLAB/ TITLE, ANAM, BNAM, BFLAB

        REAL(KIND=dp), DIMENSION(3,MXATM) :: CR
        REAL(KIND=dp), DIMENSION(MXATM) :: ZAN
        INTEGER, DIMENSION(MXATM) :: IAN
        INTEGER :: ICH, MUL, NA, NAT, NB, NE, NQMT, NUM
      COMMON /INFOA / NAT, ICH, MUL, NUM, NQMT, NE, NA, NB, ZAN, CR, IAN

        REAL(KIND=dp) :: E1, E2, ECORE, EDISP, EELCT, EERD, EKIN, ENUCR,&
     &                   EPOT, ESCF, ETOT, STATN, SZ, SZZ, VEE, VEN
        REAL(KIND=dp), DIMENSION(2) :: EDFT
        REAL(KIND=dp), DIMENSION(MXRT) :: ESTATE
      COMMON /ENRGYS/ ENUCR, EELCT, ETOT, SZ, SZZ, ECORE, ESCF, EERD,   &
     &                E1, E2, VEN, VEE, EPOT, EKIN, ESTATE, STATN,EDFT, &
     &                EDISP

        LOGICAL :: ALPHKWD, BETAKWD, MRDEA, MREKT, SG1T, TAMMD, TPA,    &
     &             TRIPLET
        REAL(KIND=dp) :: CNVTOL
        INTEGER :: IRECTD, ITDFG, ITDPRP, JANST, MAXVEC, MODTD, MTHST,  &
     &             MULTD, NLEBT, NONEQR, NPHIT, NRADT, NSTAT, NTHET,    &
     &             NTHST, NTRIAL
        INTEGER, DIMENSION(4) :: IFEDAT
        REAL(KIND=dp), DIMENSION(2) :: PFREQ
        REAL(KIND=dp), DIMENSION(3) :: SPCP
      COMMON /INFOTD/ CNVTOL, PFREQ, MODTD, JANST, NRADT, NTHET, NPHIT, &
     &                NLEBT, NSTAT, NTRIAL, MAXVEC, NTHST, IRECTD,      &
     &                ITDFG, ITDPRP, TRIPLET, SG1T, NONEQR, TAMMD, TPA, &
     &                ALPHKWD, BETAKWD, SPCP, MULTD, MREKT, MRDEA,      &
     &                MTHST, IFEDAT

        REAL(KIND=dp) :: ENERGY, EG
      COMMON /FUNCT / ENERGY, EG(3,MXATM)

        INTEGER :: ISTEP, NDNINT, NTHSTOLD
        REAL(KIND=dp), DIMENSION(NSTEPS+1) :: RANDOM
        REAL(KIND=dp), DIMENSION(2,NSTAT) :: TDECOE, TDECOEO
        INTEGER :: NSTEPS
      COMMON /MDSIM4/ NSTEPS


      CHARACTER*60 :: CARD
      INTEGER :: IWR,IST, I,J,K, ISP, LEC
      REAL(KIND=dp) :: DTAU, DT
      REAL(KIND=dp), DIMENSION(NUM,NUM) :: VECA, VECB, VECOLD, VECBLD
      REAL(KIND=dp), DIMENSION(NUM) :: MOE
      REAL(KIND=dp), DIMENSION(:,:), ALLOCATABLE :: VRO, VROLD
      REAL(KIND=dp), DIMENSION(NSTAT) :: EE, EOLD
      REAL(KIND=dp), DIMENSION(3,MXATM) :: EGOLD
      REAL(KIND=dp), DIMENSION(NSTAT*NSTAT) :: STAS, NACT
C
      ALLOCATE(VRO(NA*(NQMT-NB),NSTAT), VROLD(NA*(NQMT-NB),NSTAT))
      IWR=35
      IF (MASWRK) THEN
         READ(IWR,IOSTAT=IST)
         CALL SEQOPN(IWR,'QOPWRF','UNKNOWN',.FALSE.,'FORMATTED')
      ENDIF
C     DTAU  = DT/AU2SEC              !! DT expectred to be TIME DIFFERENCE or STEP !!!
      DTAU = 20.67D+00
      ISP = 0
      IST = 0
      LEC = NA*(NQMT-NB)-1                        ! for MRSFS
      CALL FNDGRP(IWR,' $UMSTEP',IST)
C
      IF ( IST .NE. 1) THEN
         CALL SEQREW(IWR)
         READ(IWR,'(10X,I8)') ISP
         READ(IWR,'(8X,F8.4)') DTAU
         ISP = ISP+1
         CALL FNDGRP(IWR,' $DATA  ',IST)
         DO I=1,NAT
            READ (IWR,'(16X,3F20.10)') ( COLD(J,I), J=1,3 )
         ENDDO
         CALL FNDGRP(IWR,' $VEC   ',IST)
         DO I=1,NUM
            READ(IWR,'(5x,5E15.8)')(VECOLD(J,I),J=1,NUM)
         ENDDO
         CALL DAWRIT(IDAF,IODA,VECOLD,NUM*NUM,700,0)        !  ALPHA MO
         CALL DAWRIT(IDAF,IODA,VECOLD,NUM*NUM,701,0)        !  BETA = ALPHA MO !!!
         CALL FNDGRP(IWR,' $XVEC  ',IST)
         DO I=1,NSTAT
            READ(IWR,'(22X,F20.10)') EOLD(I)
            READ(IWR,'(5(1x,F14.8))')(VROLD(J,I),J=1,LEC)
         ENDDO
         CALL DAWRIT(IDAF,IODA,VROLD,LEC*NSTAT,704,0)   ! old XVEC from previous step
      ELSE IF ( IST .EQ. 1) THEN
         CALL SEQREW(IWR)
      ENDIF
C
      WRITE (IWR,'(A,I8)') " $UMSTEP =", ISP
      WRITE (IWR,'(A,F8.4,/)') " $DTAU =", DTAU
      IF (MASWRK) WRITE (IW,'(/,A,I8)') "EXTERNAL MD STEP =", ISP
C
      WRITE (IWR,'(A)') " $DATA"
      DO I=1,NAT
         WRITE (IWR,'(1X,A8,A2,F5.1,3F20.10)')
     *      ANAM(I), BNAM(I), ZAN(I), ( CR(J,I)*UNITS, J=1,3 )
      ENDDO
      WRITE (IWR,'(A,/)') " $END"
C
      CALL DAREAD(IDAF,IODA,VECA,NUM*NUM,15,0)    !  new ALPHA MO -> old 700
      CALL DAREAD(IDAF,IODA,MOE,NUM,17,0)
      CALL DCOPY(NUM*NUM,VECA,1,VECB,1)           !  new BETA=ALPHA MO -> 701
      WRITE (IWR,'(A)') " $VEC"
      CALL PUSQLF(IWR,VECA,NUM,NUM,NUM,0)
      WRITE (IWR,'(5(1x,F14.8))')(MOE(J),J=1,NUM)
      WRITE (IWR,'(A,/)') " $END"
C
      CALL DAREAD(IDAF,IODA,VRO,LEC*NSTAT,705,0)
      CALL DAREAD(IDAF,IODA,EE,NSTAT   ,706,0)
      write (IWR,'(A)') " $XVEC"
      DO I=1,NSTAT
         WRITE (IWR,'(A,I4,A,F20.10)')
     *     " STATE #", I, "  ENERGY =", ESCF+EE(I)
         WRITE (IWR,'(5(1x,F14.8))')(VRO(J,I),J=1,LEC)
      ENDDO
      WRITE (IWR,'(A,/)') " $END"
C
      WRITE (IWR,'(A)') " $GRAD  "
      DO I=1,NAT
         WRITE (IWR,'(1X,A8,3F16.8)') ANAM(I), ( EG(J,I), J=1,3 )
      ENDDO
      WRITE (IWR,'(A,/)') " $END"
      WRITE (IWR,'(A,2F14.8)') " TOTAL ENERGY =",ESCF
C
      IF ( IST .NE. 1) THEN
         CALL DSCAL(3*NAT,1/UNITS,COLD,1)
C        COMPUTE NUMERICAL NACVS BY FINITE DIFFERENCE METHOD OVERLAP
         I = NUM*NUM
         J = (I+NUM)/2
         K = NA*(NQMT-NB)
         CALL MRSFOV(STAS,VROLD,VRO,NQMT,NA,NB,NUM,J,I,NSTAT,K)
         CALL NACVFD(NACT,STAS,DTAU,NSTAT)
         WRITE (IW,'(/,5X,A)') "---------------------------------"
         WRITE (IW,'(5X,A)')   "NONADIABATIC COUPLING TERM (A.U.)"
         WRITE (IW,'(5X,A)')   "BY USING FINITE DIFFERENCE APPROX"
         WRITE (IW,'(9X,A)')   "(<PHI^{I}|D/DT|PHI^{J}>)"
         WRITE (IW,'(5X,A,/)') "---------------------------------"
         WRITE (IWR,'(/,A)') "$NACT"
         K = 0
         DO J = 1,NSTAT
            K = K+1
            WRITE(IW,'(I5,10F11.6)') J,(NACT(K+I-1),I=1,
     *        NSTAT*NSTAT,NSTAT)
            WRITE(IWR,'(I5,10F11.6)') J,(NACT(K+I-1),I=1,
     *        NSTAT*NSTAT,NSTAT)
         ENDDO
         WRITE (IWR,'(/,A)') "$END"
         WRITE (IW,*) " "
      ENDIF
C
      CALL SEQCLO(IWR,'keep')
C
      RETURN
      END
