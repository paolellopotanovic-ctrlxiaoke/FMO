C     The AAMBS may be loaded by atom:
C
C     
C     SUBROUTINE AAMBS_ATOM(E,C,NUCZ,MCP,IVV)
C      - Stores the atomic AAMBS to an exponent array and a coefficient 
C        array
C
C     Alternatively, individual AAMBS orbitals may be loaded:
C
C     INTEGER FUNCTION AAMBS_CONTRACTION_LENGTH(NUCZ,N,L,IVVTYP)
C      - Returns the contraction length for a particular orbital
C
C     SUBROUTINE AAMBS_EXPONENT(NUCZ,N,L,L1,E,INDEX,IVVTYP)
C      - Appends the exponents for a particular orbital to the end
C        of the exponent array
C
C     SUBROUTINE AAMBS_COEFFICIENT(NUCZ,N,L,L1,C,INDEX,IVVTYP)
C     - Appends the coefficients for a particular orbital to the end
C       of the coefficient array
C
C
C     The arguments have the following definitions:
C
C     NUCZ   = Integer value of the true nuclear charge
C     N      = Principal quantum number of the orbital
C     L      = Orbital angular momentum quantum number of the orbital
C     L1     = Contraction length
C     E      = Exponent array
C     C      = Coefficient array
C     INDEX  = Counter; set to zero for the first orbital read
C     IVVTYP = Selects the non-relativistic (0) or relativistic (1)
C              version of the AAMBS and errors out if any other number
C              is chosen.
C     IVV    = Same as IVVTYP, but used locally only in AAMBS_ATOM
C     MCP    = Flag that specifies whether the MCP core should be
C              excluded from the atomic AAMBS
C
C
C     EXAMPLE #1: Load the non-relativistic AAMBS for carbon:
C
C     CALL AAMBS_ATOM(E,C,6,.FALSE.,0)
C
C
C     EXAMPLE #2: Retrieve the 2p orbital of carbon
C
C        Step 1: Get the contraction length. 
C
C                LP = AAMBS_CONTRACTION_LENGTH(6,2,1,IVV)
C
C        Step 2: Initialize a array index. This will be non-zero if 
C                appending the AAMBS for this single orbital to a 
C                larger array, e.g., building up an atomic AAMBS one
C                orbital at a time.
C
C           INDEX = 0
C
C        Step 3: Load the AAMBS exponents for the selected orbital into
C                the exponent array, E. INDEX is the starting position 
C                for this orbital in the array, E.
C
C                CALL AAMBS_EXPONENT(6,2,1,LP,E,INDEX,IVV)
C
C        Step 4: Load the AAMBS coefficients for the selected orbital
C                into the coefficient array, C.  INDEX should be the
C                same for AAMBS_EXPONENT and AAMBS_COEFFICIENT.
C
C                CALL AAMBS_COEFFICIENT(6,2,1,LP,E,INDEX,IVV)
C
C
CGS
C*MODULE AAMBS   *DECK AAMBS_ATOM
!>   @brief
!>
!>   @author   George Schoendorff
!>    - May 22, 2019
!>
!>   @param E     is the exponent array
!>   @param C     is the coefficient arry
!>   @param NUCZ  is the nuclear charge of the atom
!>   @param MCP   indicates if a model core potential is being used
!>   @param IVV   selects the AAMBS type,i.e., a local version of
!>                IVVTYP
!>
      SUBROUTINE AAMBS_ATOM(E,C,NUCZ,MCP,IVV)
      use prec, only:dp
      use comm_vvopar, only: mxgmbs
      use constants, only: atom
C-----------------------------------------------------------------------
      IMPLICIT NONE
C      
      INTEGER, INTENT(IN) :: NUCZ,IVV
      LOGICAL, INTENT(IN) :: MCP
      REAL(KIND=DP), INTENT(INOUT) :: E(MXGMBS),C(MXGMBS)
C
C-----------------------------------------------------------------------
      INTEGER I,J,LS,LP,LD,LF,AAMBS_CONTRACTION_LENGTH
      LOGICAL PBLOCK
C
      INTEGER IW,IDAF,IODA,IP,IPK,IR,IS,NAV
C
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C-----------------------------------------------------------------------
C
C     Check AAMBS version
C
      IF(IVV.LT.0.OR.IVV.GT.1) THEN
         IF(MASWRK) WRITE(IW,9999) IVV
         CALL ABRT
      ENDIF
C
C     Initialize contraction lengths
C
      LS = 0
      LP = 0
      LD = 0
      LF = 0
C
C     Get contraction lengths
C
      LS = AAMBS_CONTRACTION_LENGTH(NUCZ,1,0,IVV)
      IF(IVV.EQ.0) THEN
         IF(NUCZ.GE.5) LP = AAMBS_CONTRACTION_LENGTH(NUCZ,2,1,IVV)
      ELSE
         IF(NUCZ.GE.3) LP = AAMBS_CONTRACTION_LENGTH(NUCZ,2,1,IVV)
      ENDIF
      IF(IVV.EQ.0) THEN
         IF(NUCZ.GE.21) LD = AAMBS_CONTRACTION_LENGTH(NUCZ,3,2,IVV)
      ELSE
         IF(NUCZ.GE.19) LD = AAMBS_CONTRACTION_LENGTH(NUCZ,3,2,IVV)
      ENDIF
      IF(NUCZ.GE.57) LF = AAMBS_CONTRACTION_LENGTH(NUCZ,4,3,IVV)
C
C     Coefficients
C
C     Hydrogen - Helium
C
      IF(NUCZ.LE.2) THEN
         J = 0
         CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
C
C     Lithium - Neon
C
      ELSEIF(NUCZ.GE.3.AND.NUCZ.LE.10) THEN
         J = 0
         IF(.NOT.MCP) CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
         IF((NUCZ.GE.3.AND.IVV.EQ.1).OR.(NUCZ.GE.5.AND.IVV.EQ.0))
     *      CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
C
C     Sodium - Argon
C
      ELSEIF(NUCZ.GE.11.AND.NUCZ.LE.18) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
         IF((IVV.EQ.1.AND.(NUCZ.EQ.11.OR.NUCZ.EQ.12)).OR.
     *      NUCZ.GE.13)
     *      CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
C
C     Potassium - Calcium
C
      ELSEIF(NUCZ.GE.19.AND.NUCZ.LE.20) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
         IF(IVV.EQ.1) THEN
            IF(.NOT.MCP) CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
         ENDIF
C
C     Scandium - Zinc
C
      ELSEIF(NUCZ.GE.21.AND.NUCZ.LE.30) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
         IF(IVV.EQ.1) CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
C
C     Gallium - Krypton
C
      ELSEIF(NUCZ.GE.31.AND.NUCZ.LE.36) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
C
C     Rubidium - Strontium
C
      ELSEIF(NUCZ.GE.37.AND.NUCZ.LE.38) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
         IF(.NOT.MCP) CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
         IF(IVV.EQ.1) THEN
            IF(.NOT.MCP) CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
         ENDIF
C
C     Ytrrium - Cadmium
C
      ELSEIF(NUCZ.GE.39.AND.NUCZ.LE.48) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
         IF(.NOT.MCP) CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
         IF(IVV.EQ.1) CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
C
C     Indium - Xenon
C
      ELSEIF(NUCZ.GE.49.AND.NUCZ.LE.54) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
C
CGS
C     IVV = 1 only from Cesium to the end of the periodic table
C     (It makes no sense to have a non-relativistic AAMBS for the
C     heavy elements.)
CGS
C
C     Cesium - Barium
C
      ELSEIF(NUCZ.GE.55.AND.NUCZ.LE.56) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,0,LS,C,J,IVV)
         IF(.NOT.MCP) CALL AAMBS_COEFFICIENT(NUCZ,5,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,1,LP,C,J,IVV)
C
C     Lanthanum - Ytterbium
C
      ELSEIF(NUCZ.GE.57.AND.NUCZ.LE.70) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)           
            CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,4,3,LF,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,1,LP,C,J,IVV)
C
C     Lutetium - Mercury
C
      ELSEIF(NUCZ.GE.71.AND.NUCZ.LE.80) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,3,LF,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,5,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,1,LP,C,J,IVV)
C
C     Thallium - Radon
C
      ELSEIF(NUCZ.GE.81.AND.NUCZ.LE.86) THEN
         J = 0
         IF(.NOT.MCP) THEN
            CALL AAMBS_COEFFICIENT(NUCZ,1,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,2,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,3,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,2,LD,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,5,0,LS,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,5,1,LP,C,J,IVV)
            CALL AAMBS_COEFFICIENT(NUCZ,4,3,LF,C,J,IVV)
         ENDIF
         CALL AAMBS_COEFFICIENT(NUCZ,5,2,LD,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,0,LS,C,J,IVV)
         CALL AAMBS_COEFFICIENT(NUCZ,6,1,LP,C,J,IVV)
      ELSE
         IF(MASWRK) WRITE(IW,9998) ATOM(NUCZ)
         CALL ABRT
      ENDIF
C
C     Exponents
C
C     --- The PBLOCK logical checks if we need to reorder the
C         exponents from SPD to DSP for the model core potentials
C
      PBLOCK = .FALSE.
      IF(NUCZ.GE.31.AND.NUCZ.LE.36) PBLOCK=.TRUE.
      IF(NUCZ.GE.49.AND.NUCZ.LE.54) PBLOCK=.TRUE.
      IF(NUCZ.GE.81.AND.NUCZ.LE.86) PBLOCK=.TRUE.
      J = 0
      IF(MCP.AND.PBLOCK) THEN
         CALL AAMBS_EXPONENT(NUCZ,3,2,LD,E,J,IVV)
         CALL AAMBS_EXPONENT(NUCZ,1,0,LS,E,J,IVV)
         CALL AAMBS_EXPONENT(NUCZ,2,1,LP,E,J,IVV)
      ELSE
         CALL AAMBS_EXPONENT(NUCZ,1,0,LS,E,J,IVV)
         IF(IVV.NE.1) THEN
            IF(NUCZ.GE.5) CALL AAMBS_EXPONENT(NUCZ,2,1,LP,E,J,IVV)
         ELSE
            IF(NUCZ.GE.3) CALL AAMBS_EXPONENT(NUCZ,2,1,LP,E,J,IVV)
         ENDIF
         IF(IVV.NE.1) THEN
            IF(NUCZ.GE.21) CALL AAMBS_EXPONENT(NUCZ,3,2,LD,E,J,IVV)
         ELSE
            IF(NUCZ.GE.19) CALL AAMBS_EXPONENT(NUCZ,3,2,LD,E,J,IVV)
         ENDIF
         IF(NUCZ.GE.57) CALL AAMBS_EXPONENT(NUCZ,4,3,LF,E,J,IVV)
      ENDIF
C
 9999 FORMAT(/1X,'INVALID AAMBS VERSION=',I1,' IN AAMBS_ORBITAL')
 9998 FORMAT(/1X,'INVALID ATOM IN AAMBS_ORBITAL: ',A2)
C
      RETURN
      END
C
C*MODULE AAMBS   *DECK AAMBS_CONTRACTION_LENGTH
!>
!>    @brief   Retrieves the contraction length of an AAMBS orbital 
!>
!>    @author  George Schoendorff
!>     - August 5, 2022
!>
!>    @param NUCZ     is the true nuclear charge of the atom
!>    @param N        is the principal quantum number
!>    @param L        is the orbital angular momentum quantum number
!>    @parem IVVTYP   selects the relativistic or non-relativisitc AAMBS 
!>
      INTEGER FUNCTION AAMBS_CONTRACTION_LENGTH(NUCZ,N,L,IVVTYP) 
     *                 RESULT(LC)
      use mx_limits, only: mxatm
      use aambs_limits
      use constants, only: atom
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: NUCZ,N,L,IVVTYP
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      INTEGER IW,IDAF,IODA,IP,IPK,IR,IS,NAV
C
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C-----------------------------------------------------------------------
C
      IF(IVVTYP.EQ.0) THEN
         SELECT CASE (NUCZ)
            CASE(1)
               IF(N.EQ.1.AND.L.EQ.0) THEN
                  LC = NEHS
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(2)
               IF(N.EQ.1.AND.L.EQ.0) THEN
                  LC = NEHES
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(3)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NELIS
C               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
C                  LC = NELIP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(4)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEBES
C               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
C                  LC = NEBEP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(5)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEBS
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEBP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(6)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NECS
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NECP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(7)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NENS
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NENP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(8)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEOS
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEOP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(9)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEFS
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEFP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(10)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NENES
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NENEP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(11)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NENAS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NENAP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(12)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEMGS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEMGP
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(13)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEALS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEALP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(14)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NESIS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NESIP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(15)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEPS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEPP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(16)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NESS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NESP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(17)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NECLS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NECLP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(18)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEARS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEARP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(19)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEKS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEKP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = NEKD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(20)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECAS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECAP
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  LC = NECAD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(21)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NESCS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NESCP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NESCD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(22)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NETIS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NETIP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NETID
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(23)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEVS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEVP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEVD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(24)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECRS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECRP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECRD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(25)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEMNS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEMNP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEMND
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(26)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEFES
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEFEP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEFED
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(27)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECOS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECOP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECOD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(28)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NENIS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NENIP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NENID
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(29)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECUS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECUP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECUD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(30)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEZNS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEZNP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEZND
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(31)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEGAS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEGAP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEGAD
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(32)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEGES
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEGEP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEGED
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(33)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEASS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEASP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEASD
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(34)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NESES
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NESEP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NESED
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(35)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEBRS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEBRP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEBRD
C               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
C                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(36)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEKRS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEKRP
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEKRD
C               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
C                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(37)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NERBS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NERBP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NERBD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(38)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NESRS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NESRP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NESRD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(39)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEYS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEYP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEYD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(40)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEZRS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEZRP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEZRD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(41)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NENBS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NENBP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NENBD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(42)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEMOS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEMOP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEMOD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(43)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NETCS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NETCP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NETCD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(44)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NERUS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NERUP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NERUD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(45)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NERHS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NERHP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NERHD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(46)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEPDS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEPDP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEPDD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(47)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEAGS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEAGP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEAGD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(48)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NECDS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NECDP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NECDD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(49)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEINS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEINP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEIND
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(50)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NESNS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NESNP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NESND
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(51)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NESBS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NESBP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NESBD
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(52)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NETES
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NETEP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NETED
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(53)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEIS
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEIP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEID
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(54)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEXES
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEXEP
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEXED
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE DEFAULT
               IF(MASWRK) WRITE(IW,9998)
               CALL ABRT
         END SELECT
      ELSEIF(IVVTYP.EQ.1) THEN
         SELECT CASE (NUCZ)
            CASE(1)
               IF(N.EQ.1.AND.L.EQ.0) THEN
                  LC = 8
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(2)
               IF(N.EQ.1.AND.L.EQ.0) THEN
                  LC = 8
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(3)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NELISR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NELIPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(4)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEBESR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEBEPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(5)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEBSR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEBPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(6)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NECSR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NECPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(7)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NENSR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NENPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(8)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEOSR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEOPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(9)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NEFSR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NEFPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(10)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  LC = NENESR
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  LC = NENEPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(11)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NENASR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NENAPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(12)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEMGSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEMGPR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT      
               ENDIF
            CASE(13)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN   
                  LC = NEALSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEALPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT            
               ENDIF
            CASE(14)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NESISR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NESIPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(15)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEPSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEPPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(16)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NESSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NESPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(17)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NECLSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NECLPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(18)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  LC = NEARSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  LC = NEARPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = 12
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(19)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEKSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEKPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEKDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(20)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECASR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECAPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECADR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(21)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NESCSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NESCPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NESCDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(22)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NETISR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NETIPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NETIDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(23)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEVSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEVPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEVDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(24)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECRPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECRDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(25)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEMNSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEMNPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEMNDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(26)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEFESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEFEPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEFEDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(27)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECOSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECOPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECODR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(28)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NENISR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NENIPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NENIDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(29)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NECUSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NECUPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NECUDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(30)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEZNSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEZNPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEZNDR
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(31)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEGASR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEGAPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEGADR
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(32)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEGESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEGEPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEGEDR
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(33)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEASSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEASPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEASDR
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(34)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NESESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NESEPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NESEDR
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(35)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEBRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEBRPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEBRDR
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(36)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  LC = NEKRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  LC = NEKRPR
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  LC = NEKRDR
               ELSEIF(L.EQ.2.AND.N.EQ.4) THEN
                  LC = 16
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(37)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NERBSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NERBPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NERBDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(38)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NESRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NESRPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NESRDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(39)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEYSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEYPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEYDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(40)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEZRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEZRPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEZRDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(41)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NENBSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NENBPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NENBDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(42)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEMOSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEMOPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEMODR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(43)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NETCSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NETCPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NETCDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(44)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NERUSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NERUPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NERUDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(45)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NERHSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NERHPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NERHDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(46)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEPDSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEPDPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEPDDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(47)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEAGSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEAGPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEAGDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(48)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NECDSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NECDPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NECDDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(49)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEINSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEINPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEINDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(50)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NESNSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NESNPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NESNDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(51)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NESBSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NESBPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NESBDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(52)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NETESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NETEPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NETEDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(53)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEISR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEIPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEIDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(54)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  LC = NEXESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  LC = NEXEPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  LC = NEXEDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(55)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN
                  LC = NECSSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NECSPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NECSDR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(56)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN
                  LC = NEBASR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEBAPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEBADR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(57)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN          
                  LC = NELASR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NELAPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NELADR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NELAFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(58)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NECESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NECEPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NECEDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NECEFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(59)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEPRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEPRPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEPRDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEPRFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(60)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NENDSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NENDPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NENDDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NENDFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(61)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEPMSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEPMPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEPMDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEPMFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(62)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NESMSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NESMPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NESMDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NESMFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(63)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEEUSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEEUPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEEUDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEEUFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(64)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEGDSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEGDPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEGDDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEGDFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(65)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NETBSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NETBPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NETBDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NETBFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(66)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEDYSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEDYPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEDYDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEDYFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(67)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEHOSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEHOPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEHODR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEHOFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(68)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEERSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEERPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEERDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEERFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(69)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NETMSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NETMPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NETMDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NETMFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(70)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEYBSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEYBPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEYBDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEYBFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(71)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NELUSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NELUPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NELUDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NELUFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(72)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEHFSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEHFPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEHFDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEHFFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(73)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NETASR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NETAPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NETADR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NETAFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(74)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEWSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEWPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEWDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEWFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(75)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NERESR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEREPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEREDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEREFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(76)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEOSSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEOSPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEOSDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEOSFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(77)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEIRSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEIRPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEIRDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEIRFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(78)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEPTSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEPTPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEPTDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEPTFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(79)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEAUSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEAUPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEAUDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEAUFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(80)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEHGSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEHGPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEHGDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEHGFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(81)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NETLSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NETLPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NETLDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NETLFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(82)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEPBSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEPBPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEPBDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEPBFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(83)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEBISR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEBIPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEBIDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEBIFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(84)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEPOSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEPOPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEPODR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEPOFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(85)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NEATSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NEATPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NEATDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NEATFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(86)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  LC = NERNSR
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  LC = NERNPR
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  LC = NERNDR
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  LC = NERNFR
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE DEFAULT
               IF(MASWRK) WRITE(IW,9998)
               CALL ABRT
         END SELECT
      ELSE
         IF(MASWRK) WRITE(IW,9996) IVVTYP
      ENDIF
C
 9999 FORMAT(/1X,'*** ERROR: INVALID QUANTUM NUMBER FOR ELEMENT ',A2,
     *       ' N=',I1,' L=',I1)
 9998 FORMAT(/1X,'*** ERROR: THE AAMBS DOES NOT YET SUPPORT FR-OG')
 9997 FORMAT(/1X,'*** ERROR: ONLY IVVTYP=1 IS SUPPORTED FOR CS-RN')
 9996 FORMAT(/1X,'*** ERROR: INVALID CHOICE FOR IVVTYP=',I2)
C
      RETURN
      END
C
C*MODULE AAMBS   *DECK AAMBS_EXPONENT
!>
!>    @brief   Retrieves the contraction length of an AAMBS orbital 
!>
!>    @author  George Schoendorff
!>     - August 7, 2022
!>
!>    @param NUCZ     is the true nuclear charge of the atom
!>    @param N        is the principal quantum number
!>    @param L        is the orbital angular momentum quantum number
!>    @param L1       is the length of the coefficient
!>    @param E        is the exponent array on exit
!>    @param INDEX    is the zero point of the indexing array
!>    @param IVVTYP   selects the AAMBS type
!>
      SUBROUTINE AAMBS_EXPONENT(NUCZ,N,L,L1,E,INDEX,IVVTYP)
      use constants, only: atom
      use aambsrel
      use aambsnorel
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: NUCZ,N,L,L1,IVVTYP
      INTEGER, INTENT(INOUT) :: INDEX
      DOUBLE PRECISION, INTENT(OUT) :: E(*)
C
      INTEGER I,NGAU,AAMBS_CONTRACTION_LENGTH,LEXP
      DOUBLE PRECISION :: A,B,D,G
C
C-----------------------------------------------------------------------
C
      INTEGER IW,IDAF,IODA,IP,IPK,IR,IS,NAV
C
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C-----------------------------------------------------------------------
CGS
C     The non-relativistic coefficients were generated with exponents
C     that have different precision than those provided in
C     MODULE AAMBSNOREL (modf77_aambs.src). In a few cases this results
C     deviations from 1.0 for the normalization that are outside
C     the acceptable tolerance of 1.0D-6. In most cases recomputing
C     the coefficients with the fixed precision exponents solves this
C     problem, though the 3p orbital for ruthenium remains just outside
C     the tolerance. So in this case we simply loosen the tolerance
C     for ruthenium only in ETGTO.
CGS
C
      IF(IVVTYP.EQ.0) THEN
C         IF(MASWRK) WRITE(IW,9997)
C         CALL ABRT
         SELECT CASE (NUCZ)
            CASE(1)
               IF(L.EQ.0.AND.N.EQ.1) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHS(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(2)
               IF(L.EQ.0.AND.N.EQ.1) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHES(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(3)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELIS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELIP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(4)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBEP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(5)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(6)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(7)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(8)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(9)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(10)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENEP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(11)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENAS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENAP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(12)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMGS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMGP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(13)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EALS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EALP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(14)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESIS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESIP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(15)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(16)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(17)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECLS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECLP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(18)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EARS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EARP(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(19)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKP(I)
                  ENDDO
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  DO I = 1,L1
C                     INDEX = INDEX + 1
C                     E(INDEX) = EKD(I)
C                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(20)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECAS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECAP(I)
                  ENDDO
C               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
C                  DO I = 1,L1
C                     INDEX = INDEX + 1
C                     E(INDEX) = ECAD(I)
C                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(21)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESCS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESCP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESCD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(22)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETIS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETIP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETID(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(23)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EVS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EVP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EVD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(24)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECRS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECRP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECRD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(25)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMNS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMNP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMND(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(26)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFEP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFED(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(27)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECOS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECOP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECOD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(28)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENIS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENIP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENID(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(29)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECUS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECUP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECUD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(30)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZNS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZNP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZND(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(31)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGAS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGAP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGAD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(32)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGEP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGED(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(33)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EASS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EASP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EASD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(34)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESEP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESED(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(35)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBRS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBRP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBRD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(36)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKRS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKRP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKRD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(37)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERBS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERBP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERBD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(38)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESRS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESRP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESRD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(39)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(40)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZRS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZRP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZRD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(41)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENBS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENBP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENBD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(42)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMOS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMOP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMOD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(43)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETCS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETCP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETCD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(44)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERUS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERUP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERUD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(45)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERHS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERHP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERHD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(46)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPDS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPDP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPDD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(47)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAGS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAGP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAGD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(48)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECDS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECDP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECDD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(49)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EINS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EINP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIND(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(50)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESNS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESNP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESND(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(51)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESBS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESBP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESBD(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(52)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETEP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETED(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(53)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIS(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EID(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(54)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EXES(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EXEP(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EXED(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE DEFAULT 
               IF(MASWRK) WRITE(IW,9998)
               CALL ABRT
         END SELECT
      ELSEIF(IVVTYP.EQ.1) THEN
         SELECT CASE (NUCZ)
            CASE(1)
               IF(L.EQ.0.AND.N.EQ.1) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHSR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(2)
               IF(L.EQ.0.AND.N.EQ.1) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHESR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(3)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELISR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELIPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(4)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBEPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(5)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(6)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(7)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(8)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(9)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(10)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.2)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.N.EQ.2) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENEPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT       
               ENDIF
            CASE(11)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENASR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENAPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(12)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMGSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMGPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT      
               ENDIF
            CASE(13)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN   
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EALSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EALPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT            
               ENDIF
            CASE(14)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESISR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESIPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(15)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(16)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(17)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECLSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECLPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(18)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EARSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.3)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EARPR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(19)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(20)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECASR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECAPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECADR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(21)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESCSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESCPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESCDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(22)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETISR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETIPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETIDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(23)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EVSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EVPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EVDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(24)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECRDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(25)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMNSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMNPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMNDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(26)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFEPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EFEDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(27)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECOSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECOPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECODR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(28)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENISR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENIPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENIDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(29)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECUSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECUPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECUDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(30)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZNSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZNPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZNDR(I)
                  ENDDO
               ELSE        
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(31)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGASR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGAPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGADR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(32)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGEPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGEDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(33)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EASSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EASPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EASDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(34)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESEPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESEDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(35)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBRDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(36)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.N.EQ.3) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EKRDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L    
                  CALL ABRT
               ENDIF
            CASE(37)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERBPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERBDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(38)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESRDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(39)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(40)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EZRDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(41)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENBPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENBDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(42)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMOSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMOPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EMODR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(43)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETCSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETCPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETCDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(44)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERUSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERUPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERUDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(45)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERHSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERHPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERHDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(46)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPDSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPDPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPDDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(47)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAGSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAGPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAGDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(48)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECDSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECDPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECDDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(49)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EINSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EINPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EINDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(50)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESNSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESNPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESNDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(51)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESBPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESBDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(52)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETEPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETEDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(53)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EISR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(54)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EXESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EXEPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.4)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EXEDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(55)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECSSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECSPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECSDR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(56)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBASR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBAPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBADR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(57)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN          
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELASR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELAPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELADR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELAFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF
            CASE(58)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECEPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECEDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ECEFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(59)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPRDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPRFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(60)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENDSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENDPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENDDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ENDFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(61)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPMSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPMPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPMDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPMFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(62)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESMSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESMPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESMDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ESMFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(63)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EEUSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EEUPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EEUDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EEUFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(64)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGDSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGDPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGDDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EGDFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(65)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETBPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETBDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETBFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(66)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EDYSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EDYPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EDYDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EDYFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(67)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHOSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHOPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHODR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHOFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(68)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EERSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EERPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EERDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EERFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(69)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETMSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETMPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETMDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETMFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(70)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYBPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYBDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EYBFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(71)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELUSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELUPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELUDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ELUFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(72)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHFSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHFPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHFDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHFFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(73)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETASR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETAPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETADR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETAFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(74)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EWSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EWPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EWDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EWFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(75)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERESR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EREPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EREDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EREFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(76)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOSSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOSPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOSDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EOSFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(77)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIRSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIRPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIRDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EIRFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(78)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPTSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPTPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPTDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPTFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(79)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAUSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAUPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAUDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EAUFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(80)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHGSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHGPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHGDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EHGFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(81)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETLSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETLPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETLDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ETLFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(82)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPBSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPBPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPBDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPBFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(83)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBISR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBIPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBIDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EBIFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(84)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPOSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPOPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPODR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EPOFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(85)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EATSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EATPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EATDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = EATFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE(86)
               IF(L.EQ.0.AND.(N.GE.1.AND.N.LE.6)) THEN                     
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERNSR(I)
                  ENDDO
               ELSEIF(L.EQ.1.AND.(N.GE.2.AND.N.LE.6)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERNPR(I)
                  ENDDO
               ELSEIF(L.EQ.2.AND.(N.GE.3.AND.N.LE.5)) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERNDR(I)
                  ENDDO
               ELSEIF(L.EQ.3.AND.N.EQ.4) THEN
                  DO I = 1,L1
                     INDEX = INDEX + 1
                     E(INDEX) = ERNFR(I)
                  ENDDO
               ELSE
                  IF(MASWRK) WRITE(IW,9999) ATOM(NUCZ),N,L
                  CALL ABRT
               ENDIF  
            CASE DEFAULT
               IF(MASWRK) WRITE(IW,9997)
               CALL ABRT
         END SELECT
      ELSE
         IF(MASWRK) WRITE(IW,9996) IVVTYP
      ENDIF
C
 9999 FORMAT(/1X,'*** ERROR: INVALID QUANTUM NUMBER FOR ELEMENT ',A2,
     *       ' N=',I1,' L=',I1)
 9998 FORMAT(/1X,'*** ERROR: IVVTYP=0 IS VALID ONLY FOR H-XE')
 9997 FORMAT(/1X,'*** ERROR: THE AAMBS DOES NOT YET SUPPORT FR-OG')
 9996 FORMAT(/1X,'*** ERROR: INVALID CHOICE FOR IVVTYP=',I2)
C
      RETURN
      END
C
C*MODULE AAMBS   *DECK AAMBS_COEFFICIENT
!>
!>    @brief   Retrieves the coefficient array for a specified AAMBS
!>             orbital; currently works only for the relativistic
!>             AAMBS
!>
!>    @author  George Schoendorff
!>     - August 5, 2022
!>
!>    @param NUCZ     is the true nuclear charge of the atom
!>    @param N        is the principal quantum number
!>    @param L        is the orbital angular momentum quantum number
!>    @param L1       is the length of the coefficient
!>    @param C        is the coefficient array on exit
!>    @param INDEX    is the zero point of the indexing array
!>    @param IVVTYP   selects the relativistic AAMBS
!>
      SUBROUTINE AAMBS_COEFFICIENT(NUCZ,N,L,L1,C,INDEX,IVVTYP)
      use constants, only: atom
      use aambsrel
      use aambsnorel
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C
      INTEGER, INTENT(IN) :: NUCZ,N,L,L1,IVVTYP
      INTEGER, INTENT(INOUT) :: INDEX
      DOUBLE PRECISION, INTENT(OUT) :: C(*)
C
      INTEGER I
C
C-----------------------------------------------------------------------
C
      INTEGER IW,IDAF,IODA,IP,IPK,IR,IS,NAV
C
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C-----------------------------------------------------------------------
C
C     Make sure that IVVTYP is reasonable
C
      IF(IVVTYP.LT.0.OR.IVVTYP.GT.1) THEN
         IF(MASWRK) WRITE(IW,9997) IVVTYP
         CALL ABRT
      ENDIF
C
      IF(IVVTYP.EQ.0) THEN
      SELECT CASE (NUCZ)
         CASE(1)
            IF(L.EQ.0.AND.N.EQ.1) THEN
               DO I = 1,L1
                  INDEX = INDEX + 1
                  C(INDEX) = CH1S(I)
               ENDDO
            ELSE
               IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
               CALL ABRT
            ENDIF
         CASE(2)
            IF(L.EQ.0.AND.N.EQ.1) THEN
               DO I = 1,L1
                  INDEX = INDEX + 1
                  C(INDEX) = CHE1S(I)
               ENDDO
            ELSE
               IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
               CALL ABRT
            ENDIF
         CASE(3)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLI1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLI2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
C               CASE(1)
C                  SELECT CASE (N)
C                     CASE(2)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CLI2P(I)
C                        ENDDO
C                     CASE DEFAULT
C                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
C                        CALL ABRT
C                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(4)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBE2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
C               CASE(1)
C                  SELECT CASE (N)
C                     CASE(2)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CBE2P(I)
C                        ENDDO
C                     CASE DEFAULT
C                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
C                        CALL ABRT
C                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(5)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CB1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CB2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CB2P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(6)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CC1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CC2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CC2P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(7)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CN1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CN2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CN2P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(8)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CO1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CO2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CO2P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(9)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CF1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CF2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CF2P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(10)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNE2S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNE2P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(11)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA2P(I)
                        ENDDO
C                     CASE(3)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CNA3P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(12)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG2P(I)
                        ENDDO
C                     CASE(3)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CMG3P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(13)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL3P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(14)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI3P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(15)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP3P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(16)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS3P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(17)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL3P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(18)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR3S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR3P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(19)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CK4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
C               CASE(2)
C                  SELECT CASE (N)
C                     CASE(3)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CK3D(I)
C                        ENDDO
C                     CASE DEFAULT
C                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
C                        CALL ABRT
C                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(20)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CCA4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
C               CASE(2)
C                  SELECT CASE (N)
C                     CASE(3)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CCA3D(I)
C                        ENDDO
C                     CASE DEFAULT
C                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
C                        CALL ABRT
C                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(21)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CSC4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(22)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CTI4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(23)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CV4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(24)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CCR4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(25)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CMN4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(26)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CFE4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(27)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CCO4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(28)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CNI4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(29)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CCU4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(30)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN3P(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CZN4P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(31)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA4P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(32)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE4P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(33)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS4P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(34)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE4P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(35)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR4P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(36)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR4S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR4P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR3D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(37)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CRB5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB3D(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CRB4D(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(38)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CSR5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR3D(I)
                        ENDDO
C                     CASE(4)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CSR4D(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(39)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CY5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(40)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CZR5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(41)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CNB5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(42)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CMO5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(43)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CTC5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(44)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CRU5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(45)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CRH5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(46)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CPD5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(47)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CAG5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(48)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD4P(I)
                        ENDDO
C                     CASE(5)
C                        DO I = 1,L1
C                           INDEX = INDEX + 1
C                           C(INDEX) = CCD5P(I)
C                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(49)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN4P(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN5P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(50)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN4P(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN5P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(51)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB4P(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB5P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(52)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE4P(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE5P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(53)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI4P(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI5P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(54)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE1S(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE2S(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE3S(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE4S(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE5S(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE2P(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE3P(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE4P(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE5P(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE3D(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE4D(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE DEFAULT
            IF(MASWRK) WRITE(IW,9994)
            CALL ABRT
      END SELECT
C
C     End non-relativistic (IVV=0)
C
      ELSEIF(IVVTYP.EQ.1) THEN
C
      SELECT CASE (NUCZ)
         CASE(1)
            IF(L.EQ.0.AND.N.EQ.1) THEN
               DO I = 1,L1
                  INDEX = INDEX + 1
                  C(INDEX) = CH1SR(I)
               ENDDO
            ELSE
               IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
               CALL ABRT
            ENDIF
         CASE(2)
            IF(L.EQ.0.AND.N.EQ.1) THEN
               DO I = 1,L1
                  INDEX = INDEX + 1
                  C(INDEX) = CHE1SR(I)
               ENDDO
            ELSE
               IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
               CALL ABRT
            ENDIF
         CASE(3)
            SELECT CASE (L)
               CASE(0)    
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLI1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLI2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1) 
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLI2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(4)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBE2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CBE2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(5)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CB2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CB2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(6)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CC1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CC2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CC2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(7)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CN1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CN2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CN2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(8)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CO1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CO2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CO2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(9)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CF1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CF2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CF2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(10)
            SELECT CASE (L)
               CASE(0)         
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNE2SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1        
                           C(INDEX) = CNE2PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(11)
            SELECT CASE (L)
               CASE(0)     
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNA3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(12)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMG3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(13)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAL3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(14)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSI3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(15)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CP3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(16)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CS3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(17)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCL3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(18)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR3SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2) 
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAR3PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(19)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CK3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(20)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCA3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(21)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSC3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(22)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTI3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(23)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CV3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(24)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCR3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(25)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMN3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(26)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CFE3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(27)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCO3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(28)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNI3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(29)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCU3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(30)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZN3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(31)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGA3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(32)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGE3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(33)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAS3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(34)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSE3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(35)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBR3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(36)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR4SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR4PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CKR3DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(37)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRB4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(38)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CSR5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSR4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(39)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CY5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CY4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(40)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CZR5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CZR4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(41)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CNB5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CNB4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(42)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CMO5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CMO4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(43)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CTC5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTC4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(44)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CRU5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRU4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(45)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CRH5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRH4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(46)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CPD5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPD4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(47)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CAG5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAG4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(48)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CCD5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCD4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(49)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CIN5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIN4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(50)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CSN5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSN4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(51)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CSB5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSB4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(52)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CTE5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTE4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(53)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CI5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CI4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(54)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1 
                           INDEX = INDEX + 1
                           C(INDEX) = CXE5SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)  
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE5PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)  
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CXE4DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT 
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(55)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCS5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(56)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBA5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                     IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                     CALL ABRT
            END SELECT
         CASE(57)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLA4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(58)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CCE4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(59)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPR4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(60)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CND4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(61)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPM4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(62)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CSM4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(63)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CEU4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(64)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CGD4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(65)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTB4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(66)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CDY4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(67)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHO4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(68)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CER4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(69)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTM4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(70)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CYB4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(71)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CLU4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(72)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHF4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(73)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTA4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(74)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CW4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(75)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRE4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(76)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = COS4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(77)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CIR4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(78)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPT4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(79)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAU4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(80)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CHG4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(81)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CTL4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(82)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPB4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(83)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CBI4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(84)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CPO4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(85)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CAT4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
            END SELECT
         CASE(86)
            SELECT CASE (L)
               CASE(0)
                  SELECT CASE (N)
                     CASE(1)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN1SR(I)
                        ENDDO
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN2SR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN3SR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN4SR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN5SR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN6SR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(1)
                  SELECT CASE (N)
                     CASE(2)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN2PR(I)
                        ENDDO
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN3PR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN4PR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN5PR(I)
                        ENDDO
                     CASE(6)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN6PR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(2)
                  SELECT CASE (N)
                     CASE(3)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN3DR(I)
                        ENDDO
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN4DR(I)
                        ENDDO
                     CASE(5)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN5DR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE(3)
                  SELECT CASE (N)
                     CASE(4)
                        DO I = 1,L1
                           INDEX = INDEX + 1
                           C(INDEX) = CRN4FR(I)
                        ENDDO
                     CASE DEFAULT
                        IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                        CALL ABRT
                  END SELECT
               CASE DEFAULT
                  IF(MASWRK) WRITE(IW,9996) ATOM(NUCZ),N,L
                  CALL ABRT
            END SELECT
         CASE DEFAULT
            IF(MASWRK) WRITE(IW,9995)
            CALL ABRT
      END SELECT
C
C     End relativistic (IVV=1)
C
      ENDIF
C
 9997 FORMAT(/1X,'*** ERROR: INVALID VALUE OF IVVTYP =',I2)
 9996 FORMAT(/1X,'*** ERROR: INVALID QUANTUM NUMBER FOR ELEMENT ',A2,
     *       ' N=',I1,' L=',I1)
 9995 FORMAT(/1X,'*** ERROR: THE AAMBS DOES NOT YET SUPPORT FR-OG')
 9994 FORMAT(/1X,'*** ERROR: THE NON-RELATIVISITC AAMBS DOES NOT ',
     *       'ELEMENTS BEYOND Z=54 (XENON)')
C
      RETURN
      END
