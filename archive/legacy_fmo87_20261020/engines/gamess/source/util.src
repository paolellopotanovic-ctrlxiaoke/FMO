C 27 APR 21 - CF - MOVE WALL_TIME HERE
C*MODULE UTIL
C-----------------------------------------------------------------------
C     CODE THAT LIVES HERE WAS ONCE TOGETHER WITH THE NASTIEST,
C     TOTALLY UNPORTABLE CODE IN GAMESS (UNPORT). SOME CODE
C     DIDN'T DESERVE TO BE CALLED NASTY SINCE IT WAS TOTALLY
C     PORTABLE. SO THIS PLACE WAS CREATED.
C-----------------------------------------------------------------------
C
C*MODULE UTIL  *DECK WALL_TIME
C>
C>    @brief   do wall clock timings
C>
C>    @details get wall clock timings
C>
C>    @author  Dmitri Fedorov
C>
      SUBROUTINE WALL_TIME(WALL)
      IMPLICIT NONE
      DOUBLE PRECISION WALL
      INTEGER IWALL,MAX_RATE
C
      CALL SYSTEM_CLOCK(COUNT=IWALL, COUNT_RATE=MAX_RATE)
      WALL = DBLE(IWALL) / MAX_RATE
      RETURN
      END
C
C*MODULE UTIL  *DECK ABRTX
C>
C>    @brief   Print error message before calling ABRT
C>
C>    @author  ???
C>
      subroutine ABRTX(str)
      use comm_PAR, only: MASWRK
      use comm_IOFILE, only: IW
      implicit none
C     Arguments
      character*(*) str
      INTEGER :: KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
C
      if(maswrk.or.iand(modio,16384).ne.0) write(iw,*) str
      call abrt
      end
C*MODULE UTIL  *DECK FLSHBF
C>
C>    @brief   Flush the buffer for logical unit LUNIT
C>
C>    @author  probably Mike Schmidt
C>
      SUBROUTINE FLSHBF(LUNIT)
      IMPLICIT NONE
C
      INTEGER :: LUNIT
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      INTEGER :: KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      INTEGER :: ME,MASTER,NPROC,IBTYP,IPTIM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C        ----- FLUSH THE BUFFER FOR LOGICAL UNIT LUNIT -----
C        THIS ROUTINE IS MOST IMPORTANT FOR UNIX SYSTEMS,
C        WHERE OUTPUT OTHERWISE STAYS IN BUFFERS FOREVER.
C
      if(iand(modio,1).ne.0) return
      IF (MASWRK) THEN
        FLUSH(LUNIT)
      END IF
      RETURN
      END
C*MODULE UTIL  *DECK IGETGRDVECLEN
C>
C>    @brief   ???
C>
C>    @author  ???
C>
      INTEGER FUNCTION IGETGRDVECLEN(MAXVEC)
      IMPLICIT NONE
      INTEGER :: MAXVEC
      IGETGRDVECLEN = MAXVEC*3+1
      RETURN
      END
C*MODULE UTIL  *DECK PARSET
C>
C>    @brief   Set CONTROL for parallel computation
C>
C>    @author  Theresa L. Windus, Mike Schmidt
C>
      SUBROUTINE PARSET
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL DIRTRF
      COMMON /TRFOPT/ CUTTRF,NWDTRF,MPTRAN,ITRFAO,NOSYMT,IPURTF,DIRTRF
C
C     ----- SET CONTROL FOR PARALLEL COMPUTATION -----
C
C     ITRFAO CONTROLS AO INTEGRAL DISK STORAGE DURING CONVENTIONAL
C     INTEGRAL TRANSFORMATIONS, 1=DUPLICATE AO LIST IF COMMUNICATION
C     SPEED IS VERY POOR, E.G. ETHERNET, 2=DISTRIBUTE AO LIST BECAUSE
C     BROADCAST OF INTEGRALS MAY BE MORE EFFICIENT THAN DISK I/O, E.G.
C     THE SWITCH IN A SP2 MACHINE.  THIS VARIABLE CAN BE OVERRIDDEN
C     LATER BY INPUT IN $MP2 OR $TRANS GROUPS.
      ITRFAO=1
      RETURN
      END
C*MODULE UTIL  *MODULE FAST_MEMORY
C>
C>    @brief this module contains static allocatable variable for
C>           allocating dynamic memory
C>
C>    @param MEMORY - dynamic memory; must NOT be used directly
C>
C>    @author  Igor S. Gerasimov
C>
C>    @date    May, 2020 Initial release
C>
      MODULE FAST_MEMORY
         IMPLICIT NONE
         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, TARGET :: MEMORY
      END MODULE FAST_MEMORY
C*MODULE UTIL  *DECK get_ptr_addr
C>
C>    @brief   converts c_ptr to integer for pointer ariphmetics
C>
C>    @details hack of Fortran language for getting integer value of c_ptr
C>
C>    @param   ptr   pointer which address should be presented as integer
C>    @param   addr  integer value of pointer
C>
C>    @author  Igor S. Gerasimov
C>
C>    @date    May, 2020 Initial release
C>
C>    @date Jun 2021 - Peng Xu and Tosaporn Sattasathuchana
C>      - the default integer-size setting for each compiler will be
C>        different. In order to have GAMESS interfaced properly with 
C>        Verachem, it should be integer(8) ==>
C>        integer(kind=c_long).
C>
      module getptraddr
      implicit none
      contains
      elemental function get_ptr_addr(ptr) result(addr)
         use iso_c_binding, only: c_ptr, C_INTPTR_T 
         implicit none
         type(c_ptr), intent(in) :: ptr
         integer(kind=C_INTPTR_T) :: addr
         addr = transfer(ptr, addr)
      end function get_ptr_addr
      end module getptraddr
C*MODULE UTIL  *DECK SETFM
C>
C>    @brief   Obtain the total memory pool
C>
C>    @author  Steve Elbert, Mike Schmidt
C>
      SUBROUTINE SETFM(IPAR)
C-----------------------------------------------------------------------
C     ----- FAST MEMORY (FM) MANAGEMENT ROUTINES -----
C     NOTE THAT THERE ARE SIX ROUTINES, NAMELY
C        SETFM AND BIGFM ARE CALLED ONLY ONCE EACH,
C             TO OBTAIN AND RELEASE THE TOTAL MEMORY POOL.  IN UNIX,
C             THESE NORMALLY WIND UP CALLING 'MALLOC' AND 'FREE'.
C        VALFM, GETFM, RETFM ARE USED TOGETHER, AS TRIPLETS OF CALLS:
C           VALFM TELLS THE LOCATION OF THE AVAILABLE MEMORY
C           GETFM ASKS TO USE SOME OF THAT MEMORY
C           RETFM GIVES THE MEMORY BACK WHEN FINISHED
C              THESE TRIPLETS CAN BE NESTED, SO LONG AS THE CALLS
C              TO RETFM ARE IN REVERSE ORDER OF THE CALLS TO GETFM,
C              THAT IS TO SAY, THE IMPLEMENTATION IS A STACK.
C        GOTFM MAY BE USED TO LEARN THE TOTAL AMOUNT OF AVAILABLE
C              MEMORY, E.G. NOT YET COMMITTED BY GETFM'S.
C
C     ALL MEMORY CALCULATIONS HERE ARE TO BE BASED ON -WORDS-,
C     WHERE A WORD IS DEFINED AS BEING A 64 BIT QUANTITY.
C
C     SETFM(IPAR) - ON ENTRY, IPAR IS MAXIMUM MEMORY DESIRED BY USER.
C                   ON EXIT, IPAR IS MAXIMUM MEMORY ACTUALLY AVAILABLE.
C                   ALLOCATES THE MEMORY POOL FROM THE SYSTEM.
C
C     BIGFM(IPAR) - ON EXIT, IPAR IS MAXIMUM MEMORY EVER USED.
C                   FREES THE ENTIRE MEMORY POOL TO THE SYSTEM.
C
C     VALFM(IPAR) - ON EXIT, IPAR IS AN OFFSET (LTOP) TO THE HIGHEST
C                   POSITION IN MEMORY CURRENTLY USED.
C
C     GETFM(IPAR) - ON ENTRY, IPAR IS THE MEMORY PIECE REQUESTED.
C                   IT IS ALLOCATED AT LTOP+1, AFTER WHICH THE
C                   LTOP POINTER IS AUTOMATICALLY ADJUSTED UPWARDS.
C
C     RETFM(IPAR) - ON ENTRY, IPAR IS THE MEMORY PIECE TO BE FREED.
C                   THE LTOP POINTER IS AUTOMATICALLY ADJUSTED DOWN.
C
C     GOTFM(IPAR) - ON EXIT, IPAR IS THE CURRENTLY UNUSED MEMORY.
C
C--------------------------------------------------------------------
C
C     A TYPICAL MEMORY REQUEST SEQUENCE IS VALFM,GETFM,...,RETFM.
C     IT IS CRUCIAL THAT RETFM'S BE IN INVERSE ORDER OF GETFM CALLS.
C     THE DYNAMIC POOL DURING THE MIDDLE ELLIPSIS OF THE SEQUENCE
C
C         CALL VALFM,GETFM,...VALFM,GETFM,...RETFM,...RETFM
C     WITH ARG=      NEED1          NEED2    NEED2    NEED1
C
C     LOOKS LIKE THIS
C
C            .....RESERVED........      AVAILABLE
C            <-------><----------><---------------------->
C            X  NEED1    NEED2   Y                       Z
C
C     WHERE X=LOFFS, Y=LTOP, Z=MEMLIM
C
C     THERE ARE TWO IMPLEMENTATION STRATEGIES BELOW.
C
C     ONE APPROACH IS TO DECLARE A LARGE FIXED LENGTH ARRAY
C     X(MEMSIZ) JUST BELOW.  THE USER REQUESTED MEMORY -MEMLIM-
C     CANNOT EXCEED THE "STATIC" SIZE -MEMSIZ- UNLESS -MEMSIZ-
C     IS INCREASED HERE IN THE SOURCE, AND THE CODE RECOMPILED.
C
C     THE OTHER ALLOCATES ONLY A SINGLE WORD X(1), AND ACTUALLY
C     ALLOCATES THE POOL ELSEWHERE.  -LOFFS- MUST BE THE DISTANCE
C     IN WORDS FROM X(1) TO WHEREVER THIS POOL IS ALLOCATED.
C     THIS STRATEGY IS TRULY "DYNAMIC", THE CALL TO SETFM CAN REQUEST
C     ANY AMOUNT UP TO WHATEVER LIMITS THE OPERATING SYSTEM IMPOSES.
C     IN THIS CASE, -MEMSIZ- IS JUST A DEFAULT, NOT AN UPPER BOUND.
C
C-----------------------------------------------------------------------
C
      USE ISO_C_BINDING, ONLY: C_PTR, C_LOC
      USE FAST_MEMORY, ONLY: MEMORY
      use getptraddr, only: get_ptr_addr
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      PARAMETER (MEMSIZ= 1 000 000)
      COMMON /FMCOM / X(1)
      DOUBLE PRECISION, TARGET :: X
      TYPE(C_PTR) :: LOCX_PTR, LOCMEM_PTR
C
C
C     ----- INITIALIZE FAST MEMORY (FM) MANAGEMENT -----
C     -LOCX- IS THE ADDRESS OF THE FIRST ELEMENT OF X.
C     FOR STATIC IMPLEMENTATIONS, THIS SHOULD BE JUST 1.
C
      LOCX_PTR = C_LOC(X)
C
C     -MEMLIM- IS THE DESIRED MAXIMUM SIZE OF THE MEMORY POOL.
C
      MEMLIM = IPAR
      IF (MEMLIM .LE. 0) MEMLIM = MEMSIZ
C
C     ----- CHECK FOR AVAILABILITY OF MEMORY -----
C     -MEMSIZ- IS AN ABSOLUTE BOUND FOR STATIC IMPLEMENTATIONS
C
      MEMAVL = MEMLIM
C
      IF (MEMLIM .GT. MEMAVL) THEN
         IF (MASWRK)
     *   WRITE(IW,FMT='('' REQUESTED AMOUNT OF MEMORY ('',I12,
     *      '' WORDS) EXCEEDS THE AMOUNT AVAILABLE ('',I12,'')'')')
     *               MEMLIM,MEMAVL
         CALL ABRT
      END IF
C
C     ----- ALLOCATE THE DYNAMIC MEMORY POOL -----
C     -LOCMEM- IS THE STARTING ADDRESS OF THE DYNAMIC POOL.
C     FOR STATIC IMPLEMENTATIONS, THIS SHOULD JUST BE 1.
C
C   The following statement is problematic. ALLOCATE needs
C   memory in words, unlike MALLOC.
C        ALLOCATE(MEMORY((MEMLIM+2)*8))
C
      ALLOCATE(MEMORY(MEMLIM),STAT=IERR)
C
c     if(iand(LOCMEM,4).ne.0) write(6,*) 'mempnt=',LOCMEM
c     check alignment
C
c     IF (.NOT.ALLOCATED(MEMORY)) THEN
      IF(IERR .NE. 0)THEN
         IF (MASWRK) WRITE(IW,*) MEMLIM,' WORDS OF MEMORY UNAVAILABLE'
         CALL ABRT
      END IF
C
C         COMPUTE THE OFFSET -LOFFS- FROM X(1) TO THE BEGINNING
C         OF THE DYNAMIC POOL (E.G. BYTE TO WORD CONVERSION).
C         -LOFFS- WILL BE ZERO FOR STATIC IMPLEMENTATIONS.
C
      LOCMEM_PTR = C_LOC(MEMORY)
      LOCX   = get_ptr_addr(LOCX_PTR)
      LOCMEM = get_ptr_addr(LOCMEM_PTR)
      LOFFS = LOCMEM - LOCX
      LOFFS = (LOFFS+7)/8 + 1
C
      LTOP  = 0
      LENHI = 0
      MEMOK = 99
C
C     ----- REPORT MEMORY AVAILABILITY -----
C
      MEMEXT = MEMAVL - MEMLIM
      IF (MASWRK) THEN
      IF (MEMEXT .NE. 0) THEN
         WRITE(IW,*) MEMLIM,' WORDS REQUESTED ',
     *               MEMAVL,' WORDS AVAILABLE ',
     *               MEMEXT,' WORDS REMAIN FREE'
      ELSE
         WRITE(IW,FMT='(1X,I12,'' WORDS OF MEMORY AVAILABLE'')') MEMLIM
      END IF
      END IF
C
C     SET /MACHIN/ VALUES.
C        MAXFM/SM  = TOTAL FAST/SLOW MEMORY IN DYNAMIC POOL.
C        LIMFM/SM  = TOP ADDRESS OF THE FAST/SLOW DYNAMIC POOL.
C     GAMESS DOES NOT CURRENTLY USE SLOW MEMORY.
C
      MAXFM = MEMLIM
      LIMFM = MEMLIM + LOFFS
      MAXSM = 0
      LIMSM = 0
C
      IPAR = MEMAVL
      RETURN
      END
C*MODULE UTIL  *DECK BIGFM
C>
C>    @brief   Release the total memory pool
C>
C>    @author  Steve Elbert, Mike Schmidt
C>
      SUBROUTINE BIGFM(IPAR)
      USE FAST_MEMORY, ONLY: MEMORY
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
C
C     ----- RETURN MAXIMUM MEMORY USED -----
C
      IF (MEMOK .NE. 99) RETURN
      MEMOK = 0
      IPAR = LENHI
      IF (MASWRK) WRITE(IW,*) IPAR,' WORDS OF DYNAMIC MEMORY USED'
      CALL FLSHBF(IW)
      CALL FLSHBF(IP)
C
C        ----- RELEASE MEMORY TO SYSTEM -----
C     STATIC IMPLEMENTATIONS SHOULDN'T DO ANYTHING HERE.
C
      IF(ALLOCATED(MEMORY)) THEN
         DEALLOCATE(MEMORY)
      END IF
      RETURN
      END
C*MODULE UTIL  *DECK VALFM
C>
C>    @brief   Return the current top of FM array
C>
C>    @author  Steve Elbert, Mike Schmidt
C>
      SUBROUTINE VALFM(IPAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign
C
C     ----- RETURN THE CURRENT TOP OF FM ARRAY -----
C
      IPAR = LTOP + LOFFS
      RETURN
      END
C*MODULE UTIL  *DECK GETFM
C>
C>    @brief   Reserve IPAR0 words of FM
C>
C>    @author  Steve Elbert, Mike Schmidt
C>
      SUBROUTINE GETFM(IPAR0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C         UNCOMMENT NEXT LINE AND ALSO IN -RETFM- FOR PLUMBING JOBS.
CLEAK IF(MASWRK) WRITE(IW,*) 'GETFM: ALLOCATING',IPAR
C
C     ----- RESERVE IPAR WORDS OF FM -----
C
      IPAR=((IPAR0-1)/nalign+1)*nalign
c     write(6,*) 'wwwa',IPAR0,IPAR,nalign
      LTOP = LTOP + IPAR
      IF (LTOP .LE. MEMLIM ) THEN
         LENHI = MAX(LENHI,LTOP)
      ELSE
         WRITE(IW,9000) ME,LTOP,MEMLIM
         CALL ABRT
      END IF
      RETURN
C
 9000 FORMAT(1X,'***** ERROR: MEMORY REQUEST EXCEEDS AVAILABLE MEMORY'/
     *    1X,'PROCESS NO.',I5,' WORDS REQUIRED=',I12,' AVAILABLE=',I12)
      END
C*MODULE UTIL  *DECK RETFM
C>
C>    @brief   Return IPAR0 words of FM
C>
C>    @author  Steve Elbert, Mike Schmidt
C>
      SUBROUTINE RETFM(IPAR0)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign
C
C         UNCOMMENT NEXT LINES AND ALSO IN -GETFM- FOR PLUMBING JOBS.
CLEAK LOGICAL GOPARR,DSKWRK,MASWRK
CLEAK COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
CLEAK COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
CLEAK IF(MASWRK) WRITE(IW,*) 'RETFM:  RETURNING',IPAR
C
C     ----- RETURN IPAR WORDS OF FM -----
C
      IPAR=((IPAR0-1)/nalign+1)*nalign
c     write(6,*) 'wwwb',IPAR0,IPAR,nalign
      LTOP = LTOP - IPAR
      RETURN
      END
C*MODULE UTIL  *DECK GOTFM
C>
C>    @brief   Return number of free words in FM
C>
C>    @author  Steve Elbert, Mike Schmidt
C>
      SUBROUTINE GOTFM(IPAR)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /FMPARM/ LTOP,LOFFS,LENHI,LOCMEM,MEMLIM,MEMOK,nalign
C
C     ----- RETURN NUMBER OF FREE WORDS IN FM -----
C
      IPAR0 = MEMLIM - LTOP
      IPAR=((IPAR0-1)/nalign+1)*nalign-nalign
      RETURN
      END
C*MODULE UTIL  *DECK TEXIT
C>
C>    @brief   Call ABRT if total CPU time has exceeded TIMLIM
C>
C>    @author  Mike Schmidt
C>
      SUBROUTINE TEXIT(NCALL,NREST)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL PACK,GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /RESTAR/ TIMLIM,IREST,NREC,INTLOC,IST,JST,KST,LST
      COMMON /OUTPUT/ NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /INTFIL/ NINTMX,NHEX,NTUPL,PACK,INTTYP,IGRDTYP
      COMMON /TMVALS/ TI,TX,TIM
C
      CALL TIMIT(NCALL)
      IF (TIM .LT. TIMLIM) RETURN
      IF (MASWRK) THEN
         WRITE(IW,9018)
         WRITE (IW,9008) TIMLIM,NPRINT,ITOL,ICUT,NORMF,NORMP,NOPK,
     *       NREST,IST,JST,KST,LST,NREC,INTLOC,NINTMX
      END IF
      CALL ABRT
      RETURN
C
 9008 FORMAT(F10.0,11I3,I10,2I5)
 9018 FORMAT(1X,'**** JOB HAS EXHAUSTED ITS CPU ALLOTMENT ****')
      END
C*MODULE UTIL  *DECK TIMIT
C>
C>    @brief   Compute and print interval CPU time
C>
C>    @author  Mike Schmidt
C>
      SUBROUTINE TIMIT(INDEX)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /TMVALS/ TI,TX,TIM
      COMMON /TIMING/ CPU,WALL
      COMMON /SIMDAT/ NACC,NREJ,IGOMIN,NRPA,IBWM,NACCT,NREJT,NRPAT,
     *                NPRTGO,IDPUNC,IGOFLG
C
C        COMPUTE AND PRINT INTERVAL CPU TIME
C        THIS IS CALLED WITH INDEX=0 ONLY FROM THE MAIN PROGRAM
C        ALL OTHER CALLS SHOULD PRINT THE INTERVAL TIME.
C
      CALL TSECND(TIM)
      TX = TIM-TI
      TI = TIM
      IF (INDEX .EQ. 0) RETURN
C
      RATIO=100.0D+00
      IF(WALL.GT.0.0D+00) RATIO = 100.0D+00*TIM/WALL
      TMINUT = TIM/60.0D+00
      IF((IPTIM.GT.0.OR.MASWRK).AND.NPRTGO.NE.2) THEN
         IF(GOPARR) THEN
            WRITE(IW,9010) ME,TX,TIM,TMINUT,WALL,RATIO
         ELSE
            WRITE(IW,9000) TX,TIM,TMINUT,WALL,RATIO
         END IF
      END IF
      CALL FLSHBF(IW)
      CALL FLSHBF(IP)
      RETURN
C
 9000 FORMAT(1X,'STEP CPU TIME =',F9.2,
     *          ' TOTAL CPU TIME =',F13.1,' (',F9.1,' MIN)'/
     *       1X,'TOTAL WALL CLOCK TIME=',F13.1,
     *          ' SECONDS, CPU UTILIZATION IS',F9.2,'%')
 9010 FORMAT(1X,'CPU',I6,': STEP CPU TIME=',F9.2,
     *          ' TOTAL CPU TIME=',F13.1,' (',F9.1,' MIN)'/
     *       1X,'TOTAL WALL CLOCK TIME=',F13.1,
     *          ' SECONDS, CPU UTILIZATION IS',F9.2,'%')
      END
C*MODULE UTIL  *DECK TMDATE
C>
C>    @brief   RETURN REAL ARRAY CONTAINING 24 CHARACTER TIME STAMP
C>
C>    @details ANY FORMAT RESEMBLING "HH:MM:SS DD-MMM-YYYY" WILL DO.
C>
C>    @param   TIMSTR - Hollerith string of 24 characters with time stamp
C>
      SUBROUTINE TMDATE(TIMSTR)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: TIMSTR(3)
C
      CHARACTER(LEN=3), DIMENSION(12), PARAMETER :: MONTHS =
     *                (/ 'JAN','FEB','MAR','APR','MAY','JUN',
     *                   'JUL','AUG','SEP','OCT','NOV','DEC' /)
      CHARACTER(LEN=24) :: CTIME
      INTEGER, DIMENSION(8) :: values
C
      call DATE_AND_TIME(VALUES=values)
      WRITE(CTIME,"(I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,A3,A1,I4.4)")
     *           values(5), ":", values(6), ":", values(7), " ",
     *           values(3), "-", MONTHS(values(2)), "-", values(1)
      READ(UNIT=CTIME,FMT='(3A8)') TIMSTR
      RETURN
      END
C*MODULE UTIL  *DECK TSECND
C>
C>    @brief   do CPU and wall clock timings
C>
C>    @details do CPU and wall clock timings through a name TSECND
C>             that is not used on any operating system, in order
C>             to call whatever name each kind of computer uses.
C>
C>    @param   TIM   return value is CPU time since job start
C>    @param   CPU   CPU time,   since job start (in /TIMING/)
C>    @param   WALL  wall clock, since job start (in /TIMING/)
C>
C>    @author  Mike Schmidt, many moons ago
C>
      SUBROUTINE TSECND(TIM)
C-----------------------------------------------------------------------
C       ----- THIS ROUTINE PERFORMS CPU AND WALL CLOCK TIMING -----
C
C       THIS ROUTINE SHOULD SET 'CPU' AND 'WALL' VARIABLES
C       TO THE TOTAL ELAPSED CPU AND WALL CLOCK TIMES,
C       MEASURED IN SECONDS.  IN ADDITION, THE CALLING
C       ARGUMENT 'TIM' SHOULD BE SET EQUAL TO 'CPU'.
C
C       ON THE FIRST ENTRY, 'CPU0' AND 'WALL0' SHOULD BE SET
C       TO THE APPROPRIATE BASE VALUE ON JOB START.
C
C-----------------------------------------------------------------------
C
C                       >>>>> UNIX NOTE <<<<<
C
C     DEPENDING ON JUST WHAT UNIX YOU ARE USING, ETIME IS EITHER
C     SINGLE OR DOUBLE PRECISION (DOUBLE ON THE CELERITY AND ALLIANT).
C     SOME SYSTEMS ALSO DON'T RETURN THE SUM OF USER AND SYSTEM
C     TIME AS THE FUNCTION VALUE.  THEREFORE, THE MOST PORTABLE
C     WAY OF DOING THIS TIMING IS TO ADD THE USER AND SYSTEM TIMES
C     RETURNED IN TARRAY, WHICH SEEMS TO BE SINGLE PRECISION ON
C     ALL UNIX SYSTEMS THUS FAR TESTED.
C
      IMPLICIT NONE
      DOUBLE PRECISION CPU0,WALL0,CPU,WALL,TIM
      INTEGER(8) :: IWALL, MAX_RATE
      LOGICAL FIRST
      COMMON /TIMING/ CPU,WALL
      SAVE FIRST,CPU0,WALL0,MAX_RATE
      DATA FIRST/.TRUE./
      IF(FIRST) THEN
         FIRST = .FALSE.
         CALL SYSTEM_CLOCK(COUNT=IWALL, COUNT_RATE=MAX_RATE)
         WALL0 = DBLE(IWALL) / MAX_RATE
         CALL CPU_TIME(CPU0)
      END IF
      CALL SYSTEM_CLOCK(COUNT=IWALL)
      WALL = DBLE(IWALL) / MAX_RATE
      WALL = WALL - WALL0
      CALL CPU_TIME(CPU)
      CPU = CPU - CPU0
      TIM = CPU
      RETURN
      END
C
C*MODULE UTIL  *DECK SETMXSEQMTX
C>
C>    @brief   Set matrix threshold sizes to force sequential
C>             N**2 and N**3 loops
C>
C>    @author  Andrey Asadchev
C>
      SUBROUTINE SETMXSEQMTX(MXSEQ2IN,MXSEQ3IN)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3
C
C        VERY SMALL MATRIX SIZES SHOULD NOT ATTEMPT TO RUN THE
C        LINEAR ALGEBRA STEPS IN PARALLEL.  THESE THREE ROUTINES
C        RETURN SOME NOTION OF WHAT "VERY SMALL" MIGHT MEAN.
C
      IF(MXSEQ2IN.GT.0) THEN
         MXSEQ2 = MXSEQ2IN
      ELSE
         MXSEQ2 = 300
      END IF
C
      IF(MXSEQ3IN.GT.0) THEN
         MXSEQ3 = MXSEQ3IN
      ELSE
         MXSEQ3 = 150
      END IF
C
      RETURN
      END
C*MODULE UTIL  *DECK MXSQN2
C>
C>    @brief   Returns matrix threshold size to force sequential N**2
C>             loops
C>
C>    @author  Andrey Asadchev
C>
      INTEGER FUNCTION MXSQN2()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3
C       RETURNS MATRIX THRESHOLD SIZE TO FORCE SEQUENTIAL N**2 LOOPS.
      MXSQN2 = MXSEQ2
      RETURN
      END
C*MODULE UTIL  *DECK MXSQN3
C>
C>    @brief   Returns matrix threshold size to force sequential N**3
C>             loops
C>
C>    @author  Andrey Asadchev
C>
      INTEGER FUNCTION MXSQN3()
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      COMMON /MTXSIZ/ MXSEQ2,MXSEQ3
C       RETURNS MATRIX THRESHOLD SIZE TO FORCE SEQUENTIAL N**3 LOOPS.
      MXSQN3 = MXSEQ3
      RETURN
      END
C*MODULE UNPORT  *DECK ABRT
C>
C>    @brief   Generate a calling sequence trace and stop execution
C>
C>    @author  ???
C>
C>    @date   Jul, 2022 Christian Friedl
C>    - Replace NAPTIME with COMMIT if possible.
C>
      SUBROUTINE ABRT
      USE UNPORT, ONLY: COMMIT, COMMIT_SUPPORTED
C
C     ----- GENERATE A CALLING SEQUENCE TRACE AND STOP EXECUTION -----
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION TIMSTR(3)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
C
C     ----- PRINT ACCOUNTING INFO AND GOODBYE MESSAGE -----
C
      call abortdump
      CALL TMDATE(TIMSTR)
      IF (MASWRK) WRITE(IW,900) TIMSTR
  900 FORMAT(1X,'EXECUTION OF GAMESS TERMINATED -ABNORMALLY- AT ',3A8)
      CALL BIGFM(MAXFM)
      CALL TIMIT(1)
C
C     STALL FOR A SHORT WHILE TO ALLOW A CHANCE OF PRINTING
C     ANY ERROR MESSAGES, BEFORE DDIKICK.X BLOWS US AWAY.
C     MAINLY THIS IS IN CASE A NON-MASTER PROCESS ENTERS FIRST,
C     MOMENTS BEFORE THE MASTER, WHICH SHOULD BE GIVEN A CHANCE
C     TO PRINT THE ERROR MESSAGE BEFORE THE NON-MASTER EXITS.
C     Reset bit 1 of MODIO to force flushing.
C
      if(iand(modio,1).ne.0) modio=modio-1
      CALL FLSHBF(IW)
      CALL FLSHBF(IP)
      IF (COMMIT_SUPPORTED) THEN
         CALL COMMIT(IW)
         CALL COMMIT(IP)
      ELSE
         IDELAY=1
         IF(ME.GT.0) IDELAY=10
         CALL NAPTIME(IDELAY)
      ENDIF
C
C     ----- EXIT PARALLEL RUNS GRACELESSLY -----
C
      CALL DDI_PEND(1)
C
C     ----- DO A TRACEBACK AND/OR CORE DUMP -----
C
      IF (MASWRK) THEN
         IF (ICORFL.EQ.1) THEN
            WRITE(IW,1) 'GDB','BACKTRACE','QUIT'
            CALL ABORT
         ELSE
            WRITE(IW,2)
         END IF
      END IF
C
C        ----- GENERIC STOP, IN CASE YOU GET THIS FAR -----
C
      STOP 'IN ABRT'
C
    1 FORMAT(1X,'TRACEBACK CAN BE OBTAINED BY LOCATING',
     *          ' THE ''CORE'' MEMORY DUMP FILE'/
     *   1X,'AND TYPING THE FOLLOWING (IN LOWER CASE)'/
     *   1X,A,' /U.../GAMESS/GAMESS.01.X CORE (INVOKES DEBUGGER)'/
     *   1X,'   ',A,'                         (PRINTS TRACEBACK)'/
     *   1X,'   ',A,'                         (QUITS DEBUGGER)')
    2 FORMAT(1X,'IF YOU WANT A CORE FILE, SET COREFL=.TRUE.',
     *          ' IN $SYSTEM.')
      END
C*MODULE UNPORT  *DECK VNAN
C>
C>    @brief   Initializes ARRAY to quiet form of "NOT A NUMBER"
C>
C>    @author  Mike Schmidt
C>
      SUBROUTINE VNAN(ARRAY,ISTRIDE,LENGTH)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION ARRAY(*)
C
      INTEGER NAN
      EQUIVALENCE (QNAN,NAN)
      DATA  NAN/Z'7FF8000000000001'/
C
C     ---- INITIALIZES -ARRAY- TO QUIET FORM OF "NOT A NUMBER" ----
C     SET -LENGTH- TOTAL ELEMENTS, SEPARATED BY INCREMENT -ISTRIDE-.
C     THIS IS MEANT FOR DEBUGGING PURPOSES, SINCE IT IS THE ANTI-VCLR.
C
C     WE CHOOSE THE SO-CALLED QUIET FORM OF NOT-A-NUMBER, WHOSE
C     IEEE REPRESENTATION REQUIRES SETTING ALL EXPONENT BITS, WITH
C     A NON-ZERO VALUE IN THE MANTISSA.  A 'QUIET NAN', AS OPPOSED
C     TO THE SIGNALLING-NOT-A-NUMBER (Z'7FF0000000000001'), HAS THE
C     MANTISSA'S SIGN BIT SET.  NOTE: INFINITY=Z'7FF0000000000000'.
C
C     USAGE IN THE MIDDLE OF A DODGY BIT OF CODE:
C               CALL VALFM(LOADFM)
C               LA   = LOADFM + 1
C               LB   = LA     + ...
C                   ...
C               LQ   = LP     + ...
C               LAST = LQ     + ...
C               NEED = LAST - LOADFM - 1
C               CALL GETFM(NEED)
C         C------------SET BAD VALUES--------
C               CALL VNAN(XX(LOADFM+1),1,NEED)
C         C------------SET BAD VALUES--------
C               ...DO SOMETHING WITH THE MEMORY HERE...
C               CALL RETFM(NEED)
C
      IF (ISTRIDE.EQ.1) THEN
         DO L=1,LENGTH
            ARRAY(L) = QNAN
         ENDDO
      ELSE
         LA=1-ISTRIDE
         DO L=1,LENGTH
            LA=LA+ISTRIDE
            ARRAY(LA) = QNAN
         ENDDO
      END IF
      RETURN
      END
C*MODULE UNPORT  *DECK ENDING
C>
C>    @brief   Terminate execution smoothly
C>
C>    @author  ???
C>
C>    @date   Jul, 2022 Christian Friedl
C>    - Replace NAPTIME with COMMIT if possible.
C>
C>    @date   Mar, 2023 Christian Friedl
C>    - Call DDI_DLBRESET before DDI_PEND to avoid
C>      garbled output from data servers.
C>
      SUBROUTINE ENDING
      USE UNPORT, ONLY: COMMIT, COMMIT_SUPPORTED
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      COMMON /IOFILE/ IR,IW,IP,IJK,IJKT,IDAF,NAV,IODA(950)
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     ----- TERMINATE EXECUTION SMOOTHLY -----
C     IF SUPPORTED, WAIT UNTIL DATA IS COMMITED OTHERWISE
C     WAIT A FIXED AMOUNT OF TIME TO TRY TO ALLOW ALL
C     BUFFERS TO FLUSH PROPERLY
C     Reset bit 1 of MODIO to force flushing.
C
      if(iand(modio,1).ne.0) modio=modio-1
      CALL FLSHBF(IW)
      CALL FLSHBF(IP)
      IF (COMMIT_SUPPORTED) THEN
         CALL COMMIT(IW)
         CALL COMMIT(IP)
      ELSE
         IDELAY=2
         CALL NAPTIME(IDELAY)
      END IF
C
C        CLEAN UP PARALLEL EXECUTION
C
      IF (GOPARR) CALL DDI_DLBRESET()
      CALL DDI_PEND(0)
      RETURN
      END
C
C*MODULE UNPORT  *DECK NAPTIME
C>
C> @brief   Pause execution for [seconds] seconds
C>
C> @details POSIX sleep and Windows'es Sleep are used
C>          POSIX sleep      expects time in      seconds
C>          Windows'es Sleep expects time in milliseconds
C>
C> @author Igor S. Gerasimov
C>
C> @date   Sep, 2020 Initial release
C>
C> @date   Sep 2022, Christian Friedl
C> - Shift platform dependency to parameter.
C> - Move subroutine out of module.
C>
      SUBROUTINE NAPTIME(seconds)
      use iso_c_binding, only: c_int
      use unport, only: sleep, seconds_to_delay
      implicit none
      integer, intent(in) :: seconds
      integer(c_int) :: delay, tmp
      delay = seconds * seconds_to_delay
      tmp = sleep(delay)
      END SUBROUTINE NAPTIME
C
C*MODULE UNPORT  *DECK BEGING
C>
C>    @brief   Take care of run initialization
C>
C>    @author  ???
C>
C>    @date   Sep, 2022 Christian Friedl
C>    - Shift platform dependency to parameters
C>
C*MODULE UNPORT  *DECK BEGING
      SUBROUTINE BEGING(VERSN)
C
C------------------------------------------------
C     ----- TAKE CARE OF RUN INITIALIZATION -----
C------------------------------------------------
C
      USE UNPORT, ONLY: DIRECTORY_SEPARATOR, VERSION_STRING,
     *                  DEBUG_WAIT_SUPPORTED, DEBUG_WAIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      LOGICAL GOPARR,DSKWRK,MASWRK
C
      CHARACTER*40 VERSN
      CHARACTER*1  DIRSEP
C
      COMMON /MACHIN/ NWDVAR,MAXFM,MAXSM,LIMFM,LIMSM
      COMMON /MACHSW/ KDIAG,ICORFL,IXDR,modio,mem10,lpnt10,mem10m
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /FILESY/ DIRSEP
C
C-----------------------------------------------------------------------
C     ----- SET -NWDVAR- IN COMMON /MACHIN/ -----
C     THIS IS THE NUMBER OF INTEGERS IN A WORKING PRECISION
C     FLOATING POINT NUMBER, WHICH IS MOST COMMONLY 64 BITS.
C
      NWDVAR = 1
      DIRSEP = DIRECTORY_SEPARATOR
C
      MAXFM=0
      MAXSM=0
C
C-----------------------------------------------------------------------
C
C INITIALIZE PARALLEL
C
      CALL DDI_PBEG(NWDVAR)
      CALL DDI_NPROC(NPROC,ME)
      MASTER = 0
      GOPARR = NPROC.GT.1
      MASWRK = ME.EQ.MASTER
      DSKWRK = .FALSE.
      IXDR   = 1
C-----------------------------------------------------------------------
C
C        ----- DEFINE THE MACHINE VERSION -----
C
      VERSN = VERSION_STRING
      IF(DEBUG_WAIT_SUPPORTED.AND.MASWRK) CALL DEBUG_WAIT()
      RETURN
      END
