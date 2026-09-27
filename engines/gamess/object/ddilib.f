c 17 Oct 19 - DGF - tweaks for GDDI/3 (DLB)
C  6 Jun 18 - DGF - tweaks for FMO 5.3
C 18 Apr 16 - SRP - tweaks for GDDI 
C 26 Feb 15 - SRP - Moved Comm destruction to separate routine, adjusted 
C                   GDDI common, and created GDDI_REGROUP subroutine
C xx Aug 14 - AG  - updated Doxygen comments for ML-GDDI
C  6 Feb 13 - SRP - delete old communicators when making new ones
C 23 MAR 12 - DGF - PAD COMMON BLOCKS
C 28 DEC 11 - DGF - REDUCE GDDI OUTPUT
C xx Jan 11 - AG  - introduced ML-GDDI functionality
C 11 AUG 11 - MWS/GDF - GDDI_INIT: SUPPRESS SCOPE CHANGE MESSAGES
C 11 AUG 10 - DGF - ADD GDDI_ASCOPE TO CHANGE SCOPE ASYNCHRONOUSLY
C 22 DEC 06 - BN,DGF  - VSCF SUBGROUP SUPPORT, REMOVE UNUSED ROUTINES
C  6 NOV 06 - MWS - ADJUST GDDI COMMON BLOCK
C  7 SEP 06 - RMO - CHANGE INTERFACE TO DDI GROUPS
C  5 JUL 05 - MWS - SELECT NEW ATOM,BASIS,EFP,PCM,DAF DIMENSIONS
C 30 APR 05 - DGF - SYNCHRONISE SYMBLK COMMON GDDI
C 19 MAY 04 - DGF - GAMESS-SPECIFIC INTERFACE WITH DDI ROUTINES
C
C*MODULE DDILIB   *DECK DDI_NSUMI
      SUBROUTINE DDI_NSUMI(MSGTAG, INT1, INT2, INT3, INT4, MSGLEN)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NBUF=4)
      INTEGER IBUF(NBUF)
C
C     SUM UP TO 4 INTEGERS IN ONE DDI_GSUMI
C
      IF(MSGLEN.GT.NBUF) CALL ABRTX("DDI_NSUMI error")
C
C
      IF(MSGLEN.GE.1) IBUF(1)=INT1
      IF(MSGLEN.GE.2) IBUF(2)=INT2
      IF(MSGLEN.GE.3) IBUF(3)=INT3
      IF(MSGLEN.GE.4) IBUF(4)=INT4
C
      CALL DDI_GSUMI(MSGTAG,IBUF,MSGLEN)
C
      IF(MSGLEN.GE.1) INT1=IBUF(1)
      IF(MSGLEN.GE.2) INT2=IBUF(2)
      IF(MSGLEN.GE.3) INT3=IBUF(3)
      IF(MSGLEN.GE.4) INT4=IBUF(4)
C
      RETURN
      END
C*MODULE DDILIB   *DECK GDDI_INIT
C>
C>  @brief   create desired subgroup partitioning of assigned nodes
C>
C>  @details create desired subgroup partitioning of assigned nodes
C>
C>  @author  Dmitri Fedorov and/or Ryan Olson
C>  @note Modified by Spencer Pruit to destroy old MPI communicators
C>  @note Modified by Alexander Gaenko to allow splitting at GROUP level     
C>  @note Modified by Spencer R. Pruitt to move Comm destruction to a new subroutine
C>        so that ML-GDDI actually works with NGRFMO
C>
      SUBROUTINE GDDI_INIT(NEWNGROUPS,MANNOD,some,isplit)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL ISGDDI,PAROUT,INITGDDI,some,some2,wasgddi,MLGDDI,
     *        GOPARR,DSKWRK,MASWRK
      INTEGER DDI_WORLD,DDI_GROUP,DDI_MASTERS,MANNOD(0:*)
      PARAMETER(DDI_WORLD=0,DDI_GROUP=1,DDI_MASTERS=2)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
      INTEGER COMM(DDI_WORLD:DDI_MASTERS)
      COMMON /GRPCOM/ COMM
C
      LOGICAL FIRST
      DATA FIRST /.TRUE./
      SAVE FIRST
c
      IF(NSUBGR.LT.0) THEN
c       calling scope change is needed for some variable setup, just in case
        CALL GDDI_SCOPE(DDI_WORLD)
        RETURN
      ENDIF
c
      CALL GDDI_SCOPE(DDI_WORLD)
c
      some2=some.and.(parout.or.maswrk)
      if(some2) WRITE(IW,9000) MEGLOB,NGROUPS,MYGROUP
C
      IF(.NOT.FIRST .AND. COMM(DDI_WORLD).NE.COMM(DDI_GROUP)) THEN
C      IF(isplit .EQ. 0) THEN
       DO I=DDI_GROUP,DDI_MASTERS
        CALL DDI_COMM_DESTROY(COMM(I))
       ENDDO
       first=.false.
      ENDIF
C
      IF(MANNOD(0).LT.0) THEN
        CALL DDI_GROUP_CREATE(NEWNGROUPS,COMM(0),COMM(1),COMM(2))
      ELSE
        CALL DDI_GROUP_CREATE_CUSTOM(NEWNGROUPS,MANNOD,
     &                               COMM(0),COMM(1),COMM(2))
      ENDIF
C
      NGROUPS=NEWNGROUPS
      CALL GDDI_SCOPE(DDI_GROUP)
      CALL DDI_NGROUP(NGROUPS,MYGROUP)
      if(some2) WRITE(IW,9010) MEGLOB,NGROUPS,MYGROUP
      if(NGROUPS.lt.0) write(6,*) isplit
C
      RETURN
 9000 FORMAT(/1X,'GDDI IS ABOUT TO SWITCH GROUPS: ',
     *        'MEGLOB=',I4,' NGROUPS=',I3,' MYGROUP=',I3)
 9010 FORMAT(/1X,'GDDI HAS SWITCHED GROUPS: ',
     *        'MEGLOB=',I4,' NGROUPS=',I3,' MYGROUP=',I3)
      END
C*MODULE DDILIB   *DECK GDDI_DESTROY
C>
C>  @brief   Destroy DDI_GROUP and DDI_MASTERS communicators
C>
C>  @details Required to avoid COMM build-up and for ML-GDDI
C>
C>  @author  Spencer R. Pruitt
C>
      SUBROUTINE GDDI_DESTROY
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c     LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      INTEGER DDI_WORLD,DDI_GROUP,DDI_MASTERS
      PARAMETER(DDI_WORLD=0,DDI_GROUP=1,DDI_MASTERS=2)
C
      INTEGER COMM(DDI_WORLD:DDI_MASTERS)
      COMMON /GRPCOM/ COMM
C
      CALL GDDI_SCOPE(DDI_WORLD)
c     WRITE(IW,9000)
C
C  SRP: DESTROYING OLD MPI COMMUNICATORS
C
      DO I=DDI_GROUP,DDI_MASTERS
         CALL DDI_COMM_DESTROY(COMM(I))
      ENDDO
C
      RETURN
c9000 FORMAT(/1X,'GDDI IS ABOUT TO DESTROY COMMUNICATORS ')
      END
C
C*MODULE DDILIB   *DECK GDDI_COMMID
      SUBROUTINE DDI_COMMID(DDI_SCOPE,COMMID)
      IMPLICIT NONE
C
      INTEGER DDI_SCOPE,COMMID
      INTEGER COMM(0:2)
      COMMON /GRPCOM/ COMM
C
      COMMID = COMM(DDI_SCOPE)
      RETURN
      END
C
C*MODULE DDILIB   *DECK GDDI_MASTID
      SUBROUTINE GDDI_MASTID(MASTID)
      IMPLICIT NONE
      INTEGER DDI_WORLD,DDI_GROUP
      LOGICAL ISGDDI,PAROUT,INITGDDI,GOPARR,DSKWRK,MASWRK,wasgddi,MLGDDI
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM,NSUBGR,MEUNIV,NPUNIV,numdlb,
     *        myworld,nworlds,mogddi
      PARAMETER(DDI_WORLD=0,DDI_GROUP=1)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      INTEGER NGROUPS,MYGROUP,ISCOPE,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *        MASTID(0:NGROUPS-1),I
C
C     FOR MASTID TO BE FILLED IN PROPERLY THE SCOPE MUST BE DDI_GROUP.
C     MASTID WILL CONTAIN THE GLOBAL ID OF MASTERS.
C
      DO I=0,NGROUPS-1
        MASTID(I)=0
        IF(I.EQ.MYGROUP.AND.ME.EQ.MASTER) MASTID(I)=MEGLOB
      ENDDO
      CALL GDDI_SCOPE(DDI_WORLD)
      CALL DDI_GSUMI( 2900, MASTID, NGROUPS )
      CALL GDDI_SCOPE(DDI_GROUP)
      RETURN
      END
C
C*MODULE DDILIB   *DECK GDDICOUNT
      SUBROUTINE GDDICOUNT(IMODE,LGROUP,MYJOB)
      IMPLICIT NONE
      INTEGER IMODE,MINE,NEXT,NGROUPS,LGROUP,MYGROUP,ISCOPE,MEGLOB,
     *        NPGLOB,NNGLOB,JBTYP,DDI_WORLD,DDI_GROUP,NPUNIV,MEUNIV,
     *        NSUBGR,numdlb,myworld,nworlds,mogddi
      PARAMETER(DDI_WORLD=0,DDI_GROUP=1)
      LOGICAL MYJOB,NXT,ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      SAVE MINE,NEXT
C
C     STATIC AND DYNAMIC LOAD BALANCING BETWEEN GROUPS FOR GDDI.
C     NOTE THAT WORK DIVISION INSIDE A GROUP IS DONE ELSEWHERE.
C     DETERMINE IF THE CURRENT TASK MUST BE DONE BY THIS NODE.
C     BEFORE USAGE CALL WITH IMODE=-1 TO INITIALISE. THEN USE WITH
C     IMODE=0 FOR GETTING JOB STATUS. AFTER USAGE CALL WITH IMODE=1.
C     ON EXIT MYJOB IS TRUE OR FALSE, LGROUP WILL POINT WHAT GROUP WILL
C     DO THIS JOB IF STATIC LOAD BALANCING IS IN EFFECT.
C     HOPEFULLY IT IS UNDERSTOOD THAT CALLS TO GDDICOUNT CAN'T BE NESTED
C     (ONLY ONE LOOP CAN BE DIVIDED INTO GROUPS AT A TIME).
C
      NXT = JBTYP.EQ.1
      IF(ABS(IMODE).EQ.1) THEN
        NEXT = -1
        MINE = -1
        LGROUP=-1
C       NEVER, NEVER ADD .AND. GOPARR TO NXT FOR GDDI!
C       SOME GROUPS MAY HAVE ONE NODE AND GOPARR=.FALSE.
        IF(NXT) THEN
C        ACCORDING TO THE REST OF GAMESS WE CALL DDI_DLBRESET AT THE
C
C        END (NOT AT THE BEGINNING, STRANGE AS IT MAY BE).
C        CHANGING SCOPES TO THE WORLD BACK IS NECESSARY, BECAUSE
C        OTHERWISE SYNCING WILL BE DONE IN THE SAME GROUP, THUS
C        LEADING TO PROCESSES RESETTING AND ASKING FOR DLB OUT OF SYNC.
C
           IF(IMODE.EQ.1) THEN
             CALL GDDI_SCOPE(DDI_WORLD)
c            An ugly nasty hack on exigent demand.
c            DLB for GDDI/3 does not work using standard means,
c            use a haggy hack.
             if(nsubgr.gt.0) then
c              write(6,*) 'wwwbbbres',idds(2)
               CALL DDI_GDLBRESET_INNER()
c              write(6,*) 'wwwbbbres done'
             else if(nsubgr.lt.0) then
               CALL DDI_DLBRESET
             else
               CALL DDI_GDLBRESET()
             endif
             CALL GDDI_SCOPE(DDI_GROUP)
           ENDIF
        ENDIF
        MYJOB=.TRUE.
      ELSE
        MINE=MINE+1
        IF(NXT) THEN
C         DYNAMIC LOAD BALANCING
c         IF(MINE/numdlb.GT.NEXT) CALL DDI_GDLBNEXT(NEXT)
          IF(MINE/numdlb.GT.NEXT) then
             if(nsubgr.gt.0) then
               CALL DDI_GDLBNEXT_INNER(NEXT)
             else if(nsubgr.lt.0) then
               CALL DDI_DLBNEXT(NEXT)
             else
               CALL DDI_GDLBNEXT(NEXT)
             endif
c           write(6,*) 'wwwGDLBNEXT',NEXT,MINE
          endif
C         -1 SHOWS THAT THE GROUP RESPONSIBLE IS UNKNOWN.
          LGROUP=-1
c         write(6,*) 'wwwtest',NEXT,MINE,numdlb
          MYJOB=NEXT.EQ.MINE/numdlb
        ELSE
C         STATIC LOAD BALANCING
          LGROUP=MOD(MINE,NGROUPS)
          MYJOB=LGROUP.EQ.MYGROUP
        ENDIF
      ENDIF
      RETURN
      END
C
C*MODULE DDILIB   *DECK GDDI3COUNT
C>
C>  @brief   A clone of GDDICOUNT.
C>
C>  @details Handle load balancing in GDDI/3.
C>
C>  @author  unknown
C>
      SUBROUTINE GDDI3COUNT(IMODE,LGROUP,MYJOB)
      IMPLICIT NONE
      INTEGER IMODE,MINE,NEXT,NGROUPS,LGROUP,MYGROUP,ISCOPE,MEGLOB,
     *        NPGLOB,NNGLOB,JBTYP,DDI_WORLD,DDI_GROUP,NPUNIV,MEUNIV,
     *        NSUBGR,numdlb,myworld,nworlds,mogddi,numdlb3
      PARAMETER(DDI_WORLD=0,DDI_GROUP=1)
      LOGICAL MYJOB,NXT,ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      SAVE MINE,NEXT
C
C     A clone of GDDICOUNT, for the use at the universe level of GDDI/3.
c     Here, the standard GDDI/2 routines without INNER are called; whereas at 
c     the world level in GDDI/3, *INNER should be called via GDDICOUNT.
C
c     Disable batches at the universe level (use batches at the world level).
c
      numdlb3=1
c     NXT = JBTYP.EQ.1
      NXT=iand(mogddi,2).ne.0
      IF(ABS(IMODE).EQ.1) THEN
        NEXT = -1
        MINE = -1
        LGROUP=-1
C       NEVER, NEVER ADD .AND. GOPARR TO NXT FOR GDDI!
C       SOME GROUPS MAY HAVE ONE NODE AND GOPARR=.FALSE.
        IF(NXT) THEN
C        ACCORDING TO THE REST OF GAMESS WE CALL DDI_DLBRESET AT THE
C
C        END (NOT AT THE BEGINNING, STRANGE AS IT MAY BE).
C        CHANGING SCOPES TO THE WORLD BACK IS NECESSARY, BECAUSE
C        OTHERWISE SYNCING WILL BE DONE IN THE SAME GROUP, THUS
C        LEADING TO PROCESSES RESETTING AND ASKING FOR DLB OUT OF SYNC.
C
           IF(IMODE.EQ.1) THEN
             CALL GDDI_SCOPE(DDI_WORLD)
             CALL DDI_GDLBRESET()
             CALL GDDI_SCOPE(DDI_GROUP)
           ENDIF
        ENDIF
        MYJOB=.TRUE.
      ELSE
        MINE=MINE+1
        IF(NXT) THEN
C         DYNAMIC LOAD BALANCING
          IF(MINE/numdlb3.GT.NEXT) CALL DDI_GDLBNEXT(NEXT)
C         -1 SHOWS THAT THE GROUP RESPONSIBLE IS UNKNOWN.
          LGROUP=-1
          MYJOB=NEXT.EQ.MINE/numdlb3
        ELSE
C         STATIC LOAD BALANCING
          LGROUP=MOD(MINE,NGROUPS)
          MYJOB=LGROUP.EQ.MYGROUP
        ENDIF
      ENDIF
      RETURN
      END
C
C*MODULE DDILIB   *DECK GDDI_COUNT
C>
C>  @brief   A switcher for calling GDDICOUNT.
C>
C>  @details Handle load balancing in GDDI/n.
C>
C>  @author  Dmitri Fedorov
C>
      SUBROUTINE GDDI_COUNT(IMODE,LGROUP,MYJOB)
      IMPLICIT NONE
      INTEGER IMODE,MINE,NEXT,NGROUPS,LGROUP,MYGROUP,ISCOPE,MEGLOB,
     *        NPGLOB,NNGLOB,JBTYP,NPUNIV,MEUNIV,NSUBGR,numdlb,myworld,
     *        nworlds,mogddi
      LOGICAL MYJOB,NXT,ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
c
c     For GDDI/3, this subroutine should be called at the universe level,
c     and at the world level, GDDICOUNT should be called.
c
      if(nsubgr.gt.0) then
        CALL GDDI3COUNT(IMODE,LGROUP,MYJOB)
      else
        CALL GDDICOUNT(IMODE,LGROUP,MYJOB)
      endif
      RETURN
      END
C
C*MODULE DDILIB   *DECK GDDI_SCOPE
C>    @brief Switches parallelization scope or parallelization level
C>
C>    @author Ryan Olsen and/or Graham Fletcher
C>    @note Modified by Alexander Gaenko to introduce SUPERWORLD and SUBGROUP scopes
      SUBROUTINE GDDI_SCOPE(NEWSCOPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,
     *        MLGDDI
      INTEGER DDI_MASTERS
      PARAMETER(DDI_MASTERS=2)
      integer DDI_GROUP
      parameter (DDI_GROUP=1)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
C
      INTEGER COMM(0:2)
      COMMON /GRPCOM/ COMM
C
      if(NSUBGR.lt.0) then
        call NDDI_SCOPE(NEWSCOPE,0)
        return
      endif
C           POTENTIALLY USEFUL FOR DEBUGGING SUBGROUPS SHIFTS.
C     WRITE(6,1) ISCOPE,NEWSCOPE
C   1 FORMAT(1X,'GDDI_SCOPE: CHANGING SCOPE FROM',I4,' TO ', I4)
C--   CALL FLSHBF(6)
C
      call GDDI_SUBSUPER(newscope,iscope)
      CALL DDI_SCOPE(COMM(ISCOPE))
      CALL DDI_NPROC( NPROC, ME )
C
C     NGROUPS HAS THE TOTAL NUMBER OF GROUPS,
C     NOT THIS NUMBER WITHIN THE CURRENT SCOPE
C     MASWRK=IMASTER.NE.0
C
C  NOTE: DDI_NGROUP WAS NEVER USED.  THE FUNCTIONALITY HAS CHANGED, IT
C  NOW REPORTS THE TOTAL NUMBER OF GROUPS AND THE RANK OF THE GROUP OF
C  THE CALLING PROCESS.  IN THE DDI_MASTERS SCOPE, THERE ARE ALWAYS TWO
C  GROUPS: MYGRP.EQ0 = GROUP MASTERS, MYGRP.EQ.1 = ALL OTHER PROCESSES
C  THE CONVENTION IS TO USE THE LOGICAL VARIABLE MASWRK TO DECIDE WHICH
C  PROCESSES SHOULD CONTINUE ON AND PERFORM COMMUNICATION CALL AND WHICH
C  PROCESSES SHOULD BYPASS THE COMMUNICATION.
C
      MASWRK=ME.EQ.MASTER
      IF(ISCOPE.EQ.DDI_MASTERS) THEN
         CALL DDI_NGROUP(NGRP,MYGRP)
         MASWRK=MYGRP.EQ.MASTER
      else if (ISCOPE.EQ.DDI_GROUP) then
C        update number of groups and group ID
         CALL DDI_NGROUP(NGROUPS,MYGROUP)
      END IF
C
      GOPARR=NPROC.GT.1
      RETURN
      END

C
C*MODULE DDILIB   *DECK GDDI_ASCOPE
      SUBROUTINE GDDI_ASCOPE(NEWSCOPE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,
     *        MLGDDI
      INTEGER DDI_MASTERS
      PARAMETER(DDI_MASTERS=2)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
C
C     ASYNCHRONOUS CLONE OF GDDI_SCOPE (DO NOT RUN DDI_SYNC WHILE
C     SWITCHING SCOPES IN DDI_ASCOPE). 
C     DDI_NPROC AND DDI_NGROUP SHOULD NOT CALL DDI_SYNC, OR ELSE
C     THEY CAN PROBABLY BE COMMENTED OUT BELOW.
C
      INTEGER COMM(0:2)
      COMMON /GRPCOM/ COMM
C
      if(NSUBGR.lt.0) then
c       for intranode GDDI sync or async scope change are the same,
c       as there is no synchronisation.
        call NDDI_SCOPE(NEWSCOPE,1)
        return
      endif
c
      call GDDI_SUBSUPER(newscope,iscope)
c     ISCOPE=NEWSCOPE
      CALL DDI_ASCOPE(COMM(ISCOPE))
      CALL DDI_NPROC( NPROC, ME )
C
      MASWRK=ME.EQ.MASTER
      IF(ISCOPE.EQ.DDI_MASTERS) THEN
         CALL DDI_NGROUP(NGRP,MYGRP)
         MASWRK=MYGRP.EQ.MASTER
      END IF
C
      GOPARR=NPROC.GT.1
      RETURN
      END

C A.G.: Additions for multilevel grouping.
      
C>   @brief Initializes COMMON block variables used for multilevel grouping
C>      
C>   @author Alexander Gaenko
      Block Data GDDI_SPLIT_BDATA
      implicit none
      integer Level, LenCommStack
      common /COMVAR/ Level, LenCommStack
      data Level /0/, LenCommStack /0/
      End
      
C>    @brief Splits group communicator, creating new level of parallelism
C>
C>    @param[IN] newgroups: number of new subgroups
C>    @param[IN] mannod: array of nodes/group or (-1) for default equal distribution
C>
C>    @details Splits the current communicator (must be a GROUP communicator) and
C>             creates a new communicator triad (WORLD, GROUP, MASTERS);
C>             the new WORLD communicator is the old GROUP after the splitting.
C>
C>    @author Alexander Gaenko      
      subroutine GDDI_SPLIT(newgroups,mannod)
      implicit none
C     Arguments:
      integer newgroups, mannod(0:*)

C     Common blocks:
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,mlgddi
      integer ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,nsubgr,
     *        meuniv,npuniv,numdlb,myworld,nworlds,mogddi
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA (950)
C
      INTEGER DDI_WORLD,DDI_GROUP,DDI_MASTERS
      parameter (DDI_WORLD=0, DDI_GROUP=1, DDI_MASTERS=2)
      INTEGER COMM(DDI_WORLD:DDI_MASTERS)
      COMMON /GRPCOM/ COMM

C     Stack of communicators (arrays are in separate commons for robustness)
      integer MAXLEVEL
      parameter (MAXLEVEL=2048)
      integer CommStack, MCommStack, Level, LenCommStack
      common /COMSTK/ CommStack(0:MAXLEVEL-1)
      common /MCOMST/ MCommStack(0:MAXLEVEL-1)
      common /COMVAR/ Level, LenCommStack
C     Local variables
c     integer mynode
      
      if (iscope .ne. DDI_GROUP) then
         write (IW,*) 'GDDI_SPLIT(): splitting in scope=',iscope,
     $        ' is unsupported'
         call Abrt
      endif

      if (Level+2 .ge. MAXLEVEL) then
         write (IW,*) 'GDDI_SPLIT(): too many levels: Level=',Level
         call Abrt
      endif
         
      CommStack(Level)=Comm(DDI_WORLD)
      MCommStack(Level)=Comm(DDI_MASTERS)
      Level=Level+1
      CommStack(Level)=Comm(DDI_GROUP)
C     Trick to make GDDI_INIT to split the group
      Comm(DDI_WORLD)=Comm(DDI_GROUP)
      call GDDI_Init(newgroups,mannod,.true.,0)
C     Verify that splitting occured as predicted: new WORLD is old GROUP
      if (Comm(DDI_WORLD) .ne. CommStack(Level)) then
         write (IW,*) 'GDDI_SPLIT(): unexpected splitting, '//
     $        'NEW_WORLD<>OLD_GROUP'
         call Abrt
      endif
      CommStack(Level+1)=Comm(DDI_GROUP)
      MCommStack(Level)=Comm(DDI_MASTERS)
      LenCommStack=Level+2
C     Now CommStack(Level-1) is SUPERWORLD, MCommStack(Level-1) is the superworld's MASTERS,
C         CommStack(Level) is WORLD (the superworld's GROUP), MCommStack(Level) is MASTERS,
C         CommStack(Level+1) is GROUP,  MCommStack(Level+1) is undefined,
C         LenCommStack is the length of the initialized CommStack.
      call GDDI_Scope(DDI_WORLD)
      call DDI_NProc(NPGlob,MeGlob)
C???      CALL DDI_NNODE(NNGLOB,mynode)
c
c     Here, in the world scope of the split universe, each world's
c     master (0) broadcasts its data server ID idds(1), saved as idds(2)
c     on all inhabitants of the world. 
c
c     idds(2)=idds(1)
c     CALL DDI_BCAST(350,'I',idds(2),1,0)
c     call abrt
      call DDI_get_dsid(1)
c     write(6,*) 'world master id is',idds(1),idds(2)
      return
      End
      
C>    @brief Helper subroutine: Switches parallelization levels, returns new scope
C>
C>    @param swscope(IN): level or scope to switch to (SUPERWORLD,SUBGROUP,GROUP,WORLD,MASTERS)
C>    @param iscope(IN): current scope (GROUP,WORLD,MASTERS)
C>    @param iscope(OUT): scope to switch to in the new level (GROUP,WORLD,MASTERS)
C>      
C>    @details This is not supposed to be called outside ddilib.src.
C>             This helper subroutine processes SUPERWORLD and SUBGROUP pseudo-scopes
C>             and changes the parallelization level by switching to the new
C>             (WORLD, GROUP, MASTERS) triad.
C>
      subroutine GDDI_SUBSUPER(swscope,iscope)
      implicit none
C     Arguments:
      integer swscope
      integer iscope
C     Common blocks:
C     Scope constants and communicators
      INTEGER DDI_WORLD,DDI_GROUP,DDI_MASTERS
      integer DDI_SUPERWORLD,DDI_SUBGROUP
      parameter (DDI_WORLD=0, DDI_GROUP=1, DDI_MASTERS=2,
     $     DDI_SUPERWORLD=-1, DDI_SUBGROUP=4)
      INTEGER COMM(DDI_WORLD:DDI_MASTERS)
      COMMON /GRPCOM/ COMM
C     GDDI commons:
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      integer ISCOPE0,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,MEUNIV,
     *        NPUNIV,NSUBGR,numdlb,myworld,nworlds,mogddi
      COMMON /GDDI/ISCOPE0,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *             ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *             MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
C     Stack of communicators (arrays are in separate commons for robustness)
      integer MAXLEVEL
      parameter (MAXLEVEL=2048)
      integer CommStack,MCommStack,Level,LenCommStack
      common /COMSTK/ CommStack(0:MAXLEVEL-1)
      common /MCOMST/ MCommStack(0:MAXLEVEL-1)
      common /COMVAR/ Level, LenCommStack
C     GAMESS I/O
      integer IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      
      if (swscope .eq. DDI_SUPERWORLD) then
         if (iscope .eq. DDI_GROUP) then
            iscope=DDI_WORLD
            return
         endif
         if (iscope .eq. DDI_MASTERS) then
            write (IW,*) 'GDDI_UPDOWN_SCOPE(): ',
     $           'Cannot upswitch MASTERS'
            call Abrt
         endif
C        Otherwise, we're switching up from WORLD.
         if (Level .lt. 1) then
            write (IW,*) 'GDDI_UPDOWN_SCOPE(): ',
     $           'Cannot switch up, level=', Level
            call Abrt
         endif
         Level=Level-1
         iscope=DDI_GROUP
c        write(6,*) 'wwwlevup',Level,iscope
         
      else if (swscope .eq. DDI_SUBGROUP) then
         if (iscope .eq. DDI_WORLD) then
            iscope=DDI_GROUP
            return
         endif
         if (iscope .eq. DDI_MASTERS) then
            write (IW,*) 'GDDI_UPDOWN_SCOPE(): ',
     $           'Cannot downswitch MASTERS'
            call Abrt
         endif
C        Oterwise, we're switching down from GROUP.
         if (Level .ge. LenCommStack-2) then
            write (IW,*) 'GDDI_UPDOWN_SCOPE(): ',
     $           'Cannot switch down, level=', Level,
     $           ' LenCommStack=', LenCommStack
            call Abrt
         endif
         Level=Level+1
         iscope=DDI_WORLD
c        write(6,*) 'wwwlevdown',Level,iscope
      else
C        It is a traditional "explicit" switching, just assign iscope:
         iscope=swscope
         return
      endif
C     set scope communicators for the new WORLD
      Comm(DDI_WORLD)=CommStack(Level)
      Comm(DDI_GROUP)=CommStack(Level+1)
      Comm(DDI_MASTERS)=MCommStack(Level)
C     update NpGlob and MeGlob for the new WORLD
      call DDI_NProc_comm(Comm(DDI_WORLD),NpGlob,MeGlob)
c     write(6,*) 'wwwlevels',NpGlob,MeGlob
      return
      End
c

C*MODULE DDILIB   *DECK GDDI_REGROUP
C>
C>  @brief   Destroy and regroup NGROUP,NSUBGR, or both
C>
C>  @details Allows for dynamic regrouping at either level
C>           similar to NGRFMO. 
C>
C>  @author  Spencer R. Pruitt
C>
!>  @param NEWGRS : new group division
!>
!>  @param NEWSUB : new subgroup divisions
!>
!>  @param NEWMAN1 : new manual group division
!>
!>  @param NEWMAN2 : new manual subgroup division
!>
C>  @date    Feb 26, 2015 - Spencer R. Pruitt
C>           Created subroutine
C>
      SUBROUTINE GDDI_REGROUP(NEWGRS,NEWSUB,NEWMAN1,NEWMAN2)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER DDI_WORLD,DDI_GROUP,DDI_SUBGROUP,DDI_SUPERWORLD
c     PARAMETER(MXUNIT=299,MAXGR=16384)
      parameter (DDI_WORLD=0, DDI_GROUP=1, DDI_SUPERWORLD=-1,
     *          DDI_SUBGROUP=4)
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,GOPARR,DSKWRK,MASWRK
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
c     COMMON /OPNNFT/ NFTOPN(MXUNIT),NODEXT(MXUNIT),IOSMP(2)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK

c     DIMENSION NEWMAN(0:MAXGR-1)

      !ROUTINE TO ALLOW EASY DESTRUCTION AND RECREATION OF 
      !GDDI GROUPS (MPI COMMUNICATORS)

      !FOUR PASSED VARIABLES MAY BE SET:
      !NEWGRS - NEW NUMBER OF TOP LEVEL GDDI GROUPS 
      !NEWSUB - NEW NUMBER OF ML-GDDI SUBGROUPS PER GDDI GROUP
      !NEWMAN1 - MANNOD FOR TOP LEVEL GDDI GROUPS
      !NEWMAN2 - MANNOD FOR ML-GDDI SUBGROUPS

      !IT IS POSSIBLE TO ONLY REGROUP AT THE LOWEST LEVEL
      !BY CALLING WITH NEWGRS & NEWMAN = -1

      !IF ONLY REGROUPING AT THE TOP LEVEL GDDI
      !CALL WITH NEWSUB & NEWMAN2 = -1

      !YOU CAN ALSO REGROUP EVERYTHING

      IF(NEWGRS.LT.0.AND.NEWSUB.LT.0) THEN
         WRITE(IW,8000)
         RETURN
      ENDIF

      CALL GDDI_SCOPE(DDI_WORLD)

      IF(NEWGRS.LT.0) THEN
      !ONLY REGROUPING SUBGROUPS, SCOPE DOWN
         CALL GDDI_SCOPE(DDI_GROUP)
         CALL GDDI_SCOPE(DDI_SUBGROUP)
      ENDIF

      MYGROUP = 0

      !WE WANT TO CREATE OUT OWN DYNAMIC GDDI GROUPING/SUBGROUPING
      !FOR EXAMPLE, THIS CAN BE BASED ON THE NUMBER OF FRAGMENTATION
      !STATES WE HAVE FOR REACTIVE MOLECULAR DYNAMICS.

      !FIRST DESTROY THE OLD DDI_GROUP AND DDI_MASTERS AT THE TOP LEVEL

      CALL GDDI_DESTROY

      !CREATE NGROUPS=NEWGRS
      !NUMBER OF GROUPS FOR FMOX IS CONTROLLED BY NGRFMO OR MANNOD (READ LATER)

      IF(NEWGRS.GT.0) THEN

         CALL GDDI_INIT(NEWGRS,NEWMAN1,.FALSE.,0)

         IF(MASWRK.AND..NOT.PAROUT) THEN

      !SEQOPN WILL NOT OPEN FILES IF THEY ARE ALREADY OPENED.

            CALL SEQOPN(IR,'INPUT', 'OLD',.TRUE., 'FORMATTED')
            CALL SEQOPN(IW,'OUTPUT','UNKNOWN',.FALSE., 'FORMATTED')
            CALL SEQOPN(IP,'PUNCH', 'NEW',.FALSE.,'FORMATTED')
         ENDIF

         NGPROC = INT(NPUNIV / NEWGRS)

         IF(MASWRK.OR.PAROUT) WRITE(IW,9000) NEWGRS,NGPROC

      ENDIF

      !CREATE SUB-GROUPS IF USING ML-GDDI

      IF(NEWSUB.GT.0) THEN
         CALL GDDI_SPLIT(NEWSUB,NEWMAN2)
         NSGPROC = INT(NGPROC / NEWSUB)
         IF(MASWRK.OR.PAROUT) WRITE(IW,9001) NEWSUB,NSGPROC
         CALL GDDI_SCOPE(DDI_WORLD)
         CALL GDDI_SCOPE(DDI_SUPERWORLD)
      ENDIF

      CALL GDDI_SCOPE(DDI_WORLD)

      RETURN
 8000 FORMAT(//1X,55(1H-)/,
     * 1X,"GDDI_REGROUP: YOU SHOULDN'T BE HERE...",/
     * 1X,"NEWGRS AND NEWSUB BOTH EQUAL -1",1X/
     * 1X,"MAYBE A BUG, MAYBE YOU DON'T KNOW WHAT YOU'RE DOING",1X/
     * 1X,55(1H-)/)
 9000 FORMAT(//1X,55(1H-)/,
     * 1X,'GDDI: REGROUPING OF NODES SUCCESSFUL, NGROUP = ',I6,/
     * 1X,'      WITH ',I5,' PROCESSES PER GROUP',1X/
     * 1X,55(1H-)/)
 9001 FORMAT(//1X,55(1H-)/,
     * 1X,'ML-GDDI: SUBGROUPING SUCCESSFUL, NSUBGR = ',I6,/
     * 1X,'         WITH ',I5,' PROCESSES PER SUB-GROUP',1X/
     * 1X,55(1H-)/)
      END
C*MODULE DDILIB   *DECK DDI_AOP
C>
C>  @brief   Asynchronious operation.
C>
C>  @details Do operations asynchronously.
C>
C>  @author  Dmitri Fedorov
C>
      SUBROUTINE DDI_AOP(OP,HANDLE,ILO,IHI,JLO,JHI,BUFF)
      IMPLICIT NONE
C     Declarations of arguments
      CHARACTER*3 OP
      INTEGER HANDLE,ILO,IHI,JLO,JHI
      DOUBLE PRECISION BUFF(*)
C     Declarations for common blocks
      LOGICAL ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI
      INTEGER ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP
      INTEGER NSUBGR,MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      INTEGER IR,IW,IP,IS,IPK,IDAF,NAV,IODA
C     Other declarations
      integer ddi_world
      integer ITMP
      Parameter(ddi_world=0)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
C
C     Do an asynchronious DDI operation (GET, PUT or ACC)
c     in the world scope
C
      ITMP=ISCOPE
      IF(ISGDDI) CALL GDDI_ASCOPE(DDI_WORLD)
      IF(OP.EQ.'PUT') THEN
        CALL DDI_PUT(HANDLE,ILO,IHI,JLO,JHI,BUFF)
      ELSE IF(OP.EQ.'GET') THEN
        CALL DDI_GET(HANDLE,ILO,IHI,JLO,JHI,BUFF)
      ELSE IF(OP.EQ.'ACC') THEN
        CALL DDI_ACC(HANDLE,ILO,IHI,JLO,JHI,BUFF)
      ELSE
        WRITE(IW,*) 'UNKNOWN OPERATION',OP,' FOR',HANDLE
        CALL ABRT
      ENDIF
c     write(6,*) 'wwwop',iscope,OP,HANDLE,ILO,IHI,JLO,JHI,BUFF(1)
      IF(ISGDDI.AND.ITMP.NE.DDI_WORLD) CALL GDDI_ASCOPE(ITMP)
      RETURN
      END
C
C*MODULE DDILIB   *DECK NDDI_SCOPE
C>
C>  @brief   Change scope in DDI.
C>
C>  @details Change Intranode scope.
C>
C>  @author  Dmitri Fedorov
C>
      SUBROUTINE NDDI_SCOPE(NEWSCOPE,MODE)
      IMPLICIT NONE
C     Declarations of arguments
      INTEGER NEWSCOPE,MODE
C     Declarations for common blocks
      LOGICAL GOPARR,DSKWRK,MASWRK,ISGDDI,PAROUT,INITGDDI,wasgddi,
     *        MLGDDI
      INTEGER ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP
      INTEGER NSUBGR,MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      INTEGER IR,IW,IP,IS,IPK,IDAF,NAV,IODA
      INTEGER ME,MASTER,NPROC,IBTYP,IPTIM
C     Other declarations
      INTEGER DDI_WORLD,DDI_GROUP,DDI_MASTERS
      PARAMETER(DDI_WORLD=0,DDI_GROUP=1,DDI_MASTERS=2)
      COMMON /GDDI/   ISCOPE,NGROUPS,MYGROUP,MEGLOB,NPGLOB,NNGLOB,JBTYP,
     *                ISGDDI,PAROUT,INITGDDI,wasgddi,MLGDDI,NSUBGR,
     *                MeUniv,NPUniv,numdlb,myworld,nworlds,mogddi
      COMMON /IOFILE/ IR,IW,IP,IS,IPK,IDAF,NAV,IODA(950)
      COMMON /PAR   / ME,MASTER,NPROC,IBTYP,IPTIM,GOPARR,DSKWRK,MASWRK
C
C     Intranode GDDI. The scopes are logical, not real (C library is not used). 
c     It is used to split nodes into groups without defining logical nodes,
c     accomplished as GDDI with NGROUP=1 and NSUBGR=-total number of cores.
c     For example, for a node "node1" with 20 cores, run as node1:cores=20
c     $GDDI NGROUP=1 NSUBGR=-20 $END
c     At present, NSUBGR should always be equal to 20 (for a 20-core node)!
c     Now, if running on 2 20-core nodes, use 
c     node1:cores=20 node2:cores=20 and NGROUP=1 NSUBGR=-40.
c     A negative NSUBGR is used to invoke intranode GDDI.
c     MODE=0 do a SYNC to have all ranks here
c     MODE=1 do no SYNC 
c
      if(NEWSCOPE.eq.DDI_WORLD.or.NEWSCOPE.eq.DDI_MASTERS) then
c     As stated above, all groups have 1 core by construction, i.e.,
c     all cores are masters and there is no difference between world and
c     master scopes. The 1-core setup is required because no intragroup
c     parallel operations are permitted.
c
        nproc=NPGLOB
        me=MEGLOB
      else if(NEWSCOPE.eq.DDI_GROUP) then
        nproc=1
        me=0
      else
c     What the charles?
        write(iw,*) 'Unknown scope',NEWSCOPE
        call abrt
      endif
c     MYGROUP is the "GDDI" group, no matter what the scope be.
c     Likewise, NGROUPS. From the first call to NDDI_SCOPE, GDDI is 
c     "masqueraded" into thinking we have not 1 but NPGLOB groups!
      MYGROUP=MEGLOB
      NGROUPS=NPGLOB
      ISCOPE=NEWSCOPE
c
      goparr=nproc.gt.1
      maswrk=me.eq.MASTER.or.NEWSCOPE.eq.DDI_MASTERS
c     For the master scope, all ranks are masters...
c
c     It is absolutely critical that all DDI calls such as DDI_GSUM,
c     DDI_BCAST and DLB are called ONLY for GOPARR=.false. If not, a deadlock
c     will immeditaley occur.
c
c     if(maswrk) write(iw,*) '-- Intranode GDDI with',NGROUPS,' groups.'
c     write(6,9000) nproc,me,MYGROUP,NGROUPS,goparr,maswrk,iscope
c9000 format(1x,"GDDI: nproc=",I3," me=",I3,' mygroup=',I3,' ngroups=',
c    *       I3,' goparr=',L1,' maswrk=',L1,' at scope=',I1)
c     if(mode.eq.0) CALL DDI_SYNC(4567)
      if(mode.eq.0) then
c       CALL DDI_SYNC(4567)
c       It seems that a sync is not needed...
      endif
      RETURN
c9000 format(1x,"GDDI: nproc=",I3," me=",I3,' mygroup=',I3,' ngroups=',
c    *       I3,' goparr=',L1,' maswrk=',L1,' at scope=',I1)
      END
