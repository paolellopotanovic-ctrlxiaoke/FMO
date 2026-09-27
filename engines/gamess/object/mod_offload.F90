! MODULE MOD_OFFLOAD *DECK MOD_OFFLOAD
!>
!> @BRIEF   CONTAINS THE VARIABLES FOR THE OFFLOAD INPUT GROUP
!>
!> @AUTHOR  SAROM LEANG
!>
!> @DATE    JULY 9, 2023
!>
!> @PARAM RHF LOGICAL VARIABLE TO INDICATE TARGET OFFLOAD OF S,P INTEGRALS DURING RHF STEP
!> @PARAM MAKEFP LOGICAL VARIABLE TO INDICATE TARGET OFFLOAD OF S,P INTEGRALS DURING MAKEFP STEP
!> @PARAM CPHF LOGICAL VARIABLE TO INDICATE TARGET OFFLOAD  DURING CPHF STEP (FUTURE)

      MODULE MOD_OFFLOAD
      IMPLICIT NONE
      PUBLIC

      LOGICAL :: OFFLOAD_RHF = .FALSE.
      LOGICAL :: OFFLOAD_MAKEFP = .FALSE.
      LOGICAL :: OFFLOAD_CPHF = .FALSE.

      SAVE

      CONTAINS

         SUBROUTINE OFFLOADINP
            USE COMM_PAR, ONLY: MASWRK
            USE COMM_IOFILE, ONLY: IR,IW

            IMPLICIT NONE

            INTEGER, PARAMETER :: NNAM=3
            DOUBLE PRECISION   :: QNAM(NNAM)
            INTEGER            :: KQNAM(NNAM)
            INTEGER            :: JRET
            DOUBLE PRECISION   :: OFFLOAD
            LOGICAL            :: RHF,MAKEFP,CPHF

            DATA OFFLOAD/8HOFFLOAD /
            DATA QNAM   /8HRHF     ,8HMAKEFP  ,8HCPHF    /
            DATA KQNAM  /0,0,0/

            JRET=0

            RHF=.false.
            CPHF=.false.
            MAKEFP=.false.

            CALL NAMEIO(IR,JRET,OFFLOAD,NNAM,QNAM,KQNAM,                &
            RHF,MAKEFP,CPHF,0,                                          &
             0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,         &
             0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,         &
             0,0,0,0,0,    0,0,0,0,0,   0,0,0,0,0,   0,0,0,0,0,         &
             0,0,0,0,0,    0)

            IF(JRET.EQ.2) THEN
               IF (MASWRK) WRITE (IW,*) 'ERROR READING $OFFLOAD GROUP'
               CALL ABRT
            END IF

            OFFLOAD_RHF=RHF
            OFFLOAD_MAKEFP=MAKEFP
            OFFLOAD_CPHF=CPHF
            if(OFFLOAD_MAKEFP) OFFLOAD_CPHF=.true.

            IF(MASWRK) WRITE(IW,9000) OFFLOAD_RHF,OFFLOAD_CPHF,         &
                       OFFLOAD_MAKEFP

9000        FORMAT(/5X,22(1H-)/,                                        &
            5X,'TARGET OFFLOAD OPTIONS'/5X,22(1H-)/,                    &
            5X,'RHF',4X,'=',L8/,                                        &
            5X,'CPHF',4X,'=',L8/,                                       &
            5X,'MAKEFP',1X,'=',L8)

         END SUBROUTINE OFFLOADINP

      END MODULE MOD_OFFLOAD
