C>
C>     @author Simon P. Webb, Michael W. Schmidt
C>
C>     @brief A*X
C>
C>     @details Forms the product of a matrix -A- and
C>              a vector -X-
C>
C>
C*MODULE MP2GRD-DEFENSIVE  *DECK MPSPAX
      SUBROUTINE MPSPAX(A,LDA,X,AX,N)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D+00)
C
      COMMON /MP2PTR/ NFT,J00,J01,J02,IDUMMY(4)
      COMMON /FMCOM / XX(1)
C
      DIMENSION IX(1)
      EQUIVALENCE (XX(1),IX(1))
      DIMENSION A(*),X(*),AX(*)
C
      II(I) = IX(J00-1+I)
      JJ(I) = IX(J01-1+I)
C
      DO 10 I=1,N
         AX(I) = ZERO
   10 CONTINUE
C
      IJ = 0
      DO 40 I=1,N
         NUM = II(I)
         DO 20 J=1,NUM
            AX(I) = AX(I) + A(IJ+J)*X(JJ(IJ+J))
   20    CONTINUE
         DO 30 J=1,NUM-1
            AX(JJ(IJ+J)) = AX(JJ(IJ+J)) + A(IJ+J)*X(I)
   30    CONTINUE
         IJ = IJ+NUM
   40 CONTINUE
C
      DO 50 I=1,N
         AX(I)=AX(I)*XX(J02-1+I)
   50 CONTINUE
C
      RETURN
      END
C
