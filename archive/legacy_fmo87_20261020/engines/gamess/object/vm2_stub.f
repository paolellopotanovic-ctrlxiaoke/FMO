!*MODULE VM2_STUB DECK VC_CompPackage
C>  @brief Stub file for VeraChem interface
C>
C>  @author Peng Xu, Tosaporn Sattasathuchana
      SUBROUTINE VC_CompPackage()

      implicit none

      RETURN
      END

C--------------------------------------------------------------------
!*MODULE VM2_STUB DECK vchem_set_gms_control 
C>  @brief Stub file for VeraChem interface
C>
C>  @author Peng Xu, Tosaporn Sattasathuchana
      SUBROUTINE vchem_set_gms_control()

      implicit none

      RETURN
      END

C--------------------------------------------------------------------
!*MODULE VM2_STUB DECK vchem_set_gms_basis_set
C>  @brief Stub file for VeraChem interface
C>
C>  @author Peng Xu, Tosaporn Sattasathuchana
      SUBROUTINE vchem_set_gms_basis_set()

      implicit none


      RETURN
      END

C--------------------------------------------------------------------
!*MODULE VM2_STUB DECK vchem_geomopt_init
C>  @brief Stub file for VeraChem interface
C>
C>  @author Peng Xu, Tosaporn Sattasathuchana
      SUBROUTINE vchem_geomopt_init(ncval,ncoord,method,out,
     &                              goteg,nprt,npun,itbmat,
     &                              stpt,ststep,projct,snzvar,
     &                              sozvar,max,nprths,purify,
     &                              hess)

      implicit none

      INTEGER max

      DOUBLE PRECISION method
      DOUBLE PRECISION ststep
      DOUBLE PRECISION snzvar(max)
      DOUBLE PRECISION sozvar(max)
      DOUBLE PRECISION hess

      INTEGER ncval
      INTEGER ncoord
      INTEGER nprt
      INTEGER npun
      INTEGER itbmat
      INTEGER nprths

      LOGICAL out,goteg
      LOGICAL purify,stpt,projct

      RETURN
      END

C--------------------------------------------------------------------
!*MODULE VM2_STUB DECK vchem_set_gms_dft
C>  @brief Stub file for VeraChem interface
C>
C>  @author Peng Xu, Tosaporn Sattasathuchana
      SUBROUTINE vchem_set_gms_dft()  

      implicit none

      RETURN
      END

C--------------------------------------------------------------------
!*MODULE VM2_STUB DECK vchem_set_vscf_control
C>  @brief Stub file for VeraChem interface
C>
C>  @author Peng Xu, Tosaporn Sattasathuchana
      SUBROUTINE vchem_set_vscf_control(imode,mode)

      implicit none

      INTEGER imode(*)
      INTEGER mode

      RETURN
      END
    
