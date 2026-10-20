! MODULE COMM_VVOPAR  *DECK VVOS
!>
!> @brief   Contains the variables from the VVOPAR common block
!>
!> @author  George Schoendorff
!>
!> @date    July 13, 2022
!>
!> @param IVVOS  Option to select VVOs (=0) or skip VVOS (=1)
!> @param NVVOS  Number of valence virtual orbitals
!> @param IVVTYP Option to select VVOs and QUAOs based on the
!>                nonrelativistic (=0) or relativistic (=1) AAMBS
!> @param BNDDEN Density matrix
!>
      module comm_vvopar
      use mx_limits, only: mxatm
      use prec, only: dp
      implicit none
      public

      integer :: ivvos,nvvos,ivvtyp,mxgmbs

      real(kind=dp) :: bndden

      save
      end module comm_vvopar
!
! MODULE COMM_ORNTMO  *DECK LOCAL
!>
!> @brief   Contains the variables from the ORNTMO common block
!>
!> @author  George Schoendorff
!>
!> @date    July 23, 2022
!>
!> @param ISVMOR   Indexing array
!> @param ISVDOP   Indexing array
!> @param L0DEP0   The current number of final spherical MOs
!> @param EXTLOC   Options to localize the external orbitals
!> @parem EREFATM  Atomic reference energy
!> @param BOTOL    Array of bond order tolerances
!> @param ORIENT   A flag to perform the orientation of the localized
!>                  orbitals
!> @param ORIKIN   A flag that switches print tolerances for kinetic
!>                  energy bond orders
!> @param ORMFUL   A flag to select orbitals between ORMAS groups
!> @param ORMRESET A flag to reset the ormas groups
!> @param RNNTSCF  A flag to select NT SCF tolerances
!> @param SVDCOPT  A flag to select the SVD of the core orbitals
!> @param ENGATM   A flag to select the coupled NT
!> @param PADSVD   A flag to rotate localized core orbitals into the
!>                  valence space for when a valence orbital on one
!>                  atom has a lower energy than a core orbital on 
!>                  another atom
!>
      module comm_orntmo
      use mx_limits, only: mxatm
      use prec, only: dp
      implicit none
      public

      integer :: isvmor(mxatm),isvdop,l0dep0

      real(kind=dp) :: extloc,erefatm,botol(6)

      logical :: orient,orikin,ormful,ormreset,rnntscf,svdcopt
      logical :: engatm,padsvd

      save
      end module comm_orntmo
!
! MODULE AAMBS_LIMITS *DECK AAMBS
!>
!> @brief   Contains dimensions of the exponent and coefficient arrays
!>          for the relativistic AAMBS
!>
!> @author  George Schoendorff
!>
!> @date    July 13, 2022
!>
!> @details Each parameter is named as follows:
!>          First part:  NE = number of exponents
!>          Second part = Atomic symbol
!>          Third part  = Angular momentum if exponents or
!>                              orbital if coefficients
!>          Fourth part = R if relativistic (so far only the
!>                        relativistic option is implemented
!>                        except for H and He) 
!>
!> Non-relativistic paramters:
!>
!> @param NEHS     Number of s exponents for H
!> @param NEHES    Number of s exponents for He
!> @param NELIS    Number of s exponents for Li
!> @param NELIP    Number of p exponents for Li
!> @param NEBES    Number of s exponents for Be
!> @param NEBEP    Number of p exponents for Be
!> @param NEBS     Number of s exponents for B
!> @param NEBP     Number of p exponents for B
!> @param NECS     Number of s exponents for C
!> @param NECP     Number of p exponents for C
!> @param NENS     Number of s exponents for N
!> @param NENP     Number of p exponents for N
!> @param NEOS     Number of s exponents for O
!> @param NEOP     Number of p exponents for O
!> @param NEFS     Number of s exponents for F
!> @param NEFP     Number of P exponents for F
!> @param NENES    Number of s exponents for Ne
!> @param NENEP    Number of p exponents for Ne
!> @param NENAS    Number of s exponents for Na
!> @param NENAP    Number of p exponents for Na
!> @param NEMGS    Number of s exponents for Mg
!> @param NEMGP    Number of p exponents for Mg
!> @param NEALS    Number of s exponents for Al
!> @param NEALP    Number of p exponents for Al
!> @param NESIS    Number of s exponents for Si
!> @param NESIP    Number of p exponents for Si
!> @param NEPS     Number of s exponents for P
!> @param NEPP     Number of p exponents for P
!> @param NESS     Number of s exponents for S
!> @param NESP     Number of p exponents for S
!> @param NECLS    Number of s exponents for Cl
!> @param NECLP    Number of p exponents for Cl
!> @param NEARS    Number of s exponents for Ar
!> @param NEARP    Number of p exponents for Ar
!> @param NEKS     Number of s exponents for K
!> @param NEKP     Number of p exponents for K
!> @param NEKD     Number of d exponents for K
!> @param NECAS    Number of s exponents for Ca
!> @param NECAP    Number of p exponents for Ca
!> @param NECAD    Number of d exponents for Ca
!> @param NESCS    Number of s exponents for Sc
!> @param NESCP    Number of p exponents for Sc
!> @param NESCD    Number of d exponents for Sc
!> @param NETIS    Number of s exponents for Ti
!> @param NETIP    Number of p exponents for Ti
!> @param NETID    Number of d exponents for Ti
!> @param NEVS     Number of s exponents for V
!> @param NEVP     Number of p exponents for V
!> @param NEVD     Number of d exponents for V
!> @param NECRS    Number of s exponents for Cr
!> @param NECRP    Number of p exponents for Cr
!> @param NECRD    Number of d exponents for Cr
!> @param NEMNS    Number of s exponents for Mn
!> @param NEMNP    Number of p exponents for Mn
!> @param NEMND    Number of d exponents for Mn
!> @param NEFES    Number of s exponents for Fe
!> @param NEFEP    Number of p exponents for Fe
!> @param NEFED    Number of d exponents for Fe
!> @param NECOS    Number of s exponents for Co
!> @param NECOP    Number of p exponents for Co
!> @param NECOD    Number of d exponents for Co
!> @param NENIS    Number of s exponents for Ni
!> @param NENIP    Number of p exponents for Ni
!> @param NENID    Number of d exponents for Ni
!> @param NECUS    Number of s exponents for Cu
!> @param NECUP    Number of p exponents for Cu
!> @param NECUD    Number of d exponents for Cu
!> @param NEZNS    Number of s exponents for Zn
!> @param NEZNP    Number of p exponents for Zn
!> @param NEZND    Number of d exponents for Zn
!> @param NEGAS    Number of s exponents for Ga
!> @param NEGAP    Number of p exponents for Ga
!> @param NEGAD    Number of d exponents for Ga
!> @param NEGES    Number of s exponents for Ge
!> @param NEGEP    Number of p exponents for Ge
!> @param NEGED    Number of d exponents for Ge
!> @param NEASS    Number of s exponents for As
!> @param NEASP    Number of p exponents for As
!> @param NEASD    Number of d exponents for As
!> @param NESES    Number of s exponents for Se
!> @param NESEP    Number of p exponents for Se
!> @param NESED    Number of d exponents for Se
!> @param NEBRS    Number of s exponents for Br
!> @param NEBRP    Number of p exponents for Br
!> @param NEBRD    Number of d exponents for Br
!> @param NEKRS    Number of s exponents for Kr
!> @param NEKRP    Number of p exponents for Kr
!> @param NEKRD    Number of d exponents for Kr
!> @param NERBS    Number of s exponents for Rb
!> @param NERBP    Number of p exponents for Rb
!> @param NERBD    Number of d exponents for Rb
!> @param NESRS    Number of s exponents for Sr
!> @param NESRP    Number of p exponents for Sr
!> @param NESRD    Number of d exponents for Sr
!> @param NEYS     Number of s exponents for Y
!> @param NEYP     Number of p exponents for Y
!> @param NEYD     Number of d exponents for Y
!> @param NEZRS    Number of s exponents for Zr
!> @param NEZRP    Number of p exponents for Zr
!> @param NEZRD    Number of d exponents for Zr
!> @param NENBS    Number of s exponents for Nb
!> @param NENBP    Number of p exponents for Nb
!> @param NENBD    Number of d exponents for Nb
!> @param NEMOS    Number of s exponents for Mo
!> @param NEMOP    Number of p exponents for Mo
!> @param NEMOD    Number of d exponents for Mo
!> @param NETCS    Number of s exponents for Tc
!> @param NETCP    Number of p exponents for Tc
!> @param NETCD    Number of d exponents for Tc
!> @param NERUS    Number of s exponents for Ru
!> @param NERUP    Number of p exponents for Ru
!> @param NERUD    Number of d exponents for Ru
!> @param NERHS    Number of s exponents for Rh
!> @param NERHP    Number of p exponents for Rh
!> @param NERHD    Number of d exponents for Rh
!> @param NEPDS    Number of s exponents for Pd
!> @param NEPDP    Number of p exponents for Pd
!> @param NEPDD    Number of d exponents for Pd
!> @param NEAGS    Number of s exponents for Ag
!> @param NEAGP    Number of p exponents for Ag
!> @param NEAGD    Number of d exponents for Ag
!> @param NECDS    Number of s exponents for Cd
!> @param NECDP    Number of p exponents for Cd
!> @param NECDD    Number of d exponents for Cd
!> @param NEINS    Number of s exponents for In
!> @param NEINP    Number of p exponents for In
!> @param NEIND    Number of d exponents for In
!> @param NESNS    Number of s exponents for Sn
!> @param NESNP    Number of p exponents for Sn
!> @param NESND    Number of d exponents for Sn
!> @param NESBS    Number of s exponents for Sb
!> @param NESBP    Number of p exponents for Sb
!> @param NESBD    Number of d exponents for Sb
!> @param NETES    Number of s exponents for Te
!> @param NETEP    Number of p exponents for Te
!> @param NETED    Number of d exponents for Te
!> @param NEIS     Number of s exponents for I
!> @param NEIP     Number of p exponents for I
!> @param NEID     Number of d exponents for I
!> @param NEXES    Number of s exponents for Xe
!> @param NEXEP    Number of p exponents for Xe
!> @param NEXED    Number of d exponents for Xe
!>
!> Relativistic paramters:
!>
!> @param NELISR   Number of s exponents for Li
!> @param NELIPR   Number of p exponents for Li
!> @param NEBESR   Number of s exponents for Be
!> @param NEBEPR   Number of p exponents for Be
!> @param NEBSR    Number of s exponents for B
!> @param NEBPR    Number of p exponents for B
!> @param NECSR    Number of s exponents for C
!> @param NECPR    Number of p exponents for C
!> @param NENSR    Number of s exponents for N
!> @param NENPR    Number of p exponents for N
!> @param NEOSR    Number of s exponents for O
!> @param NEOPR    Number of p exponents for O
!> @param NEFSR    Number of s exponents for F
!> @param NEFPR    Number of P exponents for F
!> @param NENESR   Number of s exponents for Ne
!> @param NENEPR   Number of p exponents for Ne
!> @param NENASR   Number of s exponents for Na
!> @param NENAPR   Number of p exponents for Na
!> @param NEMGSR   Number of s exponents for Mg
!> @param NEMGPR   Number of p exponents for Mg
!> @param NEALSR   Number of s exponents for Al
!> @param NEALPR   Number of p exponents for Al
!> @param NESISR   Number of s exponents for Si
!> @param NESIPR   Number of p exponents for Si
!> @param NEPSR    Number of s exponents for P
!> @param NEPPR    Number of p exponents for P
!> @param NESSR    Number of s exponents for S
!> @param NESPR    Number of p exponents for S
!> @param NECLSR   Number of s exponents for Cl
!> @param NECLPR   Number of p exponents for Cl
!> @param NEARSR   Number of s exponents for Ar
!> @param NEARPR   Number of p exponents for Ar
!> @param NEKSR    Number of s exponents for K
!> @param NEKPR    Number of p exponents for K
!> @param NEKDR    Number of d exponents for K
!> @param NECASR   Number of s exponents for Ca
!> @param NECAPR   Number of p exponents for Ca
!> @param NECADR   Number of d exponents for Ca
!> @param NESCSR   Number of s exponents for Sc
!> @param NESCPR   Number of p exponents for Sc
!> @param NESCDR   Number of d exponents for Sc
!> @param NETISR   Number of s exponents for Ti
!> @param NETIPR   Number of p exponents for Ti
!> @param NETIDR   Number of d exponents for Ti
!> @param NEVSR    Number of s exponents for V
!> @param NEVPR    Number of p exponents for V
!> @param NEVDR    Number of d exponents for V
!> @param NECRSR   Number of s exponents for Cr
!> @param NECRPR   Number of p exponents for Cr
!> @param NECRDR   Number of d exponents for Cr
!> @param NEMNSR   Number of s exponents for Mn
!> @param NEMNPR   Number of p exponents for Mn
!> @param NEMNDR   Number of d exponents for Mn
!> @param NEFESR   Number of s exponents for Fe
!> @param NEFEPR   Number of p exponents for Fe
!> @param NEFEDR   Number of d exponents for Fe
!> @param NECOSR   Number of s exponents for Co
!> @param NECOPR   Number of p exponents for Co
!> @param NECODR   Number of d exponents for Co
!> @param NENISR   Number of s exponents for Ni
!> @param NENIPR   Number of p exponents for Ni
!> @param NENIDR   Number of d exponents for Ni
!> @param NECUSR   Number of s exponents for Cu
!> @param NECUPR   Number of p exponents for Cu
!> @param NECUDR   Number of d exponents for Cu
!> @param NEZNSR   Number of s exponents for Zn
!> @param NEZNPR   Number of p exponents for Zn
!> @param NEZNDR   Number of d exponents for Zn
!> @param NEGASR   Number of s exponents for Ga
!> @param NEGAPR   Number of p exponents for Ga
!> @param NEGADR   Number of d exponents for Ga
!> @param NEGESR   Number of s exponents for Ge
!> @param NEGEPR   Number of p exponents for Ge
!> @param NEGEDR   Number of d exponents for Ge
!> @param NEASSR   Number of s exponents for As
!> @param NEASPR   Number of p exponents for As
!> @param NEASDR   Number of d exponents for As
!> @param NESESR   Number of s exponents for Se
!> @param NESEPR   Number of p exponents for Se
!> @param NESEDR   Number of d exponents for Se
!> @param NEBRSR   Number of s exponents for Br
!> @param NEBRPR   Number of p exponents for Br
!> @param NEBRDR   Number of d exponents for Br
!> @param NEKRSR   Number of s exponents for Kr
!> @param NEKRPR   Number of p exponents for Kr
!> @param NEKRDR   Number of d exponents for Kr
!> @param NERBSR   Number of s exponents for Rb
!> @param NERBPR   Number of p exponents for Rb
!> @param NERBDR   Number of d exponents for Rb
!> @param NESRSR   Number of s exponents for Sr
!> @param NESRPR   Number of p exponents for Sr
!> @param NESRDR   Number of d exponents for Sr
!> @param NEYSR    Number of s exponents for Y
!> @param NEYPR    Number of p exponents for Y
!> @param NEYDR    Number of d exponents for Y
!> @param NEZRSR   Number of s exponents for Zr
!> @param NEZRPR   Number of p exponents for Zr
!> @param NEZRDR   Number of d exponents for Zr
!> @param NENBSR   Number of s exponents for Nb
!> @param NENBPR   Number of p exponents for Nb
!> @param NENBDR   Number of d exponents for Nb
!> @param NEMOSR   Number of s exponents for Mo
!> @param NEMOPR   Number of p exponents for Mo
!> @param NEMODR   Number of d exponents for Mo
!> @param NETCSR   Number of s exponents for Tc
!> @param NETCPR   Number of p exponents for Tc
!> @param NETCDR   Number of d exponents for Tc
!> @param NERUSR   Number of s exponents for Ru
!> @param NERUPR   Number of p exponents for Ru
!> @param NERUDR   Number of d exponents for Ru
!> @param NERHSR   Number of s exponents for Rh
!> @param NERHPR   Number of p exponents for Rh
!> @param NERHDR   Number of d exponents for Rh
!> @param NEPDSR   Number of s exponents for Pd
!> @param NEPDPR   Number of p exponents for Pd
!> @param NEPDDR   Number of d exponents for Pd
!> @param NEAGSR   Number of s exponents for Ag
!> @param NEAGPR   Number of p exponents for Ag
!> @param NEAGDR   Number of d exponents for Ag
!> @param NECDSR   Number of s exponents for Cd
!> @param NECDPR   Number of p exponents for Cd
!> @param NECDDR   Number of d exponents for Cd
!> @param NEINSR   Number of s exponents for In
!> @param NEINPR   Number of p exponents for In
!> @param NEINDR   Number of d exponents for In
!> @param NESNSR   Number of s exponents for Sn
!> @param NESNPR   Number of p exponents for Sn
!> @param NESNDR   Number of d exponents for Sn
!> @param NESBSR   Number of s exponents for Sb
!> @param NESBPR   Number of p exponents for Sb
!> @param NESBDR   Number of d exponents for Sb
!> @param NETESR   Number of s exponents for Te
!> @param NETEPR   Number of p exponents for Te
!> @param NETEDR   Number of d exponents for Te
!> @param NEISR    Number of s exponents for I
!> @param NEIPR    Number of p exponents for I
!> @param NEIDR    Number of d exponents for I
!> @param NEXESR   Number of s exponents for Xe
!> @param NEXEPR   Number of p exponents for Xe
!> @param NEXEDR   Number of d exponents for Xe
!> @param NECSSR   Number of s exponents for Cs
!> @param NECSPR   Number of p exponents for Cs
!> @param NECSDR   Number of d exponents for Cs
!> @param NEBASR   Number of s exponents for Ba
!> @param NEBAPR   Number of p exponents for Ba
!> @param NEBADR   Number of d exponents for Ba
!> @param NELASR   Number of d exponents for La
!> @param NELAPR   Number of p exponents for La
!> @param NELADR   Number of d exponents for La
!> @param NELAFR   Number of f exponents for La
!> @param NECESR   Number of s exponents for Ce
!> @param NECEPR   Number of p exponents for Ce
!> @param NECEDR   Number of d exponents for Ce
!> @param NECEFR   Number of f exponents for Ce
!> @param NEPRSR   Number of s exponents for Pr
!> @param NEPRPR   Number of p exponents for Pr
!> @param NEPRDR   Number of d exponents for Pr
!> @param NEPRFR   Number of f exponents for Pr
!> @param NENDSR   Number of s exponents for Nd
!> @param NENDPR   Number of p exponents for Nd
!> @param NENDDR   Number of d exponents for Nd
!> @param NENDFR   Number of f exponents for Nd
!> @param NEPMSR   Number of s exponents for Pm
!> @param NEPMPR   Number of p exponents for Pm
!> @param NEPMDR   Number of d exponents for Pm
!> @param NEPMFR   Number of f exponents for Pm
!> @param NESMSR   Number of s exponents for Sm
!> @param NESMPR   Number of p exponents for Sm
!> @param NESMDR   Number of d exponents for Sm
!> @param NESMFR   Number of f exponents for Sm
!> @param NEEUSR   Number of s exponents for Eu
!> @param NEEUPR   Number of p exponents for Eu
!> @param NEEUDR   Number of d exponents for Eu
!> @param NEEUFR   Number of f exponents for Eu
!> @param NEGDSR   Number of s exponents for Gd
!> @param NEGDPR   Number of p exponents for Gd
!> @param NEGDDR   Number of d exponents for Gd
!> @param NEGDFR   Number of f exponents for Gd
!> @param NETBSR   Number of s exponents for Tb
!> @param NETBPR   Number of p exponents for Tb
!> @param NETBDR   Number of d exponents for Tb
!> @param NETBFR   Number of f exponents for Tb
!> @param NEDYSR   Number of s exponents for Dy
!> @param NEDYPR   Number of p exponents for Dy
!> @param NEDYDR   Number of d exponents for Dy
!> @param NEDYFR   Number of f exponents for Dy
!> @param NEHOSR   Number of s exponents for Ho
!> @param NEHOPR   Number of p exponents for Ho
!> @param NEHODR   Number of d exponents for Ho
!> @param NEHOFR   Number of f exponents for Ho
!> @param NEERSR   Number of s exponents for Er
!> @param NEERPR   Number of p exponents for Er
!> @param NEERDR   Number of d exponents for Er
!> @param NEERFR   Number of f exponents for Er
!> @param NETMSR   Number of s exponents for Tm
!> @param NETMPR   Number of p exponents for Tm
!> @param NETMDR   Number of d exponents for Tm
!> @param NETMFR   Number of f exponents for Tm
!> @param NEYBSR   Number of s exponents for Yb
!> @param NEYBPR   Number of p exponents for Yb
!> @param NEYBDR   Number of d exponents for Yb
!> @param NEYBFR   Number of f exponents for Yb
!> @param NELUSR   Number of s exponents for Lu
!> @param NELUPR   Number of p exponents for Lu
!> @param NELUDR   Number of d exponents for Lu
!> @param NELUFR   Number of f exponents for Lu
!> @param NEHFSR   Number of s exponents for Hf
!> @param NEHFPR   Number of p exponents for Hf
!> @param NEHFDR   Number of d exponents for Hf
!> @param NEHFFR   Number of f exponents for Hf
!> @param NETASR   Number of s exponents for Ta
!> @param NETAPR   Number of p exponents for Ta
!> @param NETADR   Number of d exponents for Ta
!> @param NETAFR   Number of f exponents for Ta
!> @param NEWSR    Number of s exponents for W
!> @param NEWPR    Number of p exponents for W
!> @param NEWDR    Number of d exponents for W
!> @param NEWFR    Number of f exponents for W
!> @param NERESR   Number of s exponents for Re
!> @param NEREPR   Number of p exponents for Re
!> @param NEREDR   Number of d exponents for Re
!> @param NEREFR   Number of f exponents for Re
!> @param NEOSSR   Number of s exponents for Os
!> @param NEOSPR   Number of p exponents for Os
!> @param NEOSDR   Number of d exponents for Os
!> @param NEOSFR   Number of f exponents for Os
!> @param NEIRSR   Number of s exponents for Ir
!> @param NEIRPR   Number of p exponents for Ir
!> @param NEIRDR   Number of d exponents for Ir
!> @param NEIRFR   Number of f exponents for Ir
!> @param NEPTSR   Number of s exponents for Pt
!> @param NEPTPR   Number of p exponents for Pt
!> @param NEPTDR   Number of d exponents for Pt
!> @param NEPTFR   Number of f exponents for Pt
!> @param NEAUSR   Number of s exponents for Au
!> @param NEAUPR   Number of p exponents for Au
!> @param NEAUDR   Number of d exponents for Au
!> @param NEAUFR   Number of f exponents for Au
!> @param NEHGSR   Number of s exponents for Hg
!> @param NEHGPR   Number of p exponents for Hg
!> @param NEHGDR   Number of d exponents for Hg
!> @param NEHGFR   Number of f exponents for Hg
!> @param NETLSR   Number of s exponents for Tl
!> @param NETLPR   Number of p exponents for Tl
!> @param NETLDR   Number of d exponents for Tl
!> @param NETLFR   Number of f exponents for Tl
!> @param NEPBSR   Number of s exponents for Pb
!> @param NEPBPR   Number of p exponents for Pb
!> @param NEPBDR   Number of d exponents for Pb
!> @param NEPBFR   Number of f exponents for Pb
!> @param NEBISR   Number of s exponents for Bi
!> @param NEBIPR   Number of p exponents for Bi
!> @param NEBIDR   Number of d exponents for Bi
!> @param NEBIFR   Number of f exponents for Bi
!> @param NEPOSR   Number of s exponents for Po
!> @param NEPOPR   Number of p exponents for Po
!> @param NEPODR   Number of d exponents for Po
!> @param NEPOFR   Number of f exponents for Po
!> @param NEATSR   Number of s exponents for At
!> @param NEATPR   Number of p exponents for At
!> @param NEATDR   Number of d exponents for At
!> @param NEATFR   Number of f exponents for At
!> @param NERNSR   Number of s exponents for Rn
!> @param NERNPR   Number of p exponents for Rn
!> @param NERNDR   Number of d exponents for Rn
!> @param NERNFR   Number of f exponents for Rn
!>
      module aambs_limits
      implicit none
      public

      INTEGER, PARAMETER :: & 
               NEHS   = 8,  &
               NEHES  = 8,  &
               NELIS  = 14, &
               NELIP  = 1, &
               NEBES  = 14, &
               NEBEP  = 1, &
               NEBS  = 14, &
               NEBP  = 7, &
               NECS  = 14, &
               NECP  = 7, &
               NENS  = 14, &
               NENP  = 7, &
               NEOS  = 14, &
               NEOP  = 7, &
               NEFS  = 14, &
               NEFP  = 7, &
               NENES  = 14, &
               NENEP  = 7, &
               NENAS  = 18, &
               NENAP  = 9, &
               NEMGS  = 18, &
               NEMGP  = 9, &
               NEALS  = 18, &
               NEALP  = 12, &
               NESIS  = 18, &
               NESIP  = 12, &
               NEPS  = 18, &
               NEPP  = 12, &
               NESS  = 18, &
               NESP  = 12, &
               NECLS  = 18, &
               NECLP  = 12, &
               NEARS  = 18, &
               NEARP  = 12, &
               NEKS  = 26, &
               NEKP  = 16, &
!               NEKD  = 15, &
               NECAS  = 26, &
               NECAP  = 16, &
!               NECAD  = 15, &
               NESCS  = 26, &
               NESCP  = 17, &
               NESCD  = 13, &
               NETIS  = 26, &
               NETIP  = 17, &
               NETID  = 13, &
               NEVS  = 26, &
               NEVP  = 17, &
               NEVD  = 13, &
               NECRS  = 26, &
               NECRP  = 17, &
               NECRD  = 13, &
               NEMNS  = 26, &
               NEMNP  = 17, &
               NEMND  = 13, &
               NEFES  = 26, &
               NEFEP  = 17, &
               NEFED  = 13, &
               NECOS  = 26, &
               NECOP  = 17, &
               NECOD  = 13, &
               NENIS  = 26, &
               NENIP  = 17, &
               NENID  = 13, &
               NECUS  = 26, &
               NECUP  = 17, &
               NECUD  = 14, &
               NEZNS  = 26, &
               NEZNP  = 17, &
               NEZND  = 14, &
               NEGAS  = 26, &
               NEGAP  = 20, &
               NEGAD  = 14, &
               NEGES  = 26, &
               NEGEP  = 20, &
               NEGED  = 14, &
               NEASS  = 26, &
               NEASP  = 20, &
               NEASD  = 14, &
               NESES  = 26, &
               NESEP  = 20, &
               NESED  = 14, &
               NEBRS  = 26, &
               NEBRP  = 20, &
               NEBRD  = 14, &
               NEKRS  = 26, &
               NEKRP  = 20, &
               NEKRD  = 14, &
               NERBS  = 28, &
               NERBP  = 20, &
               NERBD  = 14, &
               NESRS  = 28, &
               NESRP  = 20, &
               NESRD  = 14, &
               NEYS  = 27, &
               NEYP  = 20, &
               NEYD  = 17, &
               NEZRS  = 27, &
               NEZRP  = 20, &
               NEZRD  = 17, &
               NENBS  = 27, &
               NENBP  = 20, &
               NENBD  = 17, &
               NEMOS  = 27, &
               NEMOP  = 20, &
               NEMOD  = 17, &
               NETCS  = 27, &
               NETCP  = 20, &
               NETCD  = 17, &
               NERUS  = 28, &
               NERUP  = 20, &
               NERUD  = 17, &
               NERHS  = 28, &
               NERHP  = 20, &
               NERHD  = 17, &
               NEPDS  = 28, &
               NEPDP  = 20, &
               NEPDD  = 17, &
               NEAGS  = 28, &
               NEAGP  = 20, &
               NEAGD  = 17, &
               NECDS  = 28, &
               NECDP  = 20, &
               NECDD  = 17, &
               NEINS  = 28, &
               NEINP  = 23, &
               NEIND  = 17, &
               NESNS  = 28, &
               NESNP  = 23, &
               NESND  = 17, &
               NESBS  = 28, &
               NESBP  = 23, &
               NESBD  = 17, &
               NETES  = 28, &
               NETEP  = 23, &
               NETED  = 17, &
               NEIS   = 28, &
               NEIP   = 23, &
               NEID   = 17, &
               NEXES  = 28, &
               NEXEP  = 23, &
               NEXED  = 17, &
!
               NELISR = 20, & 
               NELIPR = 13, & 
               NEBESR = 20, & 
               NEBEPR = 13, & 
               NEBSR = 20, & 
               NEBPR = 13, & 
               NECSR = 20, & 
               NECPR = 13, & 
               NENSR = 20, & 
               NENPR = 13, & 
               NEOSR = 20, & 
               NEOPR = 13, & 
               NEFSR = 20, & 
               NEFPR = 13, & 
               NENESR = 20, & 
               NENEPR = 13, & 
               NENASR = 24, & 
               NENAPR = 17, & 
               NEMGSR = 24, & 
               NEMGPR = 17, & 
               NEALSR = 23, & 
               NEALPR = 16, & 
               NESISR = 23, & 
               NESIPR = 16, & 
               NEPSR = 23, & 
               NEPPR = 16, & 
               NESSR = 23, & 
               NESPR = 16, & 
               NECLSR = 23, & 
               NECLPR = 16, & 
               NEARSR = 23, & 
               NEARPR = 16, & 
               NEKSR = 26, & 
               NEKPR = 20, & 
               NEKDR = 15, & 
               NECASR = 26, & 
               NECAPR = 20, & 
               NECADR = 15, & 
               NESCSR = 26, & 
               NESCPR = 20, & 
               NESCDR = 13, & 
               NETISR = 26, & 
               NETIPR = 20, & 
               NETIDR = 13, & 
               NEVSR = 26, & 
               NEVPR = 20, & 
               NEVDR = 13, & 
               NECRSR = 26, & 
               NECRPR = 20, & 
               NECRDR = 13, & 
               NEMNSR = 26, & 
               NEMNPR = 20, & 
               NEMNDR = 13, & 
               NEFESR = 26, & 
               NEFEPR = 20, & 
               NEFEDR = 13, & 
               NECOSR = 26, & 
               NECOPR = 20, & 
               NECODR = 13, & 
               NENISR = 26, & 
               NENIPR = 20, & 
               NENIDR = 13, & 
               NECUSR = 26, & 
               NECUPR = 21, & 
               NECUDR = 14, & 
               NEZNSR = 26, & 
               NEZNPR = 20, & 
               NEZNDR = 14, & 
               NEGASR = 26, & 
               NEGAPR = 20, & 
               NEGADR = 14, & 
               NEGESR = 26, & 
               NEGEPR = 20, & 
               NEGEDR = 14, & 
               NEASSR = 26, & 
               NEASPR = 20, & 
               NEASDR = 14, & 
               NESESR = 26, & 
               NESEPR = 20, & 
               NESEDR = 14, & 
               NEBRSR = 26, & 
               NEBRPR = 20, & 
               NEBRDR = 14, & 
               NEKRSR = 26, & 
               NEKRPR = 20, & 
               NEKRDR = 14, & 
               NERBSR = 30, & 
               NERBPR = 26, & 
               NERBDR = 21, & 
               NESRSR = 29, & 
               NESRPR = 25, & 
               NESRDR = 21, & 
               NEYSR = 29, & 
               NEYPR = 25, & 
               NEYDR = 21, & 
               NEZRSR = 29, & 
               NEZRPR = 25, & 
               NEZRDR = 21, & 
               NENBSR = 29, & 
               NENBPR = 25, & 
               NENBDR = 21, & 
               NEMOSR = 29, & 
               NEMOPR = 25, & 
               NEMODR = 21, & 
               NETCSR = 29, & 
               NETCPR = 25, & 
               NETCDR = 21, & 
               NERUSR = 30, & 
               NERUPR = 26, & 
               NERUDR = 21, & 
               NERHSR = 30, & 
               NERHPR = 26, & 
               NERHDR = 21, & 
               NEPDSR = 30, & 
               NEPDPR = 26, & 
               NEPDDR = 21, & 
               NEAGSR = 30, & 
               NEAGPR = 26, & 
               NEAGDR = 21, & 
               NECDSR = 30, & 
               NECDPR = 26, & 
               NECDDR = 21, & 
               NEINSR = 30, & 
               NEINPR = 25, & 
               NEINDR = 19, & 
               NESNSR = 30, & 
               NESNPR = 25, & 
               NESNDR = 19, & 
               NESBSR = 30, & 
               NESBPR = 25, & 
               NESBDR = 19, & 
               NETESR = 30, & 
               NETEPR = 25, & 
               NETEDR = 19, & 
               NEISR = 30, & 
               NEIPR = 25, & 
               NEIDR = 19, & 
               NEXESR = 30, & 
               NEXEPR = 25, & 
               NEXEDR = 19, & 
               NECSSR = 36, & 
               NECSPR = 32, & 
               NECSDR = 29, & 
               NEBASR = 37, & 
               NEBAPR = 32, & 
               NEBADR = 26, & 
               NELASR = 38, & 
               NELAPR = 35, & 
               NELADR = 25, & 
               NELAFR = 18, & 
               NECESR = 39, & 
               NECEPR = 35, & 
               NECEDR = 26, & 
               NECEFR = 17, & 
               NEPRSR = 40, & 
               NEPRPR = 36, & 
               NEPRDR = 26, & 
               NEPRFR = 19, & 
               NENDSR = 39, & 
               NENDPR = 36, & 
               NENDDR = 26, & 
               NENDFR = 18, & 
               NEPMSR = 37, & 
               NEPMPR = 36, & 
               NEPMDR = 27, & 
               NEPMFR = 17, & 
               NESMSR = 41, & 
               NESMPR = 36, & 
               NESMDR = 27, & 
               NESMFR = 18, & 
               NEEUSR = 41, & 
               NEEUPR = 38, & 
               NEEUDR = 27, & 
               NEEUFR = 18, & 
               NEGDSR = 41, & 
               NEGDPR = 38, & 
               NEGDDR = 27, & 
               NEGDFR = 19, & 
               NETBSR = 39, & 
               NETBPR = 37, & 
               NETBDR = 26, & 
               NETBFR = 18, & 
               NEDYSR = 41, & 
               NEDYPR = 37, & 
               NEDYDR = 26, & 
               NEDYFR = 17, & 
               NEHOSR = 40, & 
               NEHOPR = 37, & 
               NEHODR = 26, & 
               NEHOFR = 19, & 
               NEERSR = 38, & 
               NEERPR = 37, & 
               NEERDR = 28, & 
               NEERFR = 20, & 
               NETMSR = 41, & 
               NETMPR = 37, & 
               NETMDR = 28, & 
               NETMFR = 18, & 
               NEYBSR = 38, & 
               NEYBPR = 37, & 
               NEYBDR = 28, & 
               NEYBFR = 18, & 
               NELUSR = 39, & 
               NELUPR = 36, & 
               NELUDR = 27, & 
               NELUFR = 18, & 
               NEHFSR = 36, & 
               NEHFPR = 36, & 
               NEHFDR = 26, & 
               NEHFFR = 17, & 
               NETASR = 39, & 
               NETAPR = 37, & 
               NETADR = 26, & 
               NETAFR = 17, & 
               NEWSR = 37, & 
               NEWPR = 37, & 
               NEWDR = 26, & 
               NEWFR = 16, & 
               NERESR = 36, & 
               NEREPR = 35, & 
               NEREDR = 27, & 
               NEREFR = 16, & 
               NEOSSR = 36, & 
               NEOSPR = 38, & 
               NEOSDR = 26, & 
               NEOSFR = 16, & 
               NEIRSR = 36, & 
               NEIRPR = 35, & 
               NEIRDR = 26, & 
               NEIRFR = 16, & 
               NEPTSR = 41, & 
               NEPTPR = 36, & 
               NEPTDR = 25, & 
               NEPTFR = 16, & 
               NEAUSR = 35, & 
               NEAUPR = 36, & 
               NEAUDR = 24, & 
               NEAUFR = 17, & 
               NEHGSR = 39, & 
               NEHGPR = 35, & 
               NEHGDR = 24, & 
               NEHGFR = 16, & 
               NETLSR = 36, & 
               NETLPR = 33, & 
               NETLDR = 25, & 
               NETLFR = 17, & 
               NEPBSR = 36, & 
               NEPBPR = 34, & 
               NEPBDR = 26, & 
               NEPBFR = 16, & 
               NEBISR = 36, & 
               NEBIPR = 33, & 
               NEBIDR = 24, & 
               NEBIFR = 16, & 
               NEPOSR = 36, & 
               NEPOPR = 34, & 
               NEPODR = 25, & 
               NEPOFR = 19, & 
               NEATSR = 37, & 
               NEATPR = 35, & 
               NEATDR = 24, & 
               NEATFR = 19, & 
               NERNSR = 40, & 
               NERNPR = 35, & 
               NERNDR = 23, & 
               NERNFR = 18

      save
      end module aambs_limits

