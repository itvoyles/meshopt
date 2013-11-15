C 
C Written By R.M. Kolonay Oct 2006 
C AFRL/VA C WPAFB OH 45433
C
C Main Routine to interact with the CONMIN program.
C Operates in a pseudo server mode by interacting through input and output 
C streams (ie read from unit 5 and writes to unit 6). This
C Routine is meant to be called by a java application that uses
C a system call to start the process and then interacts with it
C via the input and output streams.
C
C     INTEGER FUNCTION CONMINSERVER( NDV0, NCON0, X, VLB, VUB, 
C    +                              G) 
      INTEGER FUNCTION CONMINSERVER(NDV0, NCON0, X, VLB, VUB, 
     +                              G, A, OBJ0, DF, IPARMS,
     +                              RPARMS, SCAL, S1, G1,
     +                              G2, B3, C4, ISC, IC, MS1,
     +                              N1, N3, IP2, RP2)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C CONMIN COMMONS
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,
     *ALPHAX,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ
     *,ITMAX,ITRM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,IPUNIT,IFUNIT,
     *I1160,I1170,I1180,I1190
C
      COMMON /CONSAV/ DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,DM10,DM11,DM12
     1,DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE,DX,DX1,FI,XI,DFTDF1,A
     2LP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1,CV2,CV3,CV4,APP,ALPCA,ALPFES,AL
     3PLN,ALPMIN,ALPNC,ALPSAV,ALPSID,ALPTOT,RSPACE,IDM1,IDM2,IDM3,JDIR,I
     4OBJ,KOBJ,KCOUNT,NCAL(2),NFEAS,MSCAL,NCOBJ,NVC,KOUNT,ICOUNT,IGOOD1,
     5IGOOD2,IGOOD3,IGOOD4,IBEST,III,NLNC,JGOTO,ISPACE(2)

C
C PARAMETERS
C MAXNDV - Maximum number of design variables that can be handled
C JAXCON - Maximum number of constraints
C N1     - MAXNDV+2
C N2     - JAXCON+2*MAXNDV
C N3     - JAXCON*2+1
C N4     - N3
C N5     - 2*N4
C     PARAMETER ( MAXNDV = 30 )
C     PARAMETER ( JAXCON = 30 )
C     PARAMETER ( N1 = MAXNDV+2)
C     PARAMETER ( N2 = JAXCON+2*MAXNDV )
C     PARAMETER ( N3 = JAXCON*2+1 )
C     PARAMETER ( N4 = N3 )
C     PARAMETER ( N5 = 2*N4 )
C     PARAMETER (NIPSIZ = 20)
C     PARAMETER (NRPSIZ = 20)
C
C     DIMENSION X(*), VLB(*), VUB(*), G(*)
      DIMENSION X(*), VLB(*), VUB(*), G(*), SCAL(*),  
     +          DF(*), A(N1,N3), S1(*), G1(*), G2(*), 
     +          B3(N3,N3), C4(*), ISC(*), IC(*), MS1(*),
     +          IP2(*), RP2(*)
  
      DIMENSION IPARMS(*), RPARMS(*)
C
C
       OPEN(7, File='optiInfoEcho.dat')
C      OPEN(8, File='conminoutput.dat',STATUS='UNKNOWN')
       WRITE(7,*)'Entering CONMIN SERVER'
C      DABFUN=0.0005
       NDV = NDV0
       NCON = NCON0
       OBJ = OBJ0
       IPUNIT = 8
       IFUNIT = 14
C       WRITE(6,*) 'NDV0=',NDV0
C      WRITE(6,*) 'NCON0=',NCON0
C      WRITE(6,*) 'FORTRAN OBJ0=',OBJ0
C
C Read in MINMAX, METHOD, IPRINT, NDV, NCON
C
C the IPARMs and RPARMs
C The order of these in the array are determined by the getINTparms get Doubleparms method
C in the ConminData.java class
C     CNMN1 info
      CALL GETCMN(IPARMS,RPARMS, IP2, RP2)
       N2     = IPARMS(14)
       N4     = IPARMS(16)
       N5     = IPARMS(17)
      WRITE(7,*)'>>>>>>>>>>>>BEFORE CONMIN CALL<<<<<<<<<<<<<<<<<'
      WRITE(7,*) ' CNMN1'
      WRITE(7,*)DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN
      WRITE(7,*)ALPHAX,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT
      WRITE(7,*)NFDG,NSCAL,LINOBJ,ITMAX,ITRM,ICNDIR,IGOTO,NAC
      WRITE(7,*)INFO,INFOG,ITER,IPUNIT
      WRITE(7,*)'CONSAV '
      WRITE(7,*)DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,DM10,DM11,DM12
      WRITE(7,*)DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE,DX,DX1
      WRITE(7,*)FI,XI,DFTDF1,ALP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1
      WRITE(7,*)CV2,CV3,CV4,APP,ALPCA,ALPFES,ALPLN,ALPMIN,ALPNC
      WRITE(7,*)ALPSAV,ALPSID,ALPTOT,RSPACE,IDM1,IDM2,IDM3,JDIR
C     WRITE(7,*)IOBJ,KOBJ,KCOUNT,NCAL(1),NCAL(2),NFEAS,MSCAL
      WRITE(7,*)NCOBJ,NVC,KOUNT,ICOUNT,IGOOD1,IGOOD2,IGOOD3
C     WRITE(7,*)IGOOD4,IBEST,III,NLNC,JGOTO,ISPACE(1),ISPACE(2)
C
       WRITE(7,*) 'in conmin driver IGOTO = ', IGOTO
C
       WRITE(7,*)'CONMIN INFO = ', info
       WRITE(7,*) 'IPRINT, NDV, NCON',IPRINT, NDV, NCON
C
C Echo the OBJ value
       WRITE(7,*)' OBJ = ', OBJ
C
C
C Echo  the constraints
C
C
C
       WRITE(7,*)'X, VLB, VUM, DF, SCAL, S1'
           WRITE(7,26)( X(I), VLB(I), VUB(I), DF(I), 
     +                  SCAL(I), S1(I), I=1,N1)
26         FORMAT(6F16.8)
       WRITE(7,*)'G, G1, G2, ISC'
           WRITE(7,27)(G(I), G1(I), G2(I), ISC(I), I=1,N2)
27         FORMAT(3F16.8,2X,I8)
       WRITE(7,*)'AMAT'
       DO 1002 J = 1, N3
            DO 1012 I = 1, N1
                WRITE(7,40)A(I,J)
1012        CONTINUE
1002   CONTINUE
40     FORMAT(F16.8)
       WRITE(7,*)'BMAT'
       DO 1009 J = 1, N3
            DO 1013 I = 1, N3
                WRITE(7,40)B3(I,J)
1013        CONTINUE
1009   CONTINUE
       WRITE(7,*)'C4'
           WRITE(7,39)(C4(I),I=1,N4)
39         FORMAT(F16.8)
       WRITE(7,*)'IC'
           WRITE(7,31)(IC(I),I=1,N3)
31         FORMAT(I8)
       WRITE(7,*)'MS1C'
           WRITE(7,31)(MS1(I),I=1,N5)
       WRITE(7,*)'N1, N2, N3, N4, N5'
       WRITE(7,*)N1, N2, N3, N4, N5

C Call CONMIN
C      WRITE(6,*)'calling conmin'
      CALL CONMIN( X, VLB, VUB, G, SCAL, DF, A, S1, G1, G2, B3, C4,
     +             ISC, IC, MS1, N1, N2, N3, N4, N5, IPUNIT )

           WRITE(7,*)'>>>>>>>>AFTER CONMIN CALL<<<<<<<<<<<<<<'
       WRITE(7,*)'returned from conmin igoto = ', igoto
       WRITE(7,*)'returned from conmin iter = ', iter
       WRITE(7,*)'conmin info = ', info
C      CLOSE(UNIT=8)
C       
C Write out the updated Objective
C
       OBJ0=OBJ
       WRITE(7,*)'obj after conmin call = ',OBJ
C
C
       WRITE(7,*)'INFO, NAC'
       WRITE(7,50)INFO, NAC
50     FORMAT(2I8)
C
C Write out the new set of design Variables
C
C Echo the DV info
       WRITE(7,*)'X, VLB, VUM, DF, SCAL, S1'
           WRITE(7,26)( X(I), VLB(I), VUB(I), DF(I), 
     +                  SCAL(I), S1(I), I=1,N1)
C Echo  the constraints
       WRITE(7,*)'G, G1, G2, ISC'
           WRITE(7,27)(G(I), G1(I), G2(I), ISC(I), I=1,N2)
       WRITE(7,*)'AMAT'
       DO 1912 J = 1, N3
            DO 1913 I = 1, N1
                WRITE(7,40)A(I,J)
1913        CONTINUE
1912   CONTINUE
       WRITE(7,*)'BMAT'
       DO 1119 J = 1, N3
            DO 1113 I = 1, N3
                WRITE(7,40)B3(I,J)
1113        CONTINUE
1119   CONTINUE
       WRITE(7,*)'C4'
           WRITE(7,39)(C4(I),I=1,N4)
       WRITE(7,*)'IC'
           WRITE(7,31)(IC(I),I=1,N3)
       WRITE(7,*)'MS1C'
           WRITE(7,31)(MS1(I),I=1,N5)
       WRITE(7,*)'N1, N2, N3, N4, N5'
       WRITE(7,*)N1, N2, N3, N4, N5
C
C
C Write out the set of active constraints
C
       WRITE(7,*)'ACTIVE CONSTRAINT INFO'
       WRITE(7,60)(IC(I),  I=1,NAC)
60     FORMAT(I8)
C
C
C Save the common
       CALL SETCMN(IPARMS,RPARMS, IP2, RP2)
C 
      WRITE(7,*)'EXITING CONMIN SERVER ITMAX = ',ITMAX
       WRITE(7,*)'CONMIN INFO = ', info
      CONMINSERVER = INFO
       END
      SUBROUTINE SETCMN(IPARMS,RPARMS, IP2, RP2)
C 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,
     *ALPHAX,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ
     *,ITMAX,ITRM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,IPUNIT, IFUNIT,
     *I1160,I1170,I1180,I1190
C
      COMMON /CONSAV/ DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,DM10,DM11,DM12
     1,DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE,DX,DX1,FI,XI,DFTDF1,A
     2LP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1,CV2,CV3,CV4,APP,ALPCA,ALPFES,AL
     3PLN,ALPMIN,ALPNC,ALPSAV,ALPSID,ALPTOT,RSPACE,IDM1,IDM2,IDM3,JDIR,I
     4OBJ,KOBJ,KCOUNT,NCAL(2),NFEAS,MSCAL,NCOBJ,NVC,KOUNT,ICOUNT,IGOOD1,
     5IGOOD2,IGOOD3,IGOOD4,IBEST,III,NLNC,JGOTO,ISPACE(2)
C23456789012345678901234567890123456789012345678901234567890123456789012
C
      DIMENSION IPARMS(*), RPARMS(*), IP2(*), RP2(*)
C the IPARMs and RPARMs
C The order of these in the array are determined by the getINTparms get Doubleparms method
C in the ConminData.java class
C     CNMN1 info
       IPARMS(1) = NSIDE  
       IPARMS(2) = IPRINT 
       IPARMS(3) = NFDG   
       IPARMS(4) = NSCAL  
       IPARMS(5) =  LINOBJ 
       IPARMS(6) = ITMAX  
       IPARMS(7) = ITRM   
       IPARMS(8) = ICNDIR 
       IPARMS(9) = IGOTO  
       IPARMS(10)= NAC    
       IPARMS(11)= IFOG   
       IPARMS(12)= ITER   
       IPARMS(18)= I1160  
       IPARMS(19)= I1170  
       IPARMS(20)= I1180  
       IPARMS(21)= I1190  
C
       RPARMS(1) = DELFUN 
       RPARMS(2) = DABFUN 
       RPARMS(3) = FDCH   
       RPARMS(4) = FDCHM  
       RPARMS(5) = CT     
       RPARMS(6) = CTMIN  
       RPARMS(7) = CTL    
       RPARMS(8) = CTLMIN 
       RPARMS(9) = ALPHAX 
       RPARMS(10) = ABOBJ1 
       RPARMS(11) = THETA  
C
C
C CONSAV info
C
       RP2(1) =  DM1
       RP2(2) =  DM2
       RP2(3) =  DM3
       RP2(4) =  DM4
       RP2(5) =  DM5
       RP2(6) =  DM6
       RP2(7) =  DM7
       RP2(8) =  DM8
       RP2(9) =  DM9
       RP2(10)=  DM10
       RP2(11)=  DM11
       RP2(12)=  DM12
       RP2(13)=  DCT
       RP2(14)=  DCTL
       RP2(15)=  PHI
       RP2(16)=  ABOBJ
       RP2(17)=  CTA
       RP2(18)=  CTAM
       RP2(19)=  CTBM
       RP2(20)=  OBJ1
       RP2(21)=  SLOPE
       RP2(22)=  DX
       RP2(23)=  DX1
       RP2(24)=  FI
       RP2(25)=  XI
       RP2(26)=  DFTDF1
       RP2(27)=  ALP
       RP2(28)=  FFF
       RP2(29)=  A1
       RP2(30)=  A2
       RP2(31)=  A3
       RP2(32)=  A4
       RP2(33)=  F1
       RP2(34)=  F2
       RP2(35)=  F3
       RP2(36)=  F4
       RP2(37)=  CV1
       RP2(38)=  CV2
       RP2(39)=  CV3
       RP2(40)=  CV4
       RP2(41)=  APP
       RP2(42)=  ALPCA
       RP2(43)=  ALPFES
       RP2(44)=  ALPLN
       RP2(45)=  ALPMIN
       RP2(46)=  ALPNC
       RP2(47)=  ALPSAV
       RP2(48)=  ALPSID
       RP2(49)=  ALPTOT
       RP2(50)=  RSPACE
C
       IP2(1) = IDM1
       IP2(2) = IDM2
       IP2(3) = IDM3
       IP2(4) = JDIR
       IP2(5) = IOBJ
       IP2(6) = KOBJ
       IP2(7) = KCOUNT
       IP2(8) = NCAL(1)
       IP2(9) = NCAL(2)
       IP2(10)= NFEAS
       IP2(11)= MSCAL
       IP2(12)= NCOBJ
       IP2(13)= NVC
       IP2(14)= KOUNT
       IP2(15)= ICOUNT
       IP2(16)= IGOOD1
       IP2(17)= IGOOD2
       IP2(18)= IGOOD3
       IP2(19)= IGOOD4
       IP2(20)= IBEST
       IP2(21)= III
       IP2(22)= NLNC
       IP2(23)= JGOTO
       IP2(24)= ISPACE(1)
       IP2(25)= ISPACE(2)
C
      RETURN 
      END
C
      SUBROUTINE GETCMN(IPARMS,RPARMS, IP2, RP2)
C 
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON /CNMN1/ DELFUN,DABFUN,FDCH,FDCHM,CT,CTMIN,CTL,CTLMIN,
C23456789012345678901234567890123456789012345678901234567890123456789012
     *ALPHAX,ABOBJ1,THETA,OBJ,NDV,NCON,NSIDE,IPRINT,NFDG,NSCAL,LINOBJ
     *,ITMAX,ITRM,ICNDIR,IGOTO,NAC,INFO,INFOG,ITER,IPUNIT, IFUNIT,
     *I1160,I1170,I1180,I1190
C
      COMMON /CONSAV/ DM1,DM2,DM3,DM4,DM5,DM6,DM7,DM8,DM9,DM10,DM11,DM12
     1,DCT,DCTL,PHI,ABOBJ,CTA,CTAM,CTBM,OBJ1,SLOPE,DX,DX1,FI,XI,DFTDF1,A
     2LP,FFF,A1,A2,A3,A4,F1,F2,F3,F4,CV1,CV2,CV3,CV4,APP,ALPCA,ALPFES,AL
     3PLN,ALPMIN,ALPNC,ALPSAV,ALPSID,ALPTOT,RSPACE,IDM1,IDM2,IDM3,JDIR,I
     4OBJ,KOBJ,KCOUNT,NCAL(2),NFEAS,MSCAL,NCOBJ,NVC,KOUNT,ICOUNT,IGOOD1,
     5IGOOD2,IGOOD3,IGOOD4,IBEST,III,NLNC,JGOTO,ISPACE(2)
C
      DIMENSION IPARMS(*), RPARMS(*), IP2(*), RP2(*)
C the IPARMs and RPARMs
C The order of these in the array are determined by the getINTparms get Doubleparms method
C in the ConminData.java class
C     CNMN1 info
       NSIDE  = IPARMS(1);
       IPRINT = IPARMS(2)
       NFDG   = IPARMS(3)
       NSCAL  = IPARMS(4)
       LINOBJ = IPARMS(5)
       ITMAX  = IPARMS(6)
       ITRM   = IPARMS(7)
       ICNDIR = IPARMS(8)
       IGOTO  = IPARMS(9)
       NAC    = IPARMS(10)
       IFOG   = IPARMS(11)
       ITER   = IPARMS(12)
C
       DELFUN = RPARMS(1)
       DABFUN = RPARMS(2)
       FDCH   = RPARMS(3)
       FDCHM  = RPARMS(4)
       CT     = RPARMS(5)
       CTMIN  = RPARMS(6)
       CTL    = RPARMS(7)
       CTLMIN = RPARMS(8)
       ALPHAX = RPARMS(9)
       ABOBJ1 = RPARMS(10)
       THETA  = RPARMS(11)
C
C CONSAV info
C
       DM1 = RP2(1)
       DM2 = RP2(2)
       DM3 = RP2(3)
       DM4 = RP2(4)
       DM5 = RP2(5)
       DM6 = RP2(6)
       DM7 = RP2(7)
       DM8 = RP2(8)
       DM9 = RP2(9)
       DM10= RP2(10)
       DM11= RP2(11)
       DM12= RP2(12)
       DCT = RP2(13)
       DCTL= RP2(14)
       PHI = RP2(15)
       ABOBJ=RP2(16)
       CTA = RP2(17)
       CTAM =RP2(18)
       CTBM =RP2(19)
       OBJ1 =RP2(20)
       SLOPE=RP2(21)
       DX  = RP2(22)
       DX1 = RP2(23)
       FI  = RP2(24)
       XI  = RP2(25)
       DFTDF1=RP2(26)
       ALP = RP2(27)
       FFF = RP2(28)
       A1  = RP2(29)
       A2  = RP2(30)
       A3  = RP2(31)
       A4  = RP2(32)
       F1  = RP2(33)
       F2  = RP2(34)
       F3  = RP2(35)
       F4  = RP2(36)
       CV1 = RP2(37)
       CV2 = RP2(38)
       CV3 = RP2(39)
       CV4 = RP2(40)
       APP = RP2(41)
       ALPCA=RP2(42)
       ALPFES=RP2(43)
       ALPLN=RP2(44)
       ALPMIN=RP2(45)
       ALPNC=RP2(46)
       ALPSAV=RP2(47)
       ALPSID=RP2(48)
       ALPTOT=RP2(49)
       RSPACE=RP2(50)
C
       IDM1 = IP2(1)
       IDM2 = IP2(2)
       IDM3 = IP2(3)
       JDIR = IP2(4)
       IOBJ = IP2(5)
       KOBJ = IP2(6)
       KCOUNT=IP2(7)
       NCAL(1)=IP2(8)
       NCAL(2)=IP2(9)
       NFEAS= IP2(10)
       MSCAL= IP2(11)
       NCOBJ= IP2(12)
       NVC  = IP2(13)
       KOUNT= IP2(14)
       ICOUNT=IP2(15)
       IGOOD1=IP2(16)
       IGOOD2=IP2(17)
       IGOOD3=IP2(18)
       IGOOD4=IP2(19)
       IBEST =IP2(20)
       III   =IP2(21)
       NLNC  =IP2(22)
       JGOTO =IP2(23)
       ISPACE(1)=IP2(24)
       ISPACE(2)=IP2(25)
C
      RETURN
      END
