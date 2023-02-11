!
! Module containing routines needed to compute cumulative
! insolation during caloric summer half-year.
!
! Author: Maduvi
! Date: 28-07-2021
!
MODULE calormod

    IMPLICIT NONE

    REAL(KIND=8), PARAMETER  ::   pi = 4.D0 * DATAN(1.D0)
    REAL(KIND=8), PARAMETER  ::   pir = pi / 180.0
    REAL(KIND=8), PARAMETER  ::   pirr = pir / 3600.0
    REAL(KIND=8), PARAMETER  ::   ttrop = 360.0         ! 365.24219876  ! tropical year
    REAL(KIND=8), PARAMETER  ::   hyear = 180.0       ! half a year
    REAL(KIND=8), PARAMETER  ::   Tmin = ttrop * 24.0 * 60.0  ! minutes
    REAL(KIND=8), PARAMETER  ::   conv = 2.06500956023 ! conversion Wm-2 to Lang day-1
CONTAINS

    SUBROUTINE TODAY_CALSUMMER(s0, lat, ecc, xob, perh, Qs)
        ! compute cumulative insolation during caloric summer
        ! for a given latitude and orbital parameters.

        REAL(KIND=8), INTENT(IN)  :: lat, ecc, xob, perh, s0
        REAL(KIND=8), INTENT(OUT) :: Qs

        ! local
        REAL(KIND=8) :: Rs, Rw, m, ss, sw
        REAL(KIND=8), DIMENSION(360) :: solarm

        ! Astronomical cumulative insolation
        CALL SINSOL(s0, lat, ecc, xob, perh, solarm, Rs, Rw, ss, sw)

        ! m coefficient
        m = 2.0 * ttrop * s0 * COS(lat * pir) / (pi**2 * SQRT(1 - ecc**2))

        ! caloric summer
        Qs = (Rs - m * ecc * SIN(perh * pir)) / hyear
        
    END SUBROUTINE TODAY_CALSUMMER

    SUBROUTINE CALSUMMER(s0, lat,  perh, ecc, xob, prec0, xob0, Qs0, Qs)
        ! compute caloric summer insolation as defined by Milankovitch

        REAL(KIND=8), INTENT(IN)  :: s0, lat,  perh, ecc, xob, prec0, xob0, Qs0
        REAL(KIND=8), INTENT(OUT) :: Qs

        ! local
        REAL(KIND=8) :: deps, m, prec
        REAL(KIND=8) :: Rs, dRs, Raux, Rw, Rwaux
        REAL(KIND=8) :: dQs, ss, sw
        REAL(KIND=8), DIMENSION(360) :: solarm
        
        ! delta epsilon
        deps = xob - xob0
        
        ! delta precession index
        prec  = (ecc * SIN(perh * pir)) - prec0
        
        ! coefficient m
        m = 2.0 * ttrop * s0 * COS(lat * pir) / (pi**2 * SQRT(1 - ecc**2))

        ! compute delta per degree in obliquity
        CALL SINSOL(s0, lat, ecc, xob, perh, solarm, Rs, Rw, ss, sw)
        CALL SINSOL(s0, lat, ecc, xob + 1, perh, solarm, Raux, Rwaux, ss, sw)
        dRs = Raux - Rs

        ! change in Q following milankovitch
        dQs = (dRs * deps - prec * m) / hyear

        ! absolute caloric summer
        Qs = Qs0 + dQs
        
    END SUBROUTINE CALSUMMER
    
    SUBROUTINE SEASONS_LENGTH(ecc, perh, Ts, Tw)
        ! estimate length of astronomical summer half-year.

        REAL(KIND=8), INTENT(IN)  :: ecc, perh
        REAL(KIND=8), INTENT(OUT) :: Ts, Tw

        Ts = (0.5D0 * ttrop) * (1.D0 + (4.D0 * ecc * SIN(perh * pir) / pi))
        Tw = ttrop - Ts
        
    END SUBROUTINE SEASONS_LENGTH

    SUBROUTINE SINSOL(s0, lat, ecc, xob, perh, solarm, astrosum, astrowin, ss, sw)
        ! original by Kubatzki and Ganopolski,
        ! computes daily solar insolation for a year
        ! given a latitude and orbital parameters.
        ! modified from original part of CLIMBER-2.
        ! modified by maduvi to accumulate astronomical seasons.
        
        IMPLICIT NONE
        
        REAL(KIND=8), INTENT(IN)                  :: lat, ecc, xob, perh, s0
        REAL(KIND=8), INTENT(OUT)                 :: astrosum, astrowin, ss, sw
        REAL(KIND=8), DIMENSION(360), INTENT(OUT) :: solarm

        ! local
        INTEGER      :: nd, m, nh
        REAL(KIND=8) :: fi
        REAL(KIND=8) :: tperi, zavexpe, pclock
        REAL(KIND=8) :: pytime, pdisse, pzen1, pzen2, pzen3, prae
        REAL(KIND=8) :: cosp, cosn, s
        REAL(KIND=8) :: Ts, Tw
        REAL(KIND=8) :: dfracs

        ! latitude in radians
        fi = lat * pir
      
        CALL BERGOR(ecc, perh, tperi, zavexpe)
        
        CALL SEASONS_LENGTH(ecc, perh, Ts, Tw)

        dfracs = MOD(Ts, 1.0)  ! get decimal parts (day fractions)
        astrosum = 0.0
        astrowin = 0.0
        
        ! Daily insolation is calculated by hourly integration for each day
        
        DO nd = 1, 360
            m=nd
            pytime = nd * 2.0 * pi / 360.0
            solarm(m) = 0.0
            DO nh = 1, 24
                pclock = nh * 2.0 * pi / 24.0

                CALL ORBIT(ecc, xob, tperi, zavexpe, pclock,&
                    &pytime, pdisse, pzen1, pzen2, pzen3, prae)

                cosp = pzen1 * SIN(fi) + pzen2 * COS(fi)
                cosn = AMAX1(cosp, 0.0)
       
                s = s0 * cosn * pdisse
                solarm(m) = solarm(m) + s
                
                ! astronomical summer calculation
                if (nd > 80 .and. nd <= 80 + FLOOR(Ts)) then
                    astrosum = astrosum + s
                end if
                ! add little fractions lost in flooring
                if (nd == 80 + FLOOR(Ts)) then
                    astrosum = astrosum + (s * dfracs)
                    astrowin = astrowin + (s * (1.D0 - dfracs))
                end if
                
                ! astronomical winter calculation
                if (nd < 80 .or. nd > 80 + FLOOR(Ts)) then
                    astrowin = astrowin + s
                end if

            END DO
        END DO
        
        ! Daily insolation and zenith angle
        DO m = 1, 360
            solarm(m) = solarm(m) / 24.0
        END DO

        ! average per day
        astrosum = astrosum / 24.0
        astrowin = astrowin / 24.0
        ss = MAXVAL(solarm)  ! solstices are extreme values
        sw = MINVAL(solarm)

    END SUBROUTINE SINSOL
    
    SUBROUTINE ORBIT(ecc, xobch, tperi, zavexpe, pclock, pytime,&
        &pdisse, pzen1, pzen2, pzen3, prae)
        ! original by Geleyn.
        ! computes orbital parameters required for solar insolation.
        ! modified from original part in CLIMBER-2.
        
        IMPLICIT NONE
        
        ! global
        REAL(KIND=8), INTENT(IN)  :: ecc, xobch, tperi, zavexpe, pclock, pytime
        REAL(KIND=8), INTENT(OUT) :: pdisse, pzen1, pzen2, pzen3, prae
        
        ! local
        INTEGER      :: niter
        REAL(KIND=8) :: zclock, zytime, time, eold, enew, eps
        REAL(KIND=8) :: cose, zeps, e, zsqecc, ztgean, znu, zlambda
        REAL(KIND=8) :: zsinde, zdecli, xobche, zdisse, zzen1, zzen2, zzen3
        REAL(KIND=8), PARAMETER :: zrae = +0.1277E-02
        
        zclock = pclock
        zytime = pytime
        
        time = zytime - 2.0 * pi * tperi / ttrop
        eold = time / (1.0 - ecc)
        enew = time
        eps = 1.E-6
        niter=0
        
        zeps = eold - enew
        DO WHILE (ABS(zeps) > eps)
            niter = niter + 1
            IF (niter >= 30) THEN
                print*, ' SUBROUTINE *ORBIT* - eccentric anomaly not found!'
                print*, ' ERROR IN   *ORBIT* -- STOP'
                stop
            END IF
          
            eold = enew
            cose = COS(enew)
            enew = (time + ecc * (SIN(enew) - enew * cose)) / (1.0 - ecc * cose)
            zeps = eold - enew
        END DO
      
        e = enew
        zdisse = (1.0 / (1.0 - ecc * COS(e)))**2
        zsqecc = SQRT((1 + ecc) / (1 - ecc))
        ztgean = TAN(e / 2)
        znu = 2.0 * ATAN(zsqecc * ztgean)
        zlambda = znu + zavexpe
        xobche = xobch * pir
        zsinde = SIN(xobche) * SIN(zlambda)
        zdecli = ASIN(zsinde)
        
        zzen1 = SIN(zdecli)
        zzen2 = COS(zdecli) * COS(zclock)
        zzen3 = COS(zdecli) * SIN(zclock)
        
        pdisse = zdisse
        pzen1 = zzen1
        pzen2 = zzen2
        pzen3 = zzen3
        prae = zrae
        
    END SUBROUTINE  ORBIT

    SUBROUTINE BERGOR(ecc, perh, tper, zan)
        ! original by Lorenz.
        ! computes day of perihelion in year and zenith angle.
        ! modified from original part in CLIMBER-2.
        
        IMPLICIT NONE
        
        REAL(KIND=8), INTENT(IN)  :: ecc, perh
        REAL(KIND=8), INTENT(OUT) :: tper, zan
        
        ! local
        REAL(KIND=8) :: ang, tgnu, sqecc, e, tian, days, epc
        
        ! calculate time of perihelion in year:
        ang = perh - 180.0
        
        IF (ang < 0.0) THEN
            ang = ang + 360.0
        END IF
        
        zan = ang * pir
        tgnu = TAN(zan / 2.0)
        sqecc = SQRT((1 - ecc) / (1 + ecc))
        e = 2.0 * ATAN(sqecc * tgnu) ! eccentric Anomaly
        
        ! time angle in radians of perihelion from vernal equinox
        tian = e - ecc * SIN(e)
        days = tian / pir !   360 day year only
        IF (days.lt.0.0) THEN
            days = days + 360.0 !   days from ver.eq. to perh.
        END IF
        
        ! time in days from begin of year: vernal eq. fixed at 3/21, 12 GMT
        ! = 80.5 days in 360 day year
        tper = days + 80.5
        IF (tper.gt.360.0) THEN
            tper = tper - 360.0
        END IF

    END SUBROUTINE BERGOR
    
    SUBROUTINE BERGER(t, perh, ecc, xob)
        ! original Berger (1978) routine. compute:
        ! eccentricity (19 coeff), obliquity (deg) (18 coeff),
        ! general precession longitude (9 coeff).
        ! modified from original part in CLIMBER-2.
        
        IMPLICIT NONE
        
        ! global
        REAL(KIND=8), INTENT(IN)    :: t
        REAL(KIND=8), INTENT(out)   :: perh, ecc, xob
        
        ! local
        INTEGER                      :: i
        INTEGER, PARAMETER           :: nef=19, nob=18, nop=9
        REAL(KIND=8), PARAMETER      :: xod=23.320556
        REAL(KIND=8), PARAMETER      :: xop=3.392506
        REAL(KIND=8), PARAMETER      :: prm = 50.439273
        REAL(KIND=8)                 :: arg, xes, xec
        REAL(KIND=8)                 :: tra, rp, prg
        REAL(KIND=8), DIMENSION(nef) :: AE, BE ,CE
        REAL(KIND=8), DIMENSION(nob) :: AOB, BOB, COB
        REAL(KIND=8), DIMENSION(nop) :: AOP, BOP, COP
        
        ! 1.earth orbital elements:
        
        ! eccentricity
        AE(1) =  0.01860798
        AE(2) =  0.01627522
        AE(3) = -0.01300660
        AE(4) =  0.00988829
        AE(5) = -0.00336700
        AE(6) =  0.00333077
        AE(7) = -0.00235400
        AE(8) =  0.00140015
        AE(9) =  0.00100700
        AE(10) =  0.00085700
        AE(11) =  0.00064990
        AE(12) =  0.00059900
        AE(13) =  0.00037800
        AE(14) = -0.00033700
        AE(15) =  0.00027600
        AE(16) =  0.00018200
        AE(17) = -0.00017400
        AE(18) = -0.00012400
        AE(19) =  0.00001250
        BE(1) =  4.2072050 * pirr
        BE(2) =  7.3460910 * pirr
        BE(3) = 17.8572630 * pirr
        BE(4) = 17.2205460 * pirr
        BE(5) = 16.8467330 * pirr
        BE(6) =  5.1990790 * pirr
        BE(7) = 18.2310760 * pirr
        BE(8) = 26.2167580 * pirr
        BE(9) =  6.3591690 * pirr
        BE(10) = 16.2100160 * pirr
        BE(11) =  3.0651810 * pirr
        BE(12) = 16.5838290 * pirr
        BE(13) = 18.4939800 * pirr
        BE(14) =  6.1909530 * pirr
        BE(15) = 18.8677930 * pirr
        BE(16) = 17.4255670 * pirr
        BE(17) =  6.1860010 * pirr
        BE(18) = 18.4174410 * pirr
        BE(19) =  0.6678630 * pirr
        CE(1) =  28.620089 * pir
        CE(2) = 193.788772 * pir
        CE(3) = 308.307024 * pir
        CE(4) = 320.199637 * pir
        CE(5) = 279.376984 * pir
        CE(6) =  87.195000 * pir
        CE(7) = 349.129677 * pir
        CE(8) = 128.443387 * pir
        CE(9) = 154.143880 * pir
        CE(10) = 291.269597 * pir
        CE(11) = 114.860583 * pir
        CE(12) = 332.092251 * pir
        CE(13) = 296.414411 * pir
        CE(14) = 145.769910 * pir
        CE(15) = 337.237063 * pir
        CE(16) = 152.092288 * pir
        CE(17) = 126.839891 * pir
        CE(18) = 210.667199 * pir
        CE(19) =  72.108838 * pir
      
        ! obliquity
        AOB(1) = -2462.2214466
        AOB(2) =  -857.3232075
        AOB(3) =  -629.3231835
        AOB(4) =  -414.2804924
        AOB(5) =  -311.7632587
        AOB(6) =   308.9408604
        AOB(7) =  -162.5533601
        AOB(8) =  -116.1077911
        AOB(9) =   101.1189923
        AOB(10) =   -67.6856209
        AOB(11) =    24.9079067
        AOB(12) =    22.5811241
        AOB(13) =   -21.1648355
        AOB(14) =   -15.6549876
        AOB(15) =    15.3936813
        AOB(16) =    14.6660938
        AOB(17) =   -11.7273029
        AOB(18) =    10.2742696
        BOB(1) = 31.609974 * pirr
        BOB(2) = 32.620504 * pirr
        BOB(3) = 24.172203 * pirr
        BOB(4) = 31.983787 * pirr
        BOB(5) = 44.828336 * pirr
        BOB(6) = 30.973257 * pirr
        BOB(7) = 43.668246 * pirr
        BOB(8) = 32.246691 * pirr
        BOB(9) = 30.599444 * pirr
        BOB(10) = 42.681324 * pirr
        BOB(11) = 43.836462 * pirr
        BOB(12) = 47.439436 * pirr
        BOB(13) = 63.219948 * pirr
        BOB(14) = 64.230478 * pirr
        BOB(15) =  1.010530 * pirr
        BOB(16) =  7.437771 * pirr
        BOB(17) = 55.782177 * pirr
        BOB(18) =  0.373813 * pirr
        COB(1) = 251.9025 * pir
        COB(2) = 280.8325 * pir
        COB(3) = 128.3057 * pir
        COB(4) = 292.7252 * pir
        COB(5) =  15.3747 * pir
        COB(6) = 263.7951 * pir
        COB(7) = 308.4258 * pir
        COB(8) = 240.0099 * pir
        COB(9) = 222.9725 * pir
        COB(10) = 268.7809 * pir
        COB(11) = 316.7998 * pir
        COB(12) = 319.6024 * pir
        COB(13) = 143.8050 * pir
        COB(14) = 172.7351 * pir
        COB(15) =  28.9300 * pir
        COB(16) = 123.5968 * pir
        COB(17) =  20.2082 * pir
        COB(18) =  40.8226 * pir

        ! general precession in longitude
        AOP(1) =  7391.0225890
        AOP(2) =  2555.1526947
        AOP(3) =  2022.7629188
        AOP(4) = -1973.6517951
        AOP(5) =  1240.2321818
        AOP(6) =   953.8679112
        AOP(7) =  -931.7537108
        AOP(8) =   872.3795383
        AOP(9) =   606.3544732
        BOP(1) =  31.609974 * pirr
        BOP(2) =  32.620504 * pirr
        BOP(3) =  24.172203 * pirr
        BOP(4) =   0.636717 * pirr
        BOP(5) =  31.983787 * pirr
        BOP(6) =   3.138886 * pirr
        BOP(7) =  30.973257 * pirr
        BOP(8) =  44.828336 * pirr
        BOP(9) =   0.991874 * pirr
        COP(1) = 251.9025 * pir
        COP(2) = 280.8325 * pir
        COP(3) = 128.3057 * pir
        COP(4) = 348.1074 * pir
        COP(5) = 292.7252 * pir
        COP(6) = 165.1686 * pir
        COP(7) = 263.7951 * pir
        COP(8) =  15.3747 * pir
        COP(9) =  58.5749 * pir

        ! 3.numerical value for ecc pre xob
        ! t is negative for the past      
        xes = 0.0
        xec = 0.0
        DO  i = 1, nef
            arg = BE(i) * t + CE(i)
            xes = xes + AE(i) * SIN(arg)
            xec = xec + AE(i) * COS(arg)
        END DO

        ! eccentricity
        ecc = SQRT(xes * xes + xec * xec)

        tra = ABS(xec)

        ! adjusting rp if needed
        IF (tra <= 1.0E-08) THEN
            IF (xes < 0) then
                rp = 1.5 * pi
            ELSE IF (xes == 0) THEN
                rp = 0.0
            ELSE
                rp = pi / 2.0
            END IF
        ELSE
            rp = ATAN(xes / xec)

            IF (xec < 0) THEN
                rp = rp + pi
            ELSE IF (xec == 0) THEN
                IF (xes < 0) then
                    rp = 1.5 * pi
                ELSE IF (xes == 0) THEN
                    rp = 0.0
                ELSE
                    rp = pi / 2.0
                END IF
            ELSE
                IF (xes < 0) THEN
                    rp = rp + 2.0 * pi
                END IF
            END IF
        END IF

        ! perihelion longitude
        perh = rp / pir

        ! perigee
        prg = prm * t

        DO i = 1, nop
            arg = BOP(i) * t + COP(i)
            prg = prg + AOP(i) * SIN(arg)
        END DO

        prg = (prg / 3600.0) + xop

        ! omega
        perh = perh + prg

        DO WHILE (perh < 0)
            perh = perh + 360.0
        END DO
        DO WHILE (perh > 360)
            perh = perh - 360.0
        END DO

        ! obliquity
        xob = xod
        DO i = 1, nob
            arg = BOB(i) * t + COB(i)
            xob = XOB + AOB(i) / 3600.0 * COS(arg)
        END DO

    END SUBROUTINE BERGER

    SUBROUTINE BERGER_WDATA(t, perh, ecc, xob)
        ! original Berger (1978) routine.
        ! modified from original part in CLIMBER-2.
        ! only kept here in case debugging using
        ! coefficients from text files like bre78.dat
        
        IMPLICIT NONE
        
        ! global
        REAL(KIND=8), INTENT(IN)    :: t
        REAL(KIND=8), INTENT(out)   :: perh, ecc, xob
        
        ! local
        INTEGER                :: i, dum
        INTEGER, PARAMETER     :: nef=19, nob=47, nop=78
        REAL(KIND=8), PARAMETER      :: xod=23.32054929
        REAL(KIND=8), PARAMETER      :: prm=50.43928443
        REAL(KIND=8), PARAMETER      :: xop=3.39251498
        REAL(KIND=8)                 :: arg, xes, xec
        REAL(KIND=8)                 :: tra, rp, prg
        REAL(KIND=8), DIMENSION(nef) :: AE, BE ,CE
        REAL(KIND=8), DIMENSION(nob) :: AOB, BOB, COB, DOB
        REAL(KIND=8), DIMENSION(nop) :: AOP, BOP, COP, DOP
        CHARACTER       :: dummy

        ! 1.earth orbital elements :

        OPEN(1, file='bre78_trunc.dat', access='sequential', form='formatted', action='read')

        ! header with info
        READ(1, *) dummy
        DO i = 1, nef
            READ(1, '((I4), 3(F20.10))') dum, AE(i), BE(i), CE(i)
        END DO
        BE = BE * pirr
        CE = CE * pir
        DO i = 1, nob
            READ(1, '((I5), (F15.7), (F12.6), (F12.4), (F12.0))') dum, AOB(i), BOB(i), COB(i), DOB(i)
        END DO
        BOB = BOB * pirr
        COB = COB * pir
        DO i = 1, nop
            READ(1, '((I5), (F15.7), (F12.6), (F12.4), (F12.0))') dum, AOP(i), BOP(i), COP(i), DOP(i)
        END DO
        BOP = BOP * pirr
        COP = COP * pir
        CLOSE(1)
        
        ! 3.numerical value for ecc pre xob
        ! t is negative for the past      
        xes = 0.0
        xec = 0.0
        DO  i = 1, nef
            arg = BE(i) * t + CE(i)
            xes = xes + AE(i) * SIN(arg)
            xec = xec + AE(i) * COS(arg)
        END DO

        ! eccentricity
        ecc = SQRT((xes * xes) + (xec * xec))
        tra = ABS(xec)

        ! adjusting rp if needed
        IF (tra <= 1.0E-08) THEN
            IF (xes < 0) then
                rp = 1.5 * pi
            ELSE IF (xes == 0) THEN
                rp = 0.0
            ELSE
                rp = pi / 2.0
            END IF
        ELSE
            rp = ATAN(xes / xec)

            IF (xec < 0) THEN
                rp = rp + pi
            ELSE IF (xec == 0) THEN
                IF (xes < 0) then
                    rp = 1.5 * pi
                ELSE IF (xes == 0) THEN
                    rp = 0.0
                ELSE
                    rp = pi / 2.0
                END IF
            ELSE
                IF (xes < 0) THEN
                    rp = rp + 2.0 * pi
                END IF
            END IF
        END IF

        ! perihelion longitude
        perh = rp / pir

        ! perigee
        prg = prm * t

        DO i = 1, nop
            arg = BOP(i) * t + COP(i)
            prg = prg + AOP(i) * SIN(arg)
        END DO

        prg = (prg / 3600.0) + xop

        ! omega
        perh = perh + prg

        DO WHILE (perh < 0)
            perh = perh + 360.0
        END DO
        DO WHILE (perh > 360)
            perh = perh - 360.0
        END DO

        ! obliquity
        xob = xod
        DO i = 1, nob
            arg = (BOB(i) * t) + COB(i)
            xob = XOB + (AOB(i) / 3600.0) * COS(arg)
        END DO

    END SUBROUTINE BERGER_WDATA

END MODULE calormod
