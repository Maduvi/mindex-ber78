#!/bin/sh
# simply test things

cat << EOF > test_insolation_N65.nl
! namelist file
! ------------------->,
&INPUT_INFO
 KTYP=               1,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 KUNI=               0,  ! 0=Wm-2, 1=Langd-1
 TINI=       -500000.0,  ! negative past, used if ktyp=0,1,2
 TFIN=             0.0,  ! if ktyp=0,1
 TINC=          1000.0,  ! if ktyp=0,1
 PERH=               0,  ! if ktyp=3
 ECC =               0,  ! if ktyp=3
 XOB =               0,  ! if ktyp=3
 QLAT=            65.0,  ! if ktyp=1,2,3  
 S0=          1361.371,  ! W m-2
 /
EOF

cat << EOF > test_insolation_S65.nl
! namelist file
! ------------------->,
&INPUT_INFO
 KTYP=               1,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 KUNI=               0,  ! 0=Wm-2, 1=Langd-1
 TINI=       -500000.0,  ! negative past, used if ktyp=0,1,2
 TFIN=             0.0,  ! if ktyp=0,1
 TINC=          1000.0,  ! if ktyp=0,1
 PERH=               0,  ! if ktyp=3
 ECC =               0,  ! if ktyp=3
 XOB =               0,  ! if ktyp=3
 QLAT=           -65.0,  ! if ktyp=1,2,3  
 S0=          1361.371,  ! W m-2
 /
EOF

echo "Running monsindex.exe"
../bin/monsindex.exe test_insolation_N65.nl  > test_insolation_N65.res
../bin/monsindex.exe test_insolation_S65.nl  > test_insolation_S65.res

echo "Plotting with Python"
python test_insolation.py
