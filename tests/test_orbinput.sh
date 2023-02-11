#!/bin/sh
# simply test things

cat << EOF > test_orbinput.nl
! namelist file
! ------------------->,
&INPUT_INFO
 ktyp=               3,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 kuni=               0,  ! 0=Wm-2, 1=Langd-1
 tini=          -10000,  ! negative past, used if ktyp=0,1,2
 tfin=               0,	 ! if ktyp=0,1
 tinc=               0,  ! if ktyp=0,1
 perh=           290.5,  ! if ktyp=3
 ecc =          0.0185,  ! if ktyp=3
 xob =            24.1,  ! if ktyp=3
 qlat=              65,  ! if ktyp=1,2,3
 s0=            1360.0,  ! W m-2
 /
EOF

echo "Running monsindex.exe"
../bin/monsindex.exe test_orbinput.nl  > test_orbinput.res
echo "Check output in test_orbinput.res"
