#!/bin/sh
# simply test things

cat << EOF > test_mindex.nl
! namelist file
! ------------------->,
&INPUT_INFO
 ktyp=               0,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 kuni=               1,  ! 0=Wm-2, 1=Langd-1
 tini=       -500000.0,  ! negative past, used if ktyp=0,1,2
 tfin=             0.0,  ! if ktyp=0,1
 tinc=          1000.0,  ! if ktyp=0,1
 perh=               0,  ! if ktyp=3
 ecc =               0,  ! if ktyp=3
 xob =               0,  ! if ktyp=3
 qlat=               0,  ! if ktyp=1,2,3  
 s0=            1360.0,  ! W m-2
 /
EOF

echo "Running monsindex.exe"
../bin/monsindex.exe test_mindex.nl  > test_mindex.res

echo "Plotting with Python"
python test_mindex.py
