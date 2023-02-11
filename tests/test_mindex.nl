! namelist file
! ------------------->,
&INPUT_INFO
 ktyp=               0,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 kuni=               1,  ! 0=Wm-2, 1=Langd-1
 tini=       -500000.0,  ! negative past
 tfin=             0.0,
 tinc=          1000.0,
 perh=               0,  ! when KTYP=3
 ecc =               0,  ! when KTYP=3
 xob =               0,  ! when KTYP=3
 qlat=               0,  ! when KTYP=1,2,3  
 s0=            1360.0,  ! W m-2
 /
