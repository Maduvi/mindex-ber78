! namelist file
! ------------------->,
&INPUT_INFO
 ktyp=               2,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 kuni=               0,  ! 0=Wm-2, 1=Langd-1
 tini=          -10000,  ! negative past, used if ktyp=0,1,2
 tfin=               0,	 ! if ktyp=0,1
 tinc=               0,  ! if ktyp=0,1
 perh=               0,  ! if ktyp=3
 ecc =               0,  ! if ktyp=3
 xob =               0,  ! if ktyp=3
 qlat=              65,  ! if ktyp=1,2,3  
 s0=            1360.0,  ! W m-2
 /
