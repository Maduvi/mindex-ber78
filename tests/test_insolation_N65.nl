! namelist file
! ------------------->,
&INPUT_INFO
 KTYP=               1,  ! 0=mindex, 1=insol, 2=timeslice, 3=orbinput
 KUNI=               0,  ! 0=Wm-2, 1=Langd-1
 TINI=       -500000.0,  ! negative past
 TFIN=             0.0,
 TINC=          1000.0,
 PERH=               0,  ! when KTYP=3
 ECC =               0,  ! when KTYP=3
 XOB =               0,  ! when KTYP=3
 QLAT=            65.0,  ! when KTYP=1,2,3  
 S0=          1361.371,  ! W m-2
 /
