# monsindex-fortran 

Orbital insolation monsoon index calculator in Fortran.

Compute the orbital insolation monsoon index for any time slice past, present or future. This index was defined by Rossignol-Strick (1983). It is a proxy for the orbital forcing of the Earth system in the low latitudes. It is computed as a difference in mean cumulative insolation between the Tropic of Cancer at 23.45°N and the Equator at 0° during the caloric Northern Hemisphere summer half-year (cf. Rossignol-Strick, 1983). The caloric half-years were defined by Milankovitch (1941) to be able to compare objectively summer and winter seasons of different millennia, since astronomical half-years carry the complication of having variable duration across millennia. Here we compute the secular variations of the Earth's orbital parameters using the theory and algorithm of Berger (1978) with 19 coefficients for eccentricity, 18 for obliquity and 9 for general precession longitude. Using the orbital parameters, the daily insolation is computed using code from CLIMBER-2 model (Petoukhov et al., 2000). Daily insolation is accumulated in astronomical half-years and then used for the formulae of Vernerkar (1972) to obtain the Northern Hemisphere caloric half-year summer insolation values and the monsoon index. As done by Rossignol-Strick (1983), the program returns in fact the variations in the monsoon index from its 1950 CE (reference year) value. Units are kept as in Rossignol-Strick (1983) in Langley day-1, which can be converted to W m-2 multiplying by 0.484.

## Compilation

The program is writen in Fortran 90. To compile the program a `Makefile` is provided for the GNU compiler gfortran. In a terminal go to the program folder and execute the `make` statement:

```bash
$ make
```

One executable file is created in the `bin` folder: `monsindex.exe`.

## Examples

The executable allows the user three types of applications. A Fortran namelist text file is required as input with the specifications of the calculation. Example namelist files have a header with an arrow that can help users use the number of white spaces or format that Fortran is expecting to be able to read the input parameters. If the wrong format is used in the namelist file then results can be weird. The program prints its output into `stdout`.

In the namelist files the parameter `KTYP` is the one that defines the type of application:

- `KTYP = 0`: single time slice with given orbital parameters computation of monsoon index.
- `KTYP = 1`: transient with secular variations computation of monsoon index.
- `KTYP = 2`: (or any other number for that matter) computes insolation quantities for a specific latitude. That is astronomical summer insolation, caloric summer insolation and summer solstice insolation.

Parameter `KUNI` controls the units. 0 is W m-2 and 1 is Langley day-1. Parameter `QLAT` is for the specific latitude in case `KTYP = 2`.

### Single time slice with given orbital parameters

In case the user already has orbital parameters at hand and wants to estimate the monsoon index for them, then `monsindex.exe` can be used. In the namelist file the user specifies the `PERH` the longitude of the perihelion (e.g., present day is about 102º), `ECC` the eccentricity of Earth's orbit and `XOB` its obliquity in degrees. A sample namelist file `single_21ka.namelist` contains:

```
! namelist file
! ------------------->,
&INPUT_INFO
 KTYP=               0,
 KUNI=               1,
 TINI=       -190000.0,
 TFIN=             0.0,
 TINC=          1000.0,
 PERH=    97.770733785,
 ECC =      0.01878872,
 XOB =     22.78729248,
 QLAT=            65.0,
 /
```

First notice the `KTYP` parameter is set for this application at 0, and second notice the units parameter `KUNI` is set for Langley day-1.
The other parameters that are not mentioned are ignored in this case.
To execute the program then in the command line from the `bin` folder:

```bash
$ ./monsindex.exe single_21ka.namelist
```

This returns:

```
Input:

PERH:  97.77073669
ECC:   0.01878872
XOB:  22.78729248
====================
Results:

QsT:              904.768
QsE:              823.071
Mons. index:       -8.674
```

### Transient with secular variations

The program also computes the monsoon index for a specific time window. Time 0.0 is the reference year 1950 CE. Negative numbers are for past times and positive for the future. In the namelist file the user specifies `TINI` the initial time slice in year units, `TFIN` the final time in year units and `TINC` the time increment in year units. A sample namelist file `transient_past.namelist` contains:

```
! namelist file
! ------------------->,
&INPUT_INFO
 KTYP=               1,
 KUNI=               1,
 TINI=       -190000.0,
 TFIN=             0.0,
 TINC=          1000.0,
 PERH=    97.770733785,
 ECC =      0.01878872,
 XOB =     22.78729248,
 QLAT=            65.0,
 /
```
First notice the `KTYP` parameter is set for this application at 1, and second notice the units parameter `KUNI` is set for Langley day-1.
The other parameters that are not mentioned are ignored in this case.
To execute the program then in the command line from the `bin` folder:

```bash
$ ./monsindex.exe transient_past.namelist
```

The full output of this command is put here in the `out` folder in case someone is curious. It begins like:

```
        YEAR      QSTROP      QSEQUA    MINDEX       ECC       XOB      PERH
  -190000.00      889.76      808.86 -24.48115  0.045251  22.46879  44.73389
  -189000.00      882.20      800.37 -31.11467  0.044643  22.50297  61.36694
```

It shows the year number (`YEAR`), the cumulative insolation during Northern Hemisphere caloric summer at Tropic of Cancer (`QSTROP`), at the Equator (`QSEQUA`), the monsoon index variation (`MINDEX`) and the orbital parameters (`ECC`, `XOB` and `PERH`).

### Transient insolation quantities

In case the user wants to see the astronomical, caloric and solstice summer insolation values, `KTYP` can be set to any other number than 0 or 1. For instance, 2.
The important parameter to pay attention in this case is `QLAT` that sets the latitude at which the quantities are computed. Additionally the user sets the desired time window parameters. A sample namelist file `radiation_1Ma.namelist` contains:

```
! namelist file
! ------------------->,
&INPUT_INFO
 KTYP=               2,
 KUNI=               0,
 TINI=      -1000000.0,
 TFIN=             0.0,
 TINC=          1000.0,
 PERH=    97.770733785,
 ECC =      0.01878872,
 XOB =     22.78729248,
 QLAT=            65.0,
 /
```

This will compute the insolation quantities for latitude 65N for the last million years in units of W m-2. Beginning of output looks like

```
        YEAR      RASTRO      QCALOR      RSOLST
 -1000000.00      381.41      412.25      528.32
  -999000.00      381.16      411.86      524.29
```

Where `RASTRO` is the astronomical half-year summer insolation value, `QCALOR` the caloric summer half-year insolation and `RSOLST` that at the summer solstice.

## References

Berger, A. (1978). Long-term variations of daily insolation and Quaternary climatic changes. Journal of Atmospheric Sciences, 35(12), 2362-2367. https://doi.org/10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2.

Milankovitch, M. K. (1941). Kanon der Erdbestrahlung und seine Anwendung auf das Eiszeitenproblem. Königlich Serbische Akademie Éditions Speciales, 133, 1-633.

Petoukhov, V., Ganopolski, A., Brovkin, V., Claussen, M., Eliseev, A., Kubatzki, C., & Rahmstorf, S. (2000). CLIMBER-2: a climate system model of intermediate complexity. Part I: model description and performance for present climate. Climate Dynamics, 16(1), 1-17. https://doi.org/10.1007/PL00007919.

Rossignol-Strick, M. (1983). African monsoons, an immediate climate response to orbital insolation. Nature, 304(5921), 46-49. https://doi.org/10.1038/304046a0

Vernekar, A. D. (1972). Long-period global variations of incoming solar radiation. In Long-Period Global Variations of Incoming Solar Radiation (pp. 1-128). American Meteorological Society, Boston, MA.  https://doi.org/10.1007/978-1-935704-34-8_1.
