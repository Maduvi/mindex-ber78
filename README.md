# monsindex-fortran 

Orbital insolation monsoon index calculator in Fortran.

Compute the orbital insolation monsoon index for any time slice past, present or future. This index was defined by Rossignol-Strick (1983). It is a proxy for the orbital forcing of the Earth system in the low latitudes. It is computed as a difference in mean cumulative insolation between the Tropic of Cancer at 23.45°N and the Equator at 0° during the caloric Northern Hemisphere summer half-year (cf. Rossignol-Strick, 1983). The caloric half-years were defined by Milankovitch (1941) to be able to compare objectively summer and winter seasons of different millennia, since astronomical half-years carry the complication of having variable duration across millennia. Here we compute the secular variations of the Earth's orbital parameters using the theory and algorithm of Berger (1978) with 19 coefficients for eccentricity, 18 for obliquity and 9 for general precession longitude. Using the orbital parameters, the daily insolation is computed using code from CLIMBER-2 model (Petoukhov et al., 2000). Daily insolation is accumulated in astronomical half-years and then used for the formulae of Vernerkar (1972) to obtain the Northern Hemisphere caloric half-year summer insolation values and the monsoon index. As done by Rossignol-Strick (1983), the program returns in fact the variations in the monsoon index from its 1950 CE (reference year) value. Units are kept as in Rossignol-Strick (1983) in Langley day-1, which can be converted to W m-2 multiplying by 0.484.

## Compilation

The program is writen in Fortran 90. To obtain the caloric half-year we include a sorting algorithm similar to package [palinsol](https://cran.r-project.org/web/packages/palinsol/) of the R programming language by  M. Crucifx (U. catholique de Louvain, Belgium). The sorting algorithm is a bubble sorting algorithm in Fortran by M. J. Rutter available [here](https://www.mjr19.org.uk/IT/sorts/). To compile the program a `Makefile` is provided for the GNU compiler gfortran. In a terminal go to the program folder and execute the `make` statement:

```bash
$ make
```

Two executable files are created in the `bin` folder: `single.exe` and `transient.exe`.

## Examples

The two executables allow the user two types of applications. Both of them require a Fortran namelist text file as input with the specifications of the calculation. The namelist files have a header with an arrow that can help users use the number of white spaces or format that Fortran is expecting to be able to read the input parameters. If the wrong format is used in the namelist file then results can be weird. Both programs print their output in `stdout`.

### Single time slice with given orbital parameters

In case the user already has orbital parameters at hand and wants to estimate the monsoon index for them, then `single.exe` can be used. In the namelist file the user specifies the `PERH` the longitude of the perihelion minus 180º (e.g., present day is about 102º), `ECC` the eccentricity of Earth's orbit and `XOB` its obliquity in degrees. A sample namelist file `single_21ka.namelist` contains:

```
! namelist file
! ------------------->,
&INPUT_INFO
 PERH=    97.770733785,
 ECC =      0.01878872,
 XOB =     22.78729248,
 /
```

To execute the program then in the command line from the `bin` folder:

```bash
$ ./single.exe single_21ka.namelist
```

This returns:

```
--------------------
Input:
--------------------
PERH:  97.77073669
ECC:   0.01878872
XOB:  22.78729248
====================
Results:
--------------------
Prec. index:        0.002
QsT:              900.261
QsE:              862.637
Mons. index:       -7.470
```

### Transient with secular variations

To compute the monsoon index for a specific time window there is the executable `transient.exe`. Time 0.0 is the reference year 1950 CE. Negative numbers are for past times and positive for the future. In the namelist file the user specifies `TINI` the initial time slice in year units, `TFIN` the final time in year units and `TINC` the time increment in year units. A sample namelist file `transient_past.namelist` contains:

```
! namelist file
! ------------------->,
&INPUT_INFO
 TINI=       -190000.0,
 TFIN=             0.0,
 TINC=          1000.0,
 /
```

To execute the program then in the command line from the `bin` folder:

```bash
$ ./transient.exe transient_past.namelist
```

The full output of this command is put here in the `out` folder in case someone is curious. It begins like:

```
      YEAR      QSTROP      QSEQUA    MINDEX     PRECC       ECC       XOB      PERH
-190000.00      885.32      847.94 -22.66693   0.01550  0.045251  22.46879  44.73389
-189000.00      877.74      839.51 -29.39362   0.02283  0.044643  22.50297  61.36694
```

It shows the year number (`YEAR`), the cumulative insolation during Northern Hemisphere caloric summer at Tropic of Cancer (`QSTROP`), at the Equator (`QSEQUA`), the monsoon index variation (`MINDEX`), the precession index (`PRECC`), and the orbital parameters (`ECC`, `XOB` and `PERH`). All insolation quantities are in Langley day-1 (to convert to W m-2 multiply by 0.484).

## References

Berger, A. (1978). Long-term variations of daily insolation and Quaternary climatic changes. Journal of Atmospheric Sciences, 35(12), 2362-2367. https://doi.org/10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2.

Milankovitch, M. K. (1941). Kanon der Erdbestrahlung und seine Anwendung auf das Eiszeitenproblem. Königlich Serbische Akademie Éditions Speciales, 133, 1-633.

Petoukhov, V., Ganopolski, A., Brovkin, V., Claussen, M., Eliseev, A., Kubatzki, C., & Rahmstorf, S. (2000). CLIMBER-2: a climate system model of intermediate complexity. Part I: model description and performance for present climate. Climate Dynamics, 16(1), 1-17. https://doi.org/10.1007/PL00007919.

Rossignol-Strick, M. (1983). African monsoons, an immediate climate response to orbital insolation. Nature, 304(5921), 46-49. https://doi.org/10.1038/304046a0

Vernekar, A. D. (1972). Long-period global variations of incoming solar radiation. In Long-Period Global Variations of Incoming Solar Radiation (pp. 1-128). American Meteorological Society, Boston, MA.  https://doi.org/10.1007/978-1-935704-34-8_1.
