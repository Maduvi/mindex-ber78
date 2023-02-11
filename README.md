# monsindex-fortran 

## Orbital insolation monsoon forcing index calculator.

Compute the orbital insolation monsoon forcing index for any time slice past, present or future. The index is defined in Rossignol-Strick (1983). It can be a proxy for the orbital forcing of the Earth system at low latitudes. It is computed as a difference in mean insolation between the Tropic of Cancer at 23.45°N and the Equator at 0°, during the caloric Northern Hemisphere summer half-year (cf. Rossignol-Strick, 1983). The caloric half-years (seasons) are defined in Milankovitch (1941). They allow objective comparison of summer and winter seasons across different millennia. Conversely, astronomical half-years (seasons) complicate comparisons because they have variable duration across millennia (e.g., some times astronomical summer lasts 169 days, while others 188). 

To compute the index we need know Earth's orbital configuration. Here we obtain the secular variations of the Earth's orbital parameters using the method of Berger (1978), employing 19 coefficients for eccentricity, 18 for obliquity and 9 for general precession longitude. Using the orbital parameters we compute daily insolation like in CLIMBER-2 climate model (Petoukhov et al., 2000). Daily insolation is then accumulated in astronomical half-years and used for the formulae of Vernerkar (1972) to obtain the Northern Hemisphere caloric half-year summer insolation values and ultimately the monsoon forcing index. Like in Rossignol-Strick (1983), the program returns in fact the variations in the monsoon forcing index from its 1950 CE (reference year) value. Radiation units can be chosen to be as in Rossignol-Strick (1983) (Langley day-1), or can be converted to W m-2 multiplying by about 0.484.

## Compilation

The program is writen in Fortran 90. In case the executable in `bin` folder does not work (compiled with GNU Fortran 8.5.0), then one needs to compile the source code. To compile the program a `Makefile` is provided for the GNU compiler `gfortran`. In a terminal go to the program folder and execute the `make` statement:

```bash
$ make
```

One executable file is created in the `bin` folder: `monsindex.exe`.

## Examples and testing

Please check the contents of folder `tests` for some already-made applications. Tests are run easily using the provided BASH scripts.

For the tests we compare the results of `monsindex.exe` with the output of `palinsol` R package (https://bitbucket.org/mcrucifix/insol). We can only compare insolation quantities, since `palinsol` does not calculate the monsoon forcing index. To test the calculation of the index we have digitized partly the data published in Vernekar (1972) to be able to estimate the original monsoon forcing index values of Rossignol-Strick (1983). Differences between our monsoon index and the original Rossignol-Strick (1983) should be explained in the difference in orbital parameters (we use Berger (1978), and not Vernekar's (1972)).

## Usage
The executable allows the user four types of applications. A Fortran namelist text file is required as input with the specifications of the calculation. Example namelist files have a header with an arrow that helps align white spaces for the format that the Fortran program requires to be able to reliably read input parameters. If a wrong format is used in the namelist file then results can be weird. The program prints all output into `stdout`.

In the namelist files the parameter `ktyp` is the one that defines the type of application:

- `ktyp=0`: computes monsoon forcing index in a transient time series (example: `tests/test_mindex.sh`).
- `ktyp=1`: computes insolation quantities for a chosen latitude (e.g., mean astronomical summer and summer solstice insolation) (example: `tests/test_insolation.sh`).
- `ktyp=2`: computes monsoon forcing index and insolation quantities for a chosen time slice and latitude (example: `tests/test_timeslice.sh`).
- `ktyp=3`: same as `ktyp=3` but instead of a time slice, the input are arbitrary orbital parameters (example: `tests/test_orbinput.sh`).

Parameter `kuni` controls the units. 0 is W m-2 and 1 is Langley day-1. Parameter `qlat` is for the specific latitude in cases when `ktyp=1,2,3`. Parameter `s0` sets the solar constant.

Please look into folder `tests` to see an example of each use case.

## References

Berger, A. (1978). Long-term variations of daily insolation and Quaternary climatic changes. Journal of Atmospheric Sciences, 35(12), 2362-2367. https://doi.org/10.1175/1520-0469(1978)035<2362:LTVODI>2.0.CO;2.

Milankovitch, M. K. (1941). Kanon der Erdbestrahlung und seine Anwendung auf das Eiszeitenproblem. Königlich Serbische Akademie Éditions Speciales, 133, 1-633.

Petoukhov, V., Ganopolski, A., Brovkin, V., Claussen, M., Eliseev, A., Kubatzki, C., & Rahmstorf, S. (2000). CLIMBER-2: a climate system model of intermediate complexity. Part I: model description and performance for present climate. Climate Dynamics, 16(1), 1-17. https://doi.org/10.1007/PL00007919.

Rossignol-Strick, M. (1983). African monsoons, an immediate climate response to orbital insolation. Nature, 304(5921), 46-49. https://doi.org/10.1038/304046a0

Vernekar, A. D. (1972). Long-period global variations of incoming solar radiation. In Long-Period Global Variations of Incoming Solar Radiation (pp. 1-128). American Meteorological Society, Boston, MA.  https://doi.org/10.1007/978-1-935704-34-8_1.

Palinsol (for comparison): https://cran.r-project.org/web/packages/palinsol/.
