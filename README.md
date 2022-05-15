### Description

Example algorithm for low-lunar orbit maintenance.

### Building

The program can be built and with the [Fortran Package Manager](https://github.com/fortran-lang/fpm) using the provided `fpm.toml` file like so:

```bash
fpm run -- configs/quick.json
```

Where the command line argument is the config file that sets the run settings. See the files in the `configs` directory for examples.

The program has the following dependencies which are automatically downloaded by FPM:
* [Fortran Astrodynamics Toolkit](https://github.com/jacobwilliams/Fortran-Astrodynamics-Toolkit)
* [ddeabm](https://github.com/jacobwilliams/ddeabm)
* [pyplot-fortran](https://github.com/jacobwilliams/pyplot-fortran)
* [json-fortran](https://github.com/jacobwilliams/json-fortran)
* [argv-fortran](https://github.com/jacobwilliams/argv-fortran)

### See also

 * [Low Lunar Orbit Maintenance](https://degenerateconic.com/low-lunar-orbit-maintenance.html) [degenerateconic.com] May 06, 2018
