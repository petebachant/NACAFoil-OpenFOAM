# NACA foil OpenFOAM

A [Calkit](https://calkit.org) project that simulates
NACA airfoils with OpenFOAM in a Docker container.

## Usage

To download all of the cases and their results, execute

```sh
calkit pull
```

To run all cases and aggregate the results, execute

```sh
calkit run
```

To run a single case, execute

```sh
calkit xenv python run.py {profile} {angle_of_attack}
```

For example, to simulate the flow around a NACA 0012 at 8 degrees angle
of attack, run

```sh
calkit xenv python run.py 0012 8
```

## Acknowledgements

`blockMeshDict` generation script based on work by
[Håkon Strandenes](https://www.hpc.ntnu.no/display/hpc/OpenFOAM+-+Airfoil+Calculations#OpenFOAM-AirfoilCalculations-3:Calculationofforcesandforcecoefficients).

## Experimental data

`NACA0012_6e6_Ladson*` datasets taken from http://turbmodels.larc.nasa.gov/NACA0012_validation/CLCD_Ladson_expdata.dat,
linked from http://turbmodels.larc.nasa.gov/naca0012_val.html. Header:

```
# Data from Ladson, NASA TM 4074, 1988
# Re=6 million, with transition tripped
# M=0.15
```

## License

Code is GPL licensed. See LICENSE for details.

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
<img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/4.0/88x31.png" />
</a><br />All other materials licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/"/>
Creative Commons Attribution 4.0 International License</a>.
