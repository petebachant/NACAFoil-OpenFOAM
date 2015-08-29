# NACA foil OpenFOAM case files

OpenFOAM case files for simulating NACA foils. 

## Dependencies

* Python
* NumPy
* pandas
* matplotlib
* [yPlus](https://github.com/petebachant/yPlus)

Note: All Python dependencies are included in the 
[Anaconda Python distribution](http://continuum.io/downloads).

## Usage

To run a single case, execute

    ./Allrun {foil} {alpha_deg}
    
For example, to simulate the flow around a NACA 0012 at 8 degrees angle
of attack, run

    ./Allrun 0012 8
    
To automatically simulate multiple angles of attack, execute

    python paramsweep.py
    
To plot the results for multiple angles of attack, run

    python plot.py 

## Acknowledgements

`blockMeshDict` generation script based on work by 
[Håkon Strandenes](https://www.hpc.ntnu.no/display/hpc/OpenFOAM+-+Airfoil+Calculations#OpenFOAM-AirfoilCalculations-3:Calculationofforcesandforcecoefficients).
