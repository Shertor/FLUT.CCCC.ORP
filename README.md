FLUT.SCSC
============

Console-based program for solving a Flutter of rectangular orthotropic (ORP) plate problem with SCSC boundary conditions _(link to paper coming soon)_. 

## Features

* Calculation of Critical flutter speed
* Calculation of Eigenvalues
* Drawings of Re of eigenfunction can be exported to .txt file

## Building

_coming soon_

## Limitations

### Only square square plates can be used for now

For some reason it is an error in SCSC subprogram, so it isn't stable in calculations on non square plates. Most likely, plate width coefficient is missing in one of the interpolation equations.  

## How to use

* Set the plate parameters and material properties. Recommended to set them in code rather than change in console (deprecated):

  `NUx` - Poissons modulus x direction

  `NUy` - Poissons modulus y direction

  `H` - dimensiolens plate thickness 

  `Ex` - Young modulus x direction, kg/cm**2

  `Ey` - Young modulus y direction, kg/cm**2

  `RO` - plate material density

  `BL` - half of the plate width (set 1.D0 for 2x2 plate)

  `TJ1` and `TJ2` - flow velocity vector angle in from TJ1/TJ2 * pi (eg.: TJ1=10, TJ2=180 -> 10 degrees)

* For question *'Automatic or manual finding of critical speeds (Y/N/A/Q)'* type 'A' 
* Wait for result
