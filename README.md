# pyTVDI

## Synopsis

This project contains *Python* code for computing the *Temperature-Vegetation Dryness Index* 
for assessing soil moisture and evaporative fraction with remote sensing data combining 
a vegetation index and radiometric surface temperature. 

The project consists of: 

1. lower-level modules with the basic functions needed in TVDI calculation. 

2. higher-level scripts for easily running **pyTVDI** from imagery.

## Installation

Download the project to your local system, enter the download directory and then type

`python setup.py install`. 

if you want to install **pyTVDI** and its low-level modules in your Python distribution. 

The following Python libraries will be required for running pyTVDI:

- Numpy
- GDAL

## Code Example
### High-level example

The easiest way to get a feeling of **pyTVDI** and its configuration is through the provided ipython/jupyter notebooks. 
In a terminal shell, navigate to your working folder and type

- `jupyter notebook TVDIgui.ipynb`. 

In addition, you can also run **TVDI** with the scripts *MAIN_TVDI.py*, 
which will read an input configuration file (defaults are *Config_TVDI.txt*). 
You can edit these configuration files or make a copy to fit your data and site characteristics and either run any of 
these two scripts in a Python GUI or in a terminal shell:

- `python MAIN_TVDI.py <configuration file>.`.
> where \<configuration file> points to a customized configuration file... leave it blank if you want to use the default 
file *Config_TVDI.txt*.


### Low-level example
You can obtain TVDI/EF images by importing the module *tvdi*, which contains all the methods for estimating the dry 
and wet edges and computing TVDI or Evaporative Fraction.

```python
import pyTVDI 
output=pyTVDI.tvdi(io_inf, roi_inf, alg_inf)
```

You can type
`help(pyTVDI.tvdi)`
to understand better the inputs needed and the outputs returned

   
## Basic Contents
### High-level modules
- *.src/pyTVDI.py*, class object for TSEB scripting.
- *TVDI_GUI.ipynb* notebook for using pyTVDI and configuring TSEB through a Graphical User Interface, GUI.
- *MAIN_TVDI.py*, high level scripts for running TVDI through a configuration file (*Config_TVDI.txt*).

### Low-level module
The low-level module in this project is aimed at providing customisation and more flexibility in running the Ts-VI triangle method. 
The following modules are included

- *.src/pyTVDI.py*.
> core functions for running **TVDI**. 

## API Reference
http://pytvdi.readthedocs.org/en/latest/index.html

## Main Scientific References
- Sandholt, I.; Rasmussen, K. & Andersen, J. A simple interpretation of the surface temperature/vegetation 
index space for assessment of surface moisture status. Remote Sensing of Environment , 2002, 79, 213 - 224.

- Stisen S.; Sandholt I.; Nørgaard A.; Fensholt R. & Jensen, K. Combining the triangle method with thermal inertia 
to estimate regional evapotranspiration - Applied to MSG-SEVIRI data in the Senegal River basin. Remote Sensing of Environment, 2008, 112, 1242-1255.

- de Tomás, A.; Nieto, H.; Guzinski, R.; Salas, J.; Sandholt, I. & Berliner, P. Validation and scale dependencies 
of the triangle method for the evaporative fraction estimation over heterogeneous areas. Remote Sensing of Environment , 2014, 152, 493 - 511.


## Tests
The folder *./Input* contains examples for running **pyTVDI**. Just run the high-level scripts with the configuration files 
provided by default and compare the resulting outputs with the files stored in *./Output/*

## Contributors
- **Hector Nieto** <hnieto@ias.csic.es> <hector.nieto.solana@gmail.com> main developer
- **Radoslaw Guzinski** main developer
- **Inge Sandholt** TVDI principal investigator

## License
pyTVDI: a Python Temperature-Vegetation Dryness Index

Copyright 2016 Hector Nieto and contributors.
    
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
