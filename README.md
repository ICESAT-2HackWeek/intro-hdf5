Intro to HDF5 and ICESat-2 Data Files
-------------------------------------

Instructor: Fernando Paolo (paolofer@jpl.nasa.gov)  
Institution: JPL/Caltech

## Tutorials:  

Part 1: Introduction to HDF5 data model  
Part 2: Reduction of ICESat-2 data files  

## Goals:  

- Familiarize with HDF5 data model  
- Familiarize with HDF5 basic tools  
- Download IS2 files according region of interest  
- Extract variables of interest and filter in time and space  
- Prepare data for large-scale processing  
- Learn simple/generic data parallelization strategy  
- Inspect data files from the command line and with plots  

## Libraries:  

- h5py - HDF5 handling   
- numpy - numeric routines  
- scipy - scientific routines  
- astropy - extra scientific routines
- pyproj - map projection routines   
- joblib - shared-memory parallelization  
- matplotlib - visualization routines  
- glob - pathname pattern expansion  
- argparse - arguments parsing  

## Credits

*The algorithms used in this tutorial are downscaled versions of an open-source Python altimetry package, currently being developed at JPL/Caltech: [captoolkit](https://github.com/fspaolo/captoolkit). This package provides a range of algorithms for common tasks in altimetry data processing. We also use the convenient data download interface from [icepyx](https://github.com/icesat2py/icepyx). 
