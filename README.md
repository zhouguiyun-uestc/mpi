MPI-based implementation of the parallel algorithm

**Manuscript Title**: Parallel identification and filling of depressions in raster digital elevation models

This repository contains the source codes of the MPI-based implementation of the parallel algorithm presented in the manuscript above. These codes were used in performing the tests described in the manuscript.


The codes support floating-point GeoTIFF file format through the GDAL library. Please include GDAL library into your compilation. 

Example usages:

mpiexec -n 2 mpifill  dem.tif dem_filled.tif // use two processors to fill dem.tif. Output is dem_filled.tif
