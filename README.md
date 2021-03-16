# ElongatedRiverMesh  
Currently it only supports building elongated mesh for ADCIRC. It is a working product. Take it as it is. No guarantee of the outcome.

### Prerequisite packages  
os\
numpy\
traceback\
pandas\
ogr\
gc

## 1. Prepare the channel geometry input file
(1) (Python ‘parsegeoras’ tool) Extract channel geometry from HECRAS file (e.g. Inputs/011003_NEW_RIVER.g01) and generate the output in an attribute table (e.g. ExampleFiles/XS_YF.csv file). Each cross section should keep a unique identifier such as the cross-section name or station ID.  

(2) (ArcGIS ‘Dissolve’ function) Dissolve the NHD stream centerline into one single polyline (e.g. ExampleFiles/Flowline_New_upper_main_dissolve.shp). So the segments should be merged into one polyline.  

(3) (ArcGIS ‘Join table’ function) Join the HECRAS XS cutlines shapefile (may need to round the numbers in the station ID field) with HECRAS channel geometry (the csv file generated by parsegeoras) by matching the station ID. Note: there might be multiple line segments on one HECRAS XS cutline corresponding to a unique XS ID, so you may need to manually remove the shorter segments that do not intersect the channel centerline.  

(4) (ArcGIS ‘Intersect’ function) Calculate the intersection points between HEC-RAS XS cutlines with the joined channel geometry and the NHD stream centerlines. The information of channel geometry will be preserved on the intersection points.

(5) (ArcGIS ‘Split’ function) Split the NHD stream centerlines by the intersection points.

(6) (ArcGIS ‘Spatial Join’ function) Re-join the split NHD stream centerline with the attributes of channel geometry from the intersection point shapefile.

(7) (Python ‘Sort_Trace_DS’ tool + ArcGIS ‘Sort’ function) Use the ‘Sort_Trace_DS’ tool to trace downstream and compute the order (SortID) of the stream segments. Then use ‘Sort’ function to sort segments of the updated NHD stream centerlines by SortID field.

## 2. Flowlines_BuildXS_weir24_customerized.py  
Use the sorted NHD stream centerline with the channel geometry as the input, run the python script to generate the quadrilaterals that can be read in ADCIRC.

Major parameters:  
(1) num_elems: number of elements in the cross-section direction of the channel.  
(2) flag_weir: turn on or off the weir elements.  
(3) weir_type: weir type in ADCIRC.  
(4) flag_weir_extra_rows: turn on or off building extra rows of elements next to the weir elements.  
(5) flag_smoothchannel: control of the smoothing algorithm for channel curvature.  
(6) flag_smoothbathy: control of the smoothing algorithm for channel bathymetry.  
(7) flag_smoothwidth: control of the smoothing algorithm for channel width.  
(8) flag_smoothbank: control of the smoothing algorithm for channel bank elevation.  
(9) DX: average edge length of channel elements in the along-channel direction.  
(10) DTHETALIM: maximum turning angle of a channel bend in degrees.  
(11) WIDTHLIM: minimum width threshold in meters.  
(12) H0: ADCIRC parameter. The script currently assumes that the weir elevation is 10 times of H0 higher than the bank elevation.  
(13) mindx: minimum distance between two XS in meters. Shorter elements will be removed.  

## 3. writeFort14_1direction.py  
Use the python script to generate the triangular elements that can be read in ADCIRC.  

Major parameters:  
(1) xs_start: the node IDs of the starting XS.  
(2) width_elems: the number of elements in the cross-section direction.  
(3) num_XS: number of XS that will have weirs.  
(4) num_XS_estuary: number of XS that will not have weirs.  
(5) flag_weir_extra_rows: control to have extra rows of elements adjacent to the central channel elements.  
