#!/bin/bash 

mkdir -p EPSG4326

for id_full in *.tif 
 do 
    id=$( echo "$id_full" | cut -c1-43 )
    echo $id
    # convert from pseudo mercator to WGS84
    gdalwarp -overwrite -s_srs EPSG:32639 -t_srs EPSG:4326 -of GTiff $id'.tif' $id'_4326.tif'
    # convert from geotiff to xyz
    gdal_translate -of XYZ $id'_4326.tif' 'B'$id'.xyz'
    # convert from xyz to grid
    $xyz2nc -input 'B'$id'.xyz' -output $id'_4326.nc' -smoothspval -spval 0.0
    rm *.xyz
    mv *.nc EPSG4326/
    rm *.IMD
    rm *4326.tif
 done


