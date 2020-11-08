## BACKGROUND
This repo support two-level grid planning given a shapefile of points. It creates MV and LV .shp files and also returns relevant connectivity metrics : Average LV wire length (m) per structure, Average MV wire length (m) per structure, Average NUmber of Structures per Transformer andAverage Cost per structure ($USD)

## Getting started.
Create a python 3 environment (`conda create -n network python=3.6`)
Install gdal with working numpy (`conda install gdal numpy pandas`)
Run the requirements file to install other packages (`pip install -r requirements`)
Run the `TwoLevelNetworkDesign.p` (`python TwoLevelNetworkDesign.py`)

## Things to Note
- The provided shapefiles have to be in a CRS that uses meters. There is a function (`convert_to_utm_shp`) in TwoLevelNetworkDesign.py that supports projection from WGS84 to desired crs. The function takes in a csv and outputs a shapefile in the crs.
