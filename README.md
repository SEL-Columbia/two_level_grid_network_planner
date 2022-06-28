## BACKGROUND
This repo supports two-level grid planning given a shapefile of points. It creates MV and LV .shp files and also returns relevant connectivity metrics : Average LV wire length (m) per structure, Average MV wire length (m) per structure, Average Number of Structures per Transformer andAverage Cost per structure ($USD)

@article{fobi2021scalable,
  title={A scalable framework to measure the impact of spatial heterogeneity on electrification},
  author={Fobi, Simone and Kocaman, Ayse Selin and Taneja, Jay and Modi, Vijay},
  journal={Energy for Sustainable Development},
  volume={60},
  pages={67--81},
  year={2021},
  publisher={Elsevier}
}

## Directory Structure
    
```
├── README.md               <- The top-level README for developers using this project.
├── data
│   └── raw                 <- Sample point locations to create a network design (i.e. low voltage, medium voltage and transfor)
│
├── reports                 <- This directory contains the published paper for this work
│
├── requirements.txt        <- The requirements file for reproducing the analysis environment
│
├── src                     <- Source code for use in this project.
    │
    ├── dataprep            <- Scripts to preprocess data
    │   └── split_data.py   <- Recursively split points to obtain smaller coverage areas
    │
    ├── TLND                <- Scripts to run network design model and reconstruct MV for larger areas
        ├── TwoLevelNetworkDesign.py
        └── mvReconstruction.py
    
```

## Getting started.
- Fork this repo, DO NOT clone it as it will not support pushing by noncollaborators.
- Create a python 3 environment (`conda create -n network python=3.6`)
- Install gdal with working numpy (`conda install gdal numpy pandas`)
- Run the requirements file to install other packages (`pip install -r requirements`)

## Running the Code

- Case 1 (Small number of nodes e.g. <1000): If the shapefile has a small number of points (< 1000), the TLND can be run without further splitting.
    - Run `python src/TLND/TwoLevelNetworkDesign.py`)
- Case 2 (Large number of nodes): When there are many nodes to consider within a shapefile, first split the shapefile to smaller files, then run the TLND and after reconstruct the medium voltage network.
    - Run `python src/dataprep/split_data.py`
    - Next for each generated sub-file run `python src/TLND/TwoLevelNetworkDesign.py`
    - After all the networks are generated run `mvReconstruction.py` to obtain the final network.

## Things to Note
- The provided shapefiles have to be in a CRS that uses meters. There is a function (`convert_to_utm_shp`) in TwoLevelNetworkDesign.py that supports projection from WGS84 to desired crs. The function takes in a csv and outputs a shapefile in the crs.
