## 

#1 Read postgis server details
#2 Read config file
#3 Read static driver for spatial information
#4 Convert GRETA to raster files if not done already
#5 Create proxies for downscaling emissions
#6 Downscale emission inventory
#7 Create netcdf4 drivers for PALM

import pandas as pd

## Read postgis server details
from postgres_config import *

## Read configuration file
from config import *
df           = pd.read_csv(fn_teil_fiels, delimiter=',')
teil_fields  = list(df.objectid)

## Read spatial information from static driver
from palm_output import *
origin_time, origin_x, origin_y, origin_lat, origin_lon, res, nx, ny, dx, dy = staticDriver(sim_path, sim_name)

## Import analysis functions
from tools import *

### GRETA Emission Inventory
## Calcualte sector emissions
if not greta_sector_emissions:
    sector_emissions_greta(db, username, host_ad, passp, species, sectors, teil_fields, greta_teil1, greta_teil2)
    
### Preparation for downscaling - clip & rasterize
## GRETA proxy 
if greta_sector_proxy:
    for spec in species:
        inlayer     = "greta_" + spec + "_sectors"
        geom_col    = "shape"
        field_bands = sectors
        ds_res      = 1000.0
        clip_rast_postgis(origin_x, origin_y, res, nx, ny,
                db, username, host_ad, passp,
                inlayer, geom_col, field_bands, ds_res)

## Population density proxy
if population_proxy:
    inlayer     = "population_100m_augsburg"
    geom_col    = "geom"
    field_bands = ["einwohner"]
    ds_res      = 100.0
    clip_rast_postgis(origin_x, origin_y, res, nx, ny,
            db, username, host_ad, passp,
            inlayer, geom_col, field_bands, ds_res)
## OSM Proxy
if osm_proxy:
    crs_from = 4326
    crs_to   = 25832
    buff_ext = 500
    bbox     = bbox_wgs_buffer(crs_from, crs_to, origin_x, origin_y, res, nx, ny, buff_ext)
    osm_file = "./input_data/osm_domain"
    osm_proxy_create(bbox, osm_file)

## CORINE_proxy
if corine_proxy:
    inlayer     = "clc_2018_augsburg"
    geom_col    = "geom"
    field_bands = ["code_18"]
    ds_res      = 100.0
    clip_rast_postgis(origin_x, origin_y, res, nx, ny,
            db, username, host_ad, passp,
            inlayer, geom_col, field_bands, ds_res)


### Downscale emission inventory
if apply_downscaling:
   print("Downscale emissions")
   corine = "clc_2018_augsburg_proxy" 
   greta = "greta_co2_sectors_proxy"
   # function
   downscaling_greta(corine, greta, res, origin_x, origin_y, nx, ny) 


else:
   # if downscaling is not needed
   for spec in species:
       inlayer     = "greta_" + spec + "_sectors_proxy"
       resample_greta(inlayer, res, origin_x, origin_y, nx, ny)

## Create netCDF4 emission drivers
if chemistry_output:
    chem_arr = chem_emissions(sim_origin_time, sim_end_time, ny, nx, chem_driver_spec, palm_sector, country, edgar_dir)
    chemistry_driver(sim_path, sim_name, origin_time, origin_x, origin_y, origin_lat, origin_lon, res, nx, ny, dx, dy,
            chem_driver_spec, chem_arr, sim_origin_time, sim_end_time)
    # add "nox_ratio = composition_nox or  sox_ratio = composition_sox" if using nox or sox

#if aerosol_output:
    #aerosol_driver(sim_path, sim_name, origin_lat, origin_lon, res, nx, ny, emis_categ, compo_name, dx, dy, mass_fracs_array)

print("\nEmissions driver created")
