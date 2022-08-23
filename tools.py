## 
from pyexpat.errors import XML_ERROR_INCOMPLETE_PE
import sys
import numpy as np
import sqlalchemy
import geopandas as gpd
import rasterio
import shapefile
import gdal

from rasterio.transform import Affine
from rasterio import features
from pyproj import Transformer
from OSMPythonTools.nominatim import Nominatim
from OSMPythonTools.overpass import Overpass, overpassQueryBuilder
##------------------------------------------------------------------------------
def db_connect(db, username, host_ad, passp):
    db_url = sqlalchemy.engine.url.URL(drivername='postgresql+psycopg2', 
            host = host_ad, 
            database = db,
            username = username,
            port = 5432, password = passp)
    cnet   = sqlalchemy.create_engine(db_url)

    return cnet

def NFR_to_GNFR(name):
    sectorname = name
    translation_table = {
            "A_PublicPower" : "1A1a",
            "B_Industry"    : "1A1b,1A1c,1A2a,1A2b,1A2c,1A2d,1A2e,1A2f,1A2gviii,2A1,2A2,2A3,2A5a,2A5b,2A5c,2A6,2B1,2B2,2B3,2B5,2B6,2B7,2B10a,2B10b,2C1,2C2,2C3,2C4,2C5,2C6,2C7a,2C7b,2C7c,2C7d,2D3b,2D3c,2H1,2H2,2H3,2I,2J,2K,2L",
            "C_OtherStationaryComb": "1A4ai,1A4bi,1A4ci,1A5a",
            "D_Fugitives"          : "1B1a,1B1b,1B1c,1B2ai,1B2aiv,1B2av,1B2b,1B2c,1B2d",
            "E_Solvents"           : "2D3a,2D3d,2D3e,2D3f,2D3g,2D3h,2D3i,2G",
            "F_RoadTransport"      : "1A3bi,1A3bii,1A3biii,1A3biv,1A3bv,1A3bvi,1A3bvii",
            "G_Shipping"           : "1A3di(ii),1A3dii",
            "H_Aviation"           : "1A3ai(i),1A3aii(i)",
            "I_OffRoad"            : "1A2gvii,1A3c,1A3ei,1A3eii,1A4aii,1A4bii,1A4cii,1A4ciii,1A5b",
            "J_Waste"              : "5B1,5B2,5C1a,5C1bi,5C1bii,5C1biii,5C1biv,5C1bv,5C1bvi,5C2,5D1,5D2,5D3,5E",
            "K_AgriLivestock"      : "3B1a,3B1b,3B2,3B3,3B4a,3B4d,3B4e,3B4f,3B4gi,3B4gii,3B4giii,3B4giv,3B4h",
            "L_AgriOther"          : "3Da1,3Da2a,3Da2b,3Da2c,3Da3,3Da4,3Db,3Dc,3Dd,3De,3Df,3F,3I",
            "SumAllSectors"        : "sum"
            }
    if sectorname in translation_table:
        return translation_table[sectorname]
    else:
        return None

## SQL queries for the postgis database
#*1 Create new table for each chemical speices with the same strucutre as raster_teil
#*2 Create a new column in the table for each sector
#*3 Sum the sectors and add values to table
def greta_sql(conn, greta_teil1, greta_teil2, ss, sec, greta_fld):
    greta_fld_plain = ' + '.join(greta_fld)
    if ss in greta_teil1:
        teil = "raster_emi_teil1"
    elif ss in greta_teil2:
        teil = "raster_emi_teil2"

    greta_table = "greta_" + ss + "_sectors"
    ## create new table once for each chemical species
    query  = 'CREATE TABLE IF NOT EXISTS "{}" AS (SELECT objectid, id_raster, shape FROM {}); COMMIT;'.format(greta_table, teil)
    conn.execute(query)

    ## create column in new table with sector name
    query2 = 'ALTER TABLE {} ADD COLUMN IF NOT EXISTS {} FLOAT(8); COMMIT;'.format(greta_table, sec)
    conn.execute(query2)

    ## sum sectors and add to you table
    if len(greta_fld_plain)>0:
        query3 = 'UPDATE {} AS a SET {} = ({}) FROM {} AS B WHERE b.objectid = a.objectid; COMMIT;'.format(greta_table,sec,greta_fld_plain,teil)
        conn.execute(query3)

    return greta_table

## For each species, sum GNFR from NFR sectors
def sector_emissions_greta(db, username, host_ad, passp, species, sectors, teil_fields, greta_teil1, greta_teil2):
    print("Calculating sector emissions for GRETA")
    # Connect to postgres and open GRETA inventory
    for ss in species:
        if ss not in greta_teil1 and ss not in greta_teil2:
            sys.exit("Specie:", ss," contained in Species is not included in the Emission Inventory")
        print("\tEmitted species:", ss)
        for sec in sectors:
            print("\tGNFR sector:", sec)
            nfr = NFR_to_GNFR(sec).split(",")

            # list of fields to sum
            gfd = []
            for n in nfr:
                if n == "1A3di(ii)" or n == "1A3ai(i)" or n == "1A3aii(i)":
                    prt = n.split("(")
                    prt2 = prt[1].split(")")
                    gf = 'e_' + prt[0].lower() + '_' + prt2[0] + '__' + ss
                else:
                    gf = 'e_' + n.lower() + '_' + ss

                if gf in teil_fields:
                    gfd.append(gf)
                elif gf not in teil_fields:
                    continue
            # convert to tuple from list and string
            greta_fld = " + ".join(gfd)
            greta_fld = tuple([greta_fld])
        
            ## SQL
            layer_name = greta_sql(conn, greta_teil1, greta_teil2, ss, sec, greta_fld)

    return

def bbox_dims(origin_x, origin_y, res, nx, ny):

    bbox_list = [origin_x, origin_y, origin_x + nx*res, origin_y + ny*res]

    return bbox_list

def polygon_bbox(origin_x, origin_y, res, nx, ny):
    minx = origin_x
    miny = origin_y
    maxx = origin_x + nx*res
    maxy = origin_y + ny*res

    return minx, miny, maxx, maxy

def bbox_wgs_buffer(crs_from, crs_to, origin_x, origin_y, res, nx, ny, buff_ext):
    transformer = Transformer.from_crs(crs_to, crs_from)
    xmin, ymin = transformer.transform(origin_x - buff_ext, origin_y - buff_ext)
    xmax, ymax = transformer.transform((origin_x + nx*res)+ buff_ext , (origin_y + ny*res)+ buff_ext)
    bbox_wgs = [xmin, ymin, xmax, ymax]

    return bbox_wgs

def clip_rast_postgis(origin_x, origin_y, res, nx, ny, db, username, host_ad, passp, inlayer, geom_col, field_bands, ds_res):
    print("Clipping and rasterizing: ", inlayer)
    # Get ploygon bbox
    minx, miny, maxx, maxy = polygon_bbox(origin_x, origin_y, res, nx, ny)

    # Connect to postgis
    conn = db_connect(db, username, host_ad, passp)

    # Query to clip postgis polygon
    query_clip = "SELECT * FROM {} WHERE {} && ST_GeomFromText('SRID=25832;POLYGON(({} {}, {} {}, {} {}, {} {}, {} {}))')".format(inlayer, geom_col, minx, miny, maxx, miny, maxx, maxy, minx, maxy, minx, miny)
    clip_df = gpd.read_postgis(query_clip, conn, geom_col)

    # Get extent, transform & crs  of clip geodataframe
    clip_minx, clip_miny, clip_maxx, clip_maxy = clip_df.total_bounds
    clip_crs    = rasterio.crs.CRS.from_epsg(25832)
    transform   = Affine.translation(clip_minx, clip_miny )* Affine.scale(ds_res, ds_res)

    # Write to geotiff
    with rasterio.open("./input_data/"+ inlayer + "_proxy.tif", 'w+', driver='GTiff',
            height = ((clip_maxy - clip_miny)/ds_res), width = ((clip_maxx - clip_minx)/ds_res),
            count = len(field_bands), dtype = "float64", crs = clip_crs, transform = transform) as rst:
        for nb in range(len(field_bands)):
            out_arr = rst.read(1)
            shapes = ((geom, value) for geom, value in zip(clip_df[geom_col], clip_df[field_bands[nb].lower()]))
            burned = features.rasterize(shapes = shapes, fill=0, out=out_arr, transform = transform, dtype = "float64")
            rst.write_band(nb+1, burned)
            rst.set_band_description(nb+1, field_bands[nb])
        rst.close()

    return

def osm_proxy_create(bbox, osm_file):
    print("Creating OSM Proxy")

    nominatim      = Nominatim()
    overpass       = Overpass()
    #osm_query = ['"motorway"', '"motorway_link"', '"primary"', '"primary_link"', '"secondary"', '"secondary_link"', '"trunk"', '"trunk_link"']
    osm_query = ['"motorway"', '"motorway_link"', '"primary"', '"primary_link"', '"secondary"', '"secondary_link"',
    '"trunk"', '"trunk_link"', '"tertiary"', '"tertiary_link"', '"residential"', '"living_street"']

    w = shapefile.Writer(osm_file + ".shp", shapeType = 3)
    w.field('ID', 'N')
    w.field('key', 'C')
    w.field('typeOfRoad', 'C')
    fid = 0

    for ft in osm_query:
        sql = overpassQueryBuilder(bbox=bbox, elementType = 'way', selector ='"highway"='+ ft, includeGeometry = True)
        result_sql = overpass.query(sql, timeout = 600)
        elements   = result_sql.elements()
        if len(elements)>0:
            firstElement = result_sql.elements()[0]
            geometry_sql = firstElement.geometry()
            for k in range(0,len(elements)):
                geometry = elements[k].geometry()
                coordinates = geometry.get("coordinates")
                if len(coordinates) >1:
                    final_coordinates = []
                    final_coordinates.append(coordinates)
                else:
                    final_coordinates = coordinates
                w.line(final_coordinates)
                w.record(fid, 'Highway', ft)
                fid = fid + 1

    ## Create the PRJ file
    prj  = open(osm_file +".prj", "w")
    epsg = 'GEOGCS["GCS_WGS_1984",'
    epsg += 'DATUM["D_WGS_1984",'
    epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
    epsg += ',PRIMEM["Greenwich",0],'
    epsg += 'UNIT["degree",0.0174532925199433]]'
    prj.write(epsg)
    prj.close()
    
    w.close()

    osm_filename = osm_file + ".shp"
    osm_proxy    = gpd.read_file(osm_filename)
    osm_utm25832 = osm_proxy.to_crs(25832)
    osm_utm25832.to_file(osm_filename, driver="ESRI Shapefile")

    return


def resample_greta(inlayer, res, origin_x, origin_y, nx, ny, algorithm="bilinear"):

    input_raster = "./input_data/"+ inlayer + ".tif"
    output_raster = "./input_data/"+ inlayer + "_resample.tif"
    xmin, ymin, xmax, ymax = polygon_bbox(origin_x, origin_y, res, nx, ny)

    ds = gdal.Warp(output_raster, input_raster, xRes=res, yRes=res, resampleAlg=algorithm, format="GTiff", outputBounds=(xmin,ymin,xmax,ymax),)
    return ds

def downscaling_greta(clc, greta, res, origin_x, origin_y, nx, ny):

    # clip clc to greta proxy extent 
    clc_raster = "./input_data/"+ clc + ".tif"
    clc_raster_clip = "./input_data/"+ clc + "_clip.tif"

    greta_raster = rasterio.open("./input_data/"+ greta + ".tif")
    greta_array = greta_raster.read()
    greta_bounds = greta_raster.bounds
    xmin= greta_bounds[0]
    ymax = greta_bounds[1]
    xmax = greta_bounds[2]
    ymin = greta_bounds[3] 

    clc_ds = gdal.Warp(clc_raster_clip, clc_raster, resampleAlg="bilinear", format="GTiff", outputBounds=(xmin,ymin,xmax,ymax))
    clc_array = clc_ds.ReadAsArray()

    # create empty raster 100x100m res with extent of greta proxy

    geotransform = clc_ds.GetGeoTransform()
    wkt = clc_ds.GetProjection()

    # Create new geotiff file
    driver = gdal.GetDriverByName("GTiff")
    output_file = "./input_data/"+ greta + "_downscaled.tif"

    dst_ds = driver.Create(output_file,
                        clc_ds.GetRasterBand(1).XSize,
                        clc_ds.GetRasterBand(1).YSize,
                        13,
                        gdal.GDT_Float32)

    dst_ds.SetGeoTransform(geotransform)
    dst_ds.SetProjection(wkt)
    
    for i in range(greta_array.shape[0]-1): 
        
        data = np.zeros((dst_ds.GetRasterBand(1).YSize, dst_ds.GetRasterBand(1).XSize)) 
        dst_ds.GetRasterBand(i+1).WriteArray(data)
        
    
    dst_ds = None
    empty_r = rasterio.open(output_file)
    empty_r_array = empty_r.read()

    for row in range(greta_array.shape[1]):
        for col in range(greta_array.shape[2]):

            idx_col_on = col
            idx_col_off = col * 10
            idx_row_on = 
            idx_row_off = row * 10 
            pass








    # count for each sector in greta the corresponding clc value in one cell of greta
    # calculate ratio between number of clc and greta value in one cell 


    # resample to 10x10m 

    # ds_clc = resample_greta(clc, res, origin_x, origin_y, nx, ny, algorithm="near")
    # ds_greta = resample_greta(greta, res, origin_x, origin_y, nx, ny)

    # clc_arr = ds_clc.ReadAsArray() 
    # # clc_arr = clc_arr.astype(int) # bilinear algorithm of gdal warp created some floating values - needs to be converted back in original values
    # greta_arr = ds_greta.ReadAsArray()

    # geotransform = ds_greta.GetGeoTransform()
    # wkt = ds_greta.GetProjection()

    # # Create gtif file
    # driver = gdal.GetDriverByName("GTiff")
    # output_file = "./input_data/"+ greta + "_downscaled.tif"

    # dst_ds = driver.Create(output_file,
    #                     ds_greta.GetRasterBand(1).XSize,
    #                     ds_greta.GetRasterBand(1).YSize,
    #                     13,
    #                     gdal.GDT_Float32)

    # dst_ds.SetGeoTransform(geotransform)
    # dst_ds.SetProjection(wkt)
    
    # # calculations for each sector in greta 
    # # see config file for sectors 1 = A_PublicPower 2 = B_Industry ...
    # for i in range(greta_arr.shape[0]-1): 
        
    #     data = greta_arr[i, ...]

    #     if i == 0:
    #         mask = np.where(clc_arr == 121)
    #         data[mask] = data[mask]*2
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 1: # B_Industry
    #         mask = np.where(clc_arr == 121)
    #         data[mask] = data[mask]*2
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 2: # C_OtherStationaryComb
    #         mask = np.where(clc_arr == 0)
    #         data[mask] = data[mask]
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 3: # D_Fugitives
    #         mask = np.where(clc_arr == 121)
    #         data[mask] = data[mask]*2
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 4: # E_Solvents
    #         mask = np.where(clc_arr == 0)
    #         data[mask] = data[mask]
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 5: # F_RoadTransport
    #         mask = np.where(clc_arr == 0)
    #         data[mask] = data[mask]
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 6: # G_Shipping
    #         mask = np.where(clc_arr == 123)
    #         data[mask] = data[mask]*3
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 7: # H_Aviation
    #         mask = np.where((clc_arr == 124) | (clc_arr == 142))
    #         data[mask] = data[mask]*2
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 8: # I_OffRoad
    #         mask = np.where(clc_arr == 211)
    #         data[mask] = data[mask]*4
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 9: # J_Waste
    #         mask = np.where(clc_arr == 211)
    #         data[mask] = data[mask]*4
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 10: # K_AgriLivestock
    #         mask = np.where(clc_arr == 211)
    #         data[mask] = data[mask]*4
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     elif i == 11: # L_AgriOther
    #         mask = np.where(clc_arr == 211)
    #         data[mask] = data[mask]*4
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)
    #     else:
    #         dst_ds.GetRasterBand(i+1).WriteArray(data)

    return
