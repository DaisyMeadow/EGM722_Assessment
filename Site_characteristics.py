import os
import numpy as np
import rasterio as rio
import pandas as pd
import geopandas as gpd
import rasterstats

def find_underlying_vector(starting_objects, starting_objects_identifying_column, objects_to_select,
                                         objects_to_select_identifying_column):
    # create a new GeoDataFrame using gpd.sjoin()
    joined = gpd.sjoin(starting_objects, objects_to_select, how='inner', lsuffix='left', rsuffix='right')
    # the joined GeoDatFrame contains all columns from both input GeoDataFrame so create a new GeoDataFrame with only
    # the columns required for the output (the two identifying columns)
    results = joined[[starting_objects_identifying_column, objects_to_select_identifying_column]].copy()
    # merge the required output columns back onto the starting objects GeoDataFrame based on the starting object
    # identifying column
    updated = starting_objects.merge(results, on=[starting_objects_identifying_column])
    return updated


# load the input shapefile datasets from the data_files folder using gpd.read_file(os.path.abspath('<file_path>'))
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))

# transform data files to myCRS using gdf.to_crs() (UTM(29) has an epsg of 32629 which will give measurements in metres)
# transform all input data_files to ensure all data is on the same reference system using inplace=true as we want to
# transform the datasets here and not create new ones
sites.to_crs(epsg=32629, inplace=True)
counties.to_crs(epsg=32629, inplace=True)

# find the county each site falls within
# first we need to create a new column with county names that will display better in the output - we wnat to ensure the
# name is no longer all in capitals and we want to add the word 'County' to the county name
for ind, row in counties.iterrows():  # iterate over each row in the GeoDataFrame
    counties.loc[ind, 'County'] = 'County '+row['CountyName'].title()  # assign the row's CountyName to a new column
                                                           # called County (matching the corresponding csv column)
# find the underlying county and add it's name to a new column
sites = find_underlying_vector(sites, 'Name', counties, 'County')

# calculate the areas of each site
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # calculate the area of the polygon and assign it to a new column called Areakm2
    sites.loc[ind, 'Areakm2'] = row['geometry'].area
# divide the area by 1000000 to convert from metres squared to kilometres squared and round to 2 decimal places
sites['Areakm2'] = (sites['Areakm2']/1000000).round(2)

# calculate the perimeter of each site
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # calculate the perimeter of the polygon and assign it to a new column called Perimeter
    sites.loc[ind, 'Perimeter'] = row['geometry'].length
# divide the length by 1000 to convert from metres to kilometres and round to 2 decimal places
sites['Perimeter'] = (sites['Perimeter']/1000).round(2)

# open the land cover raster and read the data
with rio.open('data_files/NI_Landcover_25m_Reclass.tif') as dataset:
    xmin, ymin, xmax, ymax = dataset.bounds
    crs = dataset.crs
    landcover = dataset.read(1)
    affine_tfm = dataset.transform

# we need to ensure the vector layer is in the same crs as the raster before performing the next step
sites.to_crs(crs, inplace=True)

landcover_names = {1: 'Broadleaf woodland',
                   2: 'Coniferous woodland',
                   3: 'Arable',
                   4: 'Improved grassland',
                   5: 'Semi-natural grassland',
                   6: 'Mountain, heath, bog',
                   7: 'Saltwater',
                   8: 'Freshwater',
                   9: 'Coastal',
                   10: 'Built-up areas and gardens'}

# calculate zonal statistics using rasterstats.zonal_stats function
site_stats = rasterstats.zonal_stats(sites, # the shapefile to use
                                       landcover, # the raster to use - here, we're using the numpy array loaded using rasterio
                                       affine=affine_tfm, # the geotransform for the raster
                                       categorical=True, # whether the data are categorical
                                       category_map=landcover_names,
                                       nodata=0 # the nodata value for the raster
                                      )

# create a dictionary with the zonal results added to each site
sites_dict = dict()
for ind, row in sites.iterrows():
    sites_dict[row['Name']] = site_stats[ind]

short_names = ['broadleaf', 'coniferous', 'arable', 'imp_grass', 'nat_grass',
               'mountain', 'saltwater', 'freshwater', 'coastal', 'built_up']
short_dict = dict(zip(landcover_names.values(), short_names)) # use dict and zip with the full names to create a
# dictionary of landcover names and the names for the columns

# add rasterstats results to columns in DataFrame
for ind, row in sites.iterrows(): # use iterrows to iterate over the rows of the table
    sites_data = sites_dict[row['Name']] # get the landcover data for this county
    for name in landcover_names.values(): # iterate over each of the landcover class names
        try:
            sites.loc[ind, short_dict[name]] = sites_data[name] # add the landcover count to a new column
        except KeyError:
            sites.loc[ind, short_dict[name]] = 0 # if name is not present, value should be 0

for ind, row in sites.iterrows(): # iterate over the rows of the table
    sites.loc[ind, short_names] = 100 * row[short_names] / row[short_names].sum()

sites = sites.rename(columns={'Name': 'Site Name'}) # for the sites we will rename the Name column to Site Name so that
                                                    # we can join this data with the proximity data

# create DataFrame to export to csv by copying only desired columns from GeoDataFrame that contains results by dropping
# unnecessary columns (geometry and Id)
site_results = sites.drop(columns=['geometry', 'Id']).copy()

# save the results DataFrame as a csv called Proximity_analysis.csv to the output_files folder with the index removed
site_results.to_csv('output_files/Site_characteristics.csv', index=False)