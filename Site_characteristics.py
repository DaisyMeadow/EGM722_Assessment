import os
import numpy as np
import rasterio as rio
import rasterio.features
import pandas as pd
import geopandas as gpd
import rasterstats


def find_underlying_vector_value(starting_objects, starting_objects_identifying_column, objects_to_select,
                                         objects_to_select_identifying_column):
    """
    Find the value of a specified attribute of a vector layer underlying another vector layer.

    - Spatial join between input GeoDataFrames
    - Create GeoDataFrame with only desired results columns
    - Merge results column back onto original input starting GeoDataFrame

    Parameters
    ----------
    starting_objects : GeoDataFrame
        Input object(s)
    starting_objects_identifying_column : column label
        Column name that uniquely identifies starting objects
    objects_to_select : GeoDataFrame
        Selecting objects
    objects_to_select_identifying_column : column label
        Desired attribute to find

    Returns
    -------
    GeoDataFrame
        Updated GeoDataFrame
    """
    # create a new GeoDataFrame using gpd.sjoin()
    joined = gpd.sjoin(starting_objects, objects_to_select, how='inner', lsuffix='left', rsuffix='right')
    # the joined GeoDatFrame contains all columns from both input GeoDataFrame so create a new GeoDataFrame with only
    # the columns required for the output (the two identifying columns)
    results = joined[[starting_objects_identifying_column, objects_to_select_identifying_column]].copy()
    # merge the required output columns back onto the starting objects GeoDataFrame based on the starting object
    # identifying column
    updated = starting_objects.merge(results, on=[starting_objects_identifying_column])
    return updated


def find_percentage_underlying_raster_categories_for_polygons(starting_polygons, starting_polygons_identifying_column,
                                                              raster, raster_category_map, column_names, affine,
                                                              nodata=0):
    """
    Generate percentage cover values based on underlying raster values for each polygon in a GeoDataFrame.

    - Calculate zonal statistics using rasterstats function giving pixel counts for each raster value
    - Create dictionary of results assigning zonal statistics to each polygon in polygon GeoDataFrame
    - Create dictionary of raster values and desired column names
    - Assign zonal stats results to new columns in GeoDataFrame using try...except block giving a 0 value for no cover
    - Convert pixel count values to percentages

    Parameters
    ----------
    starting_polygons : GeoDataFrame
        Input polygon(s)
    starting_polygons_identifying_column : column label
        Column name that uniquely identifies starting polygons
    raster : ndarray
        Input raser - must be categorical
    raster_category_map : dict
        Input raster category map dictionary
    column_names : list of str
        Desired output column names
    affine : Affine instance
        The input raster geotransform
    nodata : int, default 0
        The nodata value for the input raster

    Returns
    -------
    GeoDataFrame
        Updated GeoDataFrame
    """
    # calculate zonal statistics using rasterstats.zonal_stats function
    zonal_stats = rasterstats.zonal_stats(starting_polygons, # the shapefile to use
                                           raster, # the raster to use
                                           affine=affine, # the geotransform for the raster
                                           categorical=True, # this function only runs on categorised data
                                           category_map=raster_category_map, # the raster category map dictionary
                                           nodata=nodata # the nodata value for the raster
                                          )

    # create a dictionary with the zonal results added to each polygon in the input GeoDataFrame
    polygons_dict = dict()
    for ind, row in starting_polygons.iterrows():
        polygons_dict[row[starting_polygons_identifying_column]] = zonal_stats[ind]

    # use dict and zip with the category names to create a dictionary of category names and the names for the columns
    column_dict = dict(zip(raster_category_map.values(), column_names))

    # add rasterstats results to columns in GeoDataFrame
    for ind, row in starting_polygons.iterrows(): # use iterrows to iterate over each row in the GeoDataFrame
        results_data = polygons_dict[row[starting_polygons_identifying_column]] # get the category data for this polygon
        for category in raster_category_map.values(): # iterate over each of the category class names
            try:
                # add the category count to a new column
                starting_polygons.loc[ind, column_dict[category]] = results_data[category]
            except KeyError:
                # if category name is not present, value should be 0
                starting_polygons.loc[ind, column_dict[category]] = 0

    for ind, row in starting_polygons.iterrows(): # iterate over each row in the GeoDataFrame
        starting_polygons.loc[ind, column_names] = 100 * row[column_names] / row[column_names].sum()
    return starting_polygons # an updated version


def find_stat_values_underlying_raster_for_polygons(starting_polygons, starting_polygons_identifying_column_integer,
                                                   desired_column_names, raster, affine, fill_value=0,
                                                   starting_polygons_geometry_column='geometry'):
    """
    Generate statistics values of underlying raster for each polygon in a GeoDataFrame.

    - Create list of geometry, value pairs for polygons in GeoDataFrame
    - Rasterize vector polygons using rasterio in order to create masks for each polygon extent
    - Calculate statistics of underlying raster using masks by iterating over each polygon in GeoDataFrame
    - Add values to input polygon GeoDataFrame by creating new columns with desired column names

    Parameters
    ----------
    starting_polygons : GeoDataFrame
        Input polygon(s)
    starting_polygons_identifying_column_integer : column label
        Column name that uniquely identifies starting polygons - values in column must be of type integer and must not
        contain values equal to the fill value used with the vector geometries below
    desired_column_names : dict
        Dictionary of desired output column names - dictionary must have keys as below
                    {'mean': 'desired column name',
                    'min': 'desired column name',
                    'max': 'desired column name',
                    'range': 'desired column name',
                    'median': 'desired column name',
                    'std': 'desired column name'}
    raster : ndarray
        Input raser
    affine : Affine instance
        The input raster geotransform
    fill value : int, default 0
        The value to use for areas not covered by the polygon geometries
    starting_polygons_geometry_column : column label, default 'geometry'
        Column name of geometry column in polygon GeoDataFrame

    Returns
    -------
    GeoDataFrame
        Updated GeoDataFrame
    """
    # get a list of geometry, value pairs
    shapes = list(zip(starting_polygons[starting_polygons_geometry_column],
                      starting_polygons[starting_polygons_identifying_column_integer]))

    # create a raster based on the vector polygons
    site_mask = rio.features.rasterize(shapes=shapes,  # the list of geometry/value pairs
                                       fill=fill_value,  # the value to use for cells not covered by any geometry
                                       out_shape=raster.shape,  # the shape of the new raster
                                       transform=affine)  # the geotransform of the new raster

    for ind, row in starting_polygons.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate mean value using mask and add to column
        starting_polygons.loc[ind, desired_column_names['mean']] = np.nanmean(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]])
    for ind, row in starting_polygons.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate min value using mask and add to column
        starting_polygons.loc[ind, desired_column_names['min']] = np.nanmin(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]])
    for ind, row in starting_polygons.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate max value using mask and add to column
        starting_polygons.loc[ind, desired_column_names['max']] = np.nanmax(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]])
    for ind, row in starting_polygons.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate range value by finding max and min values and subtracting them using mask and add to column
        starting_polygons.loc[ind, desired_column_names['range']] = (np.nanmax(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]]) -
                                                          np.nanmin(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]]))
    for ind, row in starting_polygons.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate median value using mask and add to column
        starting_polygons.loc[ind, desired_column_names['median']] = np.nanmedian(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]])
    for ind, row in starting_polygons.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate standard deviation value using mask and add to column
        starting_polygons.loc[ind, desired_column_names['std']] = np.nanstd(raster[site_mask ==
                                                                row[starting_polygons_identifying_column_integer]])
    return starting_polygons # an updated version


# load the input shapefile datasets from the data_files folder using gpd.read_file(os.path.abspath('<file_path>'))
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))

# transform data files to myCRS using gdf.to_crs() (UTM(29) has an epsg of 32629 which will give measurements in metres)
# transform all input data_files to ensure all data is on the same reference system using inplace=true as we want to
# transform the datasets here and not create new ones
sites.to_crs(epsg=32629, inplace=True)
counties.to_crs(epsg=32629, inplace=True)

# find the county each site falls within
# first we need to create a new column with county names that will display better in the output - we want to ensure the
# name is no longer all in capitals and we want to add the word 'County' to the county name
for ind, row in counties.iterrows():  # iterate over each row in the GeoDataFrame
    counties.loc[ind, 'County'] = 'County '+row['CountyName'].title()  # assign the row's CountyName to a new column
                                                              # called County (matching the corresponding csv column)
# find the underlying county and add it's name to a new column
sites = find_underlying_vector_value(sites, 'Name', counties, 'County')

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
    lc_crs = dataset.crs
    landcover = dataset.read(1)
    lc_affine_tfm = dataset.transform

# we need to ensure the vector layer is in the same crs as the raster before performing the next step
sites.to_crs(lc_crs, inplace=True)

# to find the percentage land cover for each site we first need to define the land cover category map dictionary which
# maps the raster values to the land cover categories and a desired column names list:
# category map dictionary
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

# desired column names list
short_names = ['%broadleaf',
               '%coniferous',
               '%arable',
               '%imp_grass',
               '%nat_grass',
               '%mountain',
               '%saltwater',
               '%freshwater',
               '%coastal',
               '%built_up']

# the percentage land cover is the calculated using the find_percentage_underlying_raster_categories_for_polygons()
# function
find_percentage_underlying_raster_categories_for_polygons(sites, 'Name', landcover, landcover_names, short_names,
                                                          lc_affine_tfm)

# open the land cover raster and read the data
with rio.open('data_files/GBR_wind-speed_100m.tif') as dataset:
    ws_crs = dataset.crs
    wind_speed = dataset.read(1)
    ws_affine_tfm = dataset.transform

# we need to ensure the vector layer is in the same crs as the raster before performing the next step
sites.to_crs(ws_crs, inplace=True)

# before running the find_stat_values_underlying_raster_for_polygons function below we need to make sure we have an
# integer column that uniquely identifies each polygon in the GeoDataFrame that does not contain a value that is the
# same as the value we will use as our fill value (0)
# create an integer identifier column for the raster analysis called ID_RA which starts with 1 and increments by 1 for
# each row until the end of the GeoDataFrame
sites['ID_RA'] = range(1, 1+len(sites))

# we also need a dictionary of statists and the desired column names for the statistics
wind_speed_stat_columns_dict = {'mean': 'MeanWindSpeed',
                                'min': 'MinWindSpeed',
                                'max': 'MaxWindSpeed',
                                'range': 'WindSpeedRange',
                                'median': 'MedianWindSpeed',
                                'std': 'WindSpeedStdDev'}

# calculate statistics for the wind speed raster for each site
find_stat_values_underlying_raster_for_polygons(sites, 'ID_RA', wind_speed_stat_columns_dict, wind_speed, ws_affine_tfm,
                                                0, 'geometry')

sites = sites.rename(columns={'Name': 'Site Name'}) # for the sites we will rename the Name column to Site Name so that
                                                    # we can join this data with the proximity data

# create DataFrame to export to csv by copying only desired columns from GeoDataFrame that contains results by dropping
# unnecessary columns (geometry and ID_RA)
site_results = sites.drop(columns=['geometry', 'ID_RA']).copy()

# now that the results are in a DataFrame as opposed to a GeoDataFrame we can round all results to 2 decimal places
site_results = site_results.round(2)

# save the results DataFrame as a csv called Proximity_analysis.csv to the output_files folder with the index removed
site_results.to_csv('output_files/Site_characteristics.csv', index=False)