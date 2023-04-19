import os
import pandas as pd
import geopandas as gpd
from numpy import ceil
from shapely.geometry import Point, LineString, Polygon

def number_objects_within_distance(object_to_buff, distance, objects_to_select):
    buffer = object_to_buff.geometry.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    count = len(clipped.index)
    return count

def sum_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select, attribute_column):
    buffer = object_to_buff.geometry.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    sum_att = clipped[attribute_column].sum()
    return sum_att

def get_num_unique_values_of_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select,
                                                              attribute_column):
    buffer = object_to_buff.geometry.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    num_u_values = len(pd.unique(clipped[attribute_column]))
    return num_u_values

def get_unique_values_of_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select,
                                                              attribute_column):
    buffer = object_to_buff.geometry.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    u_values_n_array = (pd.unique(clipped[attribute_column])) # numpy.ndarray object
    u_values_n_array.sort()
    u_values = ', '.join(u_values_n_array) # turn into string with commas separating unique values
    return u_values

def sum_length_of_lines_within_distance(object_to_buff, distance, objects_to_select):
    buffer = object_to_buff.geometry.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    for ind, row in clipped.iterrows():  # iterate over each row in the GeoDataFrame
        clipped.loc[ind, 'Length'] = row['geometry'].length
    sum_len = clipped['Length'].sum()
    return sum_len

# function to find nearest:
    # function to find nearest
    # and then find distance of the object to polygon


# load the input shapefile datasets from the data_files folder using gpd.read_file(os.path.abspath('<file_path>'))
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))
places = gpd.read_file(os.path.abspath('data_files/NI_Places.shp'))
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))
rivers = gpd.read_file(os.path.abspath('data_files/Rivers.shp'))
lakes = gpd.read_file(os.path.abspath('data_files/Lakes.shp'))

# transform data files to myCRS using gdf.to_crs() (UTM(29) has an epsg of 32629 which will give measurements in metres)
# transform all input data_files to ensure all data is on the same reference system using inplace=true as we want to
# transform the datasets here and not create new ones
sites.to_crs(epsg=32629, inplace=True)
places.to_crs(epsg=32629, inplace=True)
counties.to_crs(epsg=32629, inplace=True)
rivers.to_crs(epsg=32629, inplace=True)
lakes.to_crs(epsg=32629, inplace=True)

# note this section is also in Locator_map.py
# to add the additional place information, firstly, we need to create a spatial join between the county and place
# GeoDataFrames - we do this using gpd.sjoin(), places will be left and counties right as we want to add the county
# information to the places GeoDataFrame
place_and_county = gpd.sjoin(places, counties, how='inner', lsuffix='left', rsuffix='right')

# we then need to create a new column with nice names for county and nice names for places so that the column names
# will match those we wish to merge with in the place information csv file as the place_and_county county names and
# place names are all capitals
for ind, row in place_and_county.iterrows():  # iterate over each row in the GeoDataFrame
    place_and_county.loc[ind, 'County'] = row['CountyName'].title()  # assign the row's CountyName to a new column
                                                           # called County (matching the corresponding csv column)
for ind, row in place_and_county.iterrows():  # iterate over each row in the GeoDataFrame
    place_and_county.loc[ind, 'Name'] = row['PLACENAME'].title()  # assign the row's PLACENAME to a new column
                                                           # called Name (matching the corresponding csv column)

# load the input csv Place_information file from the data_files folder using pd.read_csv('<file_path>')
place_info = pd.read_csv('data_files/Place_information.csv')

# join the place_and_county GeoDataFrame with the Place_information csv by merging on the shared variables/columns
# (place and county) using gdf.merge()
places_wi = place_and_county.merge(place_info, on=["Name", "County"])
# we now have an updated places GeoDataFrame with info - places_wi

# iterate over sites in gdf - run functions and save outputs to new columns
# when iterating for my purposes, make population count have no decimals when assigning to column and convert length
# to km by dividing by 1000 make it have 2 decimal places when assigning to column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.loc[ind, 'numplaces5km'] = number_objects_within_distance(row, 5000, places_wi) # assign the .. to a new
                                                                                          # column ..

for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.at[ind, 'places5km'] = get_unique_values_of_attribute_of_objects_within_distance(row, 5000, places_wi,
                                                                                             'Name')

for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.loc[ind, 'sumpop5km'] = sum_attribute_of_objects_within_distance(row, 5000, places_wi, 'Population')
sites['sumpop5km'] = ceil(sites['sumpop5km']) # remove decimal places using the numpy's ceil() function

for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.at[ind, 'numcouncils5km'] = get_num_unique_values_of_attribute_of_objects_within_distance(row, 5000, places_wi,
                                                                                             'LGD')

for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.at[ind, 'councils5km'] = get_unique_values_of_attribute_of_objects_within_distance(row, 5000, places_wi,
                                                                                             'LGD')

for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.loc[ind, 'lengthriver5km'] = sum_length_of_lines_within_distance(row, 5000, rivers)/1000
sites['lengthriver5km'] = sites['lengthriver5km'].round(2) # round to 2 decimal places (rounding to nearest 10m)

# create gdf to export to csv by copying only desired columns from gdf that contains results by dropping unnecessary
# columns
site_results = sites.drop(columns=['geometry', 'Id']).copy()

site_results.to_csv('output_files/Proximity_analysis.csv', index=False)