import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, Polygon

# function to buffer around polygon then clip gdf to buffer extent:
    # count number of objects in new gdf
    # sum attribute of objects in new gdf
    # get unique values of attribute in new gdf
    # calculate new lengths within new gdf and sum lengths

def number_objects_within_distance(object_to_buff, distance, objects_to_select):
    buffer = object_to_buff.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    count = len(clipped.index)
    return count

def sum_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select, attribute_column):
    buffer = object_to_buff.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    sum_att = clipped[attribute_column].sum()
    return sum_att

def get_unique_values_of_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select,
                                                              attribute_column):
    buffer = object_to_buff.buffer(distance)
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    u_values = list((pd.unique(clipped[attribute_column])))
    return u_values

def sum_length_of_lines_within_distance(object_to_buff, distance, objects_to_select):
    buffer = object_to_buff.buffer(distance)
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

site_test = sites.iloc[[1]] # get test site gdf

num_test = number_objects_within_distance(site_test, 5000, places)
print(num_test)

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

pop_test = sum_attribute_of_objects_within_distance(site_test, 5000, places_wi, 'Population')
print(pop_test)

LGD_test = get_unique_values_of_attribute_of_objects_within_distance(site_test, 5000, places_wi, 'LGD')
print(LGD_test)

length_test = sum_length_of_lines_within_distance(site_test, 5000, rivers)
print(length_test)

# iterate over sites in gdf - run functions and save outputs to new columns
# when iterating for my purposes, make population count have no decimals when assigning to column and convert length
# to km by dividing by 1000 make it have 2 decimal places when assigning to column


# export results?
# df1 = pd.DataFrame(gdf, copy=True)
# df1 = pd.DataFrame(gdf.drop(columns='geometry'))
# pd.DataFrame.to_csv




