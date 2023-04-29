import os
import pandas as pd
import geopandas as gpd
from numpy import ceil
from shapely.geometry import Point, LineString, Polygon


def number_objects_within_distance(object_to_buff, distance, objects_to_select):
    """
    Generate count of the number of objects within a specified distance of another specified object(s).

    - Buffer the input object(s) to the specified distance
    - Clip the selecting objects GeoDataFrame to the buffer extent
    - Count the number of selecting objects within the clipped extent

    Parameters
    ----------
    object_to_buff : GeoDataFrame
        Input object(s)
    distance : float
        Desired buffer distance
    objects_to_select : GeoDataFrame
        Selecting objects

    Returns
    -------
    float
        Count value - number of objects
    """
    # create a buffer around the input object(s) at the specified distance
    buffer = object_to_buff.geometry.buffer(distance)
    # clip the selecting objects to the extent of the buffer
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    count = len(clipped.index) # count the number of rows and therefore number of objects within the clipped extent
    return count


def sum_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select, attribute_column):
    """
    Generate sum of a specified attribute of objects within a specified distance of another specified object(s).

    - Buffer the input object(s) to the specified distance
    - Clip the selecting objects GeoDataFrame to the buffer extent
    - Sum the specified attribute within the clipped extent

    Parameters
    ----------
    object_to_buff : GeoDataFrame
        Input object(s)
    distance : float
        Desired buffer distance
    objects_to_select : GeoDataFrame
        Selecting objects
    attribute_column : column label
        Desired attribute to sum

    Returns
    -------
    float
        Sum of specified attribute
    """
    # create a buffer around the input object(s) at the specified distance
    buffer = object_to_buff.geometry.buffer(distance)
    # clip the selecting objects to the extent of the buffer
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    sum_att = clipped[attribute_column].sum() # sum the specified column within the clipped extent GeoDataFrame
    return sum_att


def get_num_unique_values_of_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select,
                                                              attribute_column):
    """
    Generate count of the number of unique values of a specified attribute of objects within a specified distance of
    another specified object(s).

    - Buffer the input object(s) to the specified distance
    - Clip the selecting objects GeoDataFrame to the buffer extent
    - Count the number of unique values of the specified attribute within the clipped extent

    Parameters
    ----------
    object_to_buff : GeoDataFrame
        Input object(s)
    distance : float
        desired buffer distance
    objects_to_select : GeoDataFrame
        Selecting objects
    attribute_column : column label
        Desired attribute to get unique values from

    Returns
    -------
    float
        Count value - number of unique values
    """
    # create a buffer around the input object(s) at the specified distance
    buffer = object_to_buff.geometry.buffer(distance)
    # clip the selecting objects to the extent of the buffer
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    # count the number of the unique values of the specified column within the clipped extent GeoDataFrame
    num_u_values = len(pd.unique(clipped[attribute_column]))
    return num_u_values


def get_unique_values_of_attribute_of_objects_within_distance(object_to_buff, distance, objects_to_select,
                                                              attribute_column):
    """
    Generate string containing the unique values of a specified attribute of objects(s) within a specified distance of
    another specified object(s).

    - Buffer the input object(s) to the specified distance
    - Clip the selecting objects GeoDataFrame to the buffer extent
    - Get list of unique values of specified attribute within the clipped extent
    - Sort alphabetically
    - Convert to improve display and allow output to be saved to GeoDataFrame cell

    Parameters
    ----------
    object_to_buff : GeoDataFrame
        Input object(s)
    distance : float
        desired buffer distance
    objects_to_select : GeoDataFrame
        Selecting objects
    attribute_column : column label
        Desired attribute to get unique values from

    Returns
    -------
    str
        String of unique values separated by commas
    """
    # create a buffer around the input object(s) at the specified distance
    buffer = object_to_buff.geometry.buffer(distance)
    # clip the selecting objects to the extent of the buffer
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    # get a list of the unique values of the specified column within the clipped extent GeoDataFrame which is of type
    # numpy.ndarray object
    u_values_n_array = (pd.unique(clipped[attribute_column]))
    u_values_n_array.sort() # sort the values into alphabetical order
    u_values = ', '.join(u_values_n_array) # turn the array into a string with commas separating unique values
    return u_values


def sum_length_of_lines_within_distance(object_to_buff, distance, objects_to_select):
    """
    Generate sum of the line lengths within a specified distance of another specified object(s).

    - Buffer the input object(s) to the specified distance
    - Clip the selecting objects GeoDataFrame to the buffer extent
    - Calculate lengths of lines within the clip extent
    - Sum the lengths within the clipped extent

    Parameters
    ----------
    object_to_buff : GeoDataFrame
        Input object(s)
    distance : float
        Desired buffer distance
    objects_to_select : GeoDataFrame
        Selecting objects

    Returns
    -------
    float
        Sum of line lengths
    """
    # create a buffer around the input object(s) at the specified distance
    buffer = object_to_buff.geometry.buffer(distance)
    # clip the selecting objects to the extent of the buffer
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    for ind, row in clipped.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate the length of the line and assign it to a new column called length
        clipped.loc[ind, 'Length'] = row['geometry'].length
    sum_len = clipped['Length'].sum() # sum the length column within the clipped extent GeoDataFrame
    return sum_len


def find_nearest_and_distance_to_nearest(starting_objects, starting_objects_identifying_column, objects_to_select,
                                         objects_to_select_identifying_column, distance_column_name='Distances'):
    """
    Update a GeoDataFrame to include two new columns related to the nearest object from another specified GeoDataFrame,
    the distance to the nearest object and an identifier.

    - Spatial join (type = nearest) between input GeoDataFrames
    - Get columns needed for output
    - Merge required output columns back onto the starting object(s) GeoDataFrame

    Parameters
    ----------
    starting_objects : GeoDataFrame
        Input object(s)
    starting_objects_identifying_column : column label
        Column name of attribute that identifies objects in GeoDataFrame
    objects_to_select : GeoDataFrame
        Selecting objects
    objects_to_select_identifying_column : column label
        Column name of attribute that identifies objects in GeoDataFrame
    distance_column_name : str, default 'Distances'
        Desired name of distance attribute column

    Returns
    -------
    GeoDataFrame
        Updated GeoDataFrame
    """
    # create a new GeoDataFrame and for each starting object find the nearest selecting object using
    # gpd.sjoing_nearest() and specifying the distance column (distance_col) to be calculated and stored with the
    # specified distance column name
    joined_distances = gpd.sjoin_nearest(starting_objects, objects_to_select, distance_col=distance_column_name)
    # the joined GeoDatFrame contains all columns from both input GeoDataFrame so create a new GeoDataFrame with only
    # the columns required for the output (the two identifying columns and the distance column)
    results = joined_distances[[starting_objects_identifying_column, distance_column_name,
                             objects_to_select_identifying_column]].copy()
    # merge the required output columns back onto the starting objects GeoDataFrame based on the starting object
    # identifying column
    updated = starting_objects.merge(results, on=[starting_objects_identifying_column])
    return updated


# load the input shapefile datasets from the data_files folder using gpd.read_file(os.path.abspath('<file_path>'))
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))
places = gpd.read_file(os.path.abspath('data_files/NI_Places.shp'))
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))
rivers = gpd.read_file(os.path.abspath('data_files/Rivers.shp'))
lakes = gpd.read_file(os.path.abspath('data_files/Lakes.shp'))
AONBs = gpd.read_file(os.path.abspath('data_files/NI_AONBs.shp'))

# transform data files to myCRS using gdf.to_crs() (UTM(29) has an epsg of 32629 which will give measurements in metres)
# transform all input data_files to ensure all data is on the same reference system using inplace=true as we want to
# transform the datasets here and not create new ones
sites.to_crs(epsg=32629, inplace=True)
places.to_crs(epsg=32629, inplace=True)
counties.to_crs(epsg=32629, inplace=True)
rivers.to_crs(epsg=32629, inplace=True)
lakes.to_crs(epsg=32629, inplace=True)
AONBs.to_crs(epsg=32629, inplace=True)

# to add the additional place information, firstly, we need to create a spatial join between the counties and places
# GeoDataFrames - we do this using gpd.sjoin(), places will be left and counties right as we want to
# add the county information to the places GeoDataFrame
place_and_county = gpd.sjoin(places, counties, how='inner', lsuffix='left', rsuffix='right')

# we then need to create a new column with nice names for counties and nice names for places so that the column names
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

# now for the analysis
# find the number of places within 5km of each site and store results to new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.loc[ind, 'numplaces5km'] = number_objects_within_distance(
                                     row, 5000, places_wi) # assign the count to a new column called numplaces5km

# find the names of places within 5km of each site and store results to new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.at[ind, 'places5km'] = get_unique_values_of_attribute_of_objects_within_distance(
                                 row, 5000, places_wi, 'Name') # assign the names to a new column called places5km

# find the sum of the populations of the places within 5km of each site and store results to new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.loc[ind, 'sumpop5km'] = sum_attribute_of_objects_within_distance(
                                  row, 5000, places_wi, 'Population') # assign the population sum to a new column
                                                                      # called sumpop5km
sites['sumpop5km'] = ceil(sites['sumpop5km']) # remove decimal places from results using the numpy's ceil() function

# find the number of different councils the places within 5km of each site are situated within and store results to
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.at[ind, 'numcouncils5km'] = get_num_unique_values_of_attribute_of_objects_within_distance(
                                      row, 5000, places_wi, 'LGD') # assign the count to a new column called
                                                                   # numcouncils5km

# find the council names the places within 5km of each site are situated within and store results to new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.at[ind, 'councils5km'] = get_unique_values_of_attribute_of_objects_within_distance(
                                   row, 5000, places_wi, 'LGD') # assign the names to a new column called councils5km

# find the lenght of rivers within 5km of each site and store results to new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    sites.loc[ind, 'lengthriver5km'] = sum_length_of_lines_within_distance(
                                       row, 5000, rivers)/1000 # divide the sum by 1000 to convert from metres to
                                       # kilometres and assign the kilometre sum to a new column called places5km
sites['lengthriver5km'] = sites['lengthriver5km'].round(2) # round results to 2 decimal places (rounding to nearest 10m)

# For clarity and to allow the find nearest function to run we need to ensure the identifier (Name) columns within all
# GeoDataFrames are unique

sites = sites.rename(columns={'Name': 'Site Name'}) # for the sites we will reaname the Name column to Site Name

# For the AONBs we want to create a new AONB Name column as we also want to add the word AONB to the AONB Name cells
for ind, row in AONBs.iterrows():  # iterate over each row in the GeoDataFrame
    AONBs.loc[ind, 'AONB Name'] = row['NAME']+' AONB'

# find thh nearest AONB to each site adding distance to nearest AONB and the AONB's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  AONBs, 'AONB Name', 'DistNearestAONB')
sites['DistNearestAONB'] = (sites['DistNearestAONB']/1000).round(2) # round results to 2 decimal places (rounding to
                                                                    # nearest 10m)

# create DataFrame to export to csv by copying only desired columns from GeoDataFrame that contains results by dropping
# unnecessary columns (geometry and Id)
site_results = sites.drop(columns=['geometry', 'Id']).copy()

# save the results DataFrame as a csv called Proximity_analysis.csv to the output_files folder with the index removed
site_results.to_csv('output_files/Proximity_analysis.csv', index=False)