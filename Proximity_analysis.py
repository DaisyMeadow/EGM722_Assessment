import os
import pandas as pd
import geopandas as gpd
from numpy import ceil


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
    count = len(clipped.index)  # count the number of rows and therefore number of objects within the clipped extent
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
    sum_att = clipped[attribute_column].sum()  # sum the specified column within the clipped extent GeoDataFrame
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
    num_u_values = int(num_u_values)  # convert to integer as these are count values and therefore whole numbers
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
    u_values_n_array.sort()  # sort the values into alphabetical order
    u_values = '; '.join(u_values_n_array)  # turn the array into a string with semicolons separating unique values
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
    sum_len = clipped['Length'].sum()  # sum the length column within the clipped extent GeoDataFrame
    return sum_len


def sum_area_of_polygons_within_distance(object_to_buff, distance, objects_to_select):
    """
    Generate sum of polygon areas within a specified distance of another specified object(s).

    - Buffer the input object(s) to the specified distance
    - Clip the selecting objects GeoDataFrame to the buffer extent
    - Calculate areas of polygons within the clip extent
    - Sum the areas within the clipped extent

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
        Sum of areas
    """
    # create a buffer around the input object(s) at the specified distance
    buffer = object_to_buff.geometry.buffer(distance)
    # clip the selecting objects to the extent of the buffer
    clipped = gpd.clip(objects_to_select, buffer, keep_geom_type=True)
    for ind, row in clipped.iterrows():  # iterate over each row in the GeoDataFrame
        # calculate the area of the polygon and assign it to a new column called Area
        clipped.loc[ind, 'Area'] = row['geometry'].area
    sum_area = clipped['Area'].sum()  # sum the area column within the clipped extent GeoDataFrame
    return sum_area


def find_nearest_and_distance_to_nearest(starting_objects, starting_objects_identifying_column, objects_to_select,
                                         objects_to_select_identifying_column,
                                         desired_distance_column_name='Distances'):
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
    desired_distance_column_name : str, default 'Distances'
        Desired name of distance attribute column

    Returns
    -------
    GeoDataFrame
        Updated GeoDataFrame
    """
    # create a new GeoDataFrame and for each starting object find the nearest selecting object using
    # gpd.sjoing_nearest() and specifying the distance column (distance_col) to be calculated and stored with the
    # specified distance column name
    joined_distances = gpd.sjoin_nearest(starting_objects, objects_to_select, distance_col=desired_distance_column_name)
    # the joined GeoDatFrame contains all columns from both input GeoDataFrame so create a new GeoDataFrame with only
    # the columns required for the output (the two identifying columns and the distance column)
    results = joined_distances[[starting_objects_identifying_column, desired_distance_column_name,
                               objects_to_select_identifying_column]].copy()
    # merge the required output columns back onto the starting objects GeoDataFrame based on the starting object
    # identifying column
    updated = starting_objects.merge(results, on=[starting_objects_identifying_column])
    return updated


def add_characters_to_end_of_object_names(objects, original_naming_column, characters_to_add,
                                          desired_new_naming_column_name='New Name'):
    """
    Update a naming column within a GeoDataFrame by adding specified characters to the end of the original object names
    saving the new name to a new column with a specified column label.

    - Iterate over each object in the GeoDataFrames
    - Add the desired characters to the end of the object's name creating a new name
    - Add the new name to the new names column

    Parameters
    ----------
    objects : GeoDataFrame
        Input objects
    original_naming_column : column label
        Column name of naming column for objects in GeoDataFrame
    characters_to_add : str
        Desired characters to add to the end of the object names
    desired_new_naming_column_name : str, default 'New Name'
        Desired name of new names column

    Returns
    -------
    GeoDataFrame
        Updated GeoDataFrame
    """
    for ind, row in objects.iterrows():  # iterate over each row in the GeoDataFrame
        # add the desired characters to the old object names and assign to a new column with the desired new naming
        # column label
        objects.loc[ind, desired_new_naming_column_name] = row[original_naming_column] + characters_to_add
    return objects  # an updated version


# load the input shapefile datasets from the data_files folder using gpd.read_file(os.path.abspath())
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))
places = gpd.read_file(os.path.abspath('data_files/NI_Places.shp'))
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))
rivers = gpd.read_file(os.path.abspath('data_files/Rivers.shp'))
ACRAs = gpd.read_file(os.path.abspath('data_files/NI_Agricultural_Critical_Risk_Areas.shp'))
peatland = gpd.read_file(os.path.abspath('data_files/NI_Peatland.shp'))
nat_reserves = gpd.read_file(os.path.abspath('data_files/NI_NNR_and_NR.shp'))
parks_gardens = gpd.read_file(os.path.abspath('data_files/NI_Parks_Gardens_Dec2022.shp'))
ramsars = gpd.read_file(os.path.abspath('data_files/NI_Ramsar_Sites.shp'))
AONBs = gpd.read_file(os.path.abspath('data_files/NI_AONBs.shp'))
ASSIs = gpd.read_file(os.path.abspath('data_files/NI_ASSIs.shp'))
SACs = gpd.read_file(os.path.abspath('data_files/NI_SACs.shp'))
SPAs = gpd.read_file(os.path.abspath('data_files/NI_SPAs.shp'))

# transform data files to Northern Ireland (NI) Universal Transverse Mercator zone (UTM(29) which has an epsg of 32629)
# which will give measurements in metres using gdf.to_crs()
# transform all input data_files to ensure all data is on the same reference system using inplace=true as we want to
# transform the datasets here and not create new ones
sites.to_crs(epsg=32629, inplace=True)
places.to_crs(epsg=32629, inplace=True)
counties.to_crs(epsg=32629, inplace=True)
rivers.to_crs(epsg=32629, inplace=True)
ACRAs.to_crs(epsg=32629, inplace=True)
peatland.to_crs(epsg=32629, inplace=True)
nat_reserves.to_crs(epsg=32629, inplace=True)
parks_gardens.to_crs(epsg=32629, inplace=True)
ramsars.to_crs(epsg=32629, inplace=True)
AONBs.to_crs(epsg=32629, inplace=True)
ASSIs.to_crs(epsg=32629, inplace=True)
SACs.to_crs(epsg=32629, inplace=True)
SPAs.to_crs(epsg=32629, inplace=True)

# to add the additional place information, firstly, we need to create a spatial join between the counties and places
# GeoDataFrames - we do this using gpd.sjoin(), places will be left and counties right as we want to
# add the county information to the places GeoDataFrame
place_and_county = gpd.sjoin(places, counties, how='inner', lsuffix='left', rsuffix='right')

# the place_and_county county names and place names are all capitals and will therefore not match the names in the place
# information csv so we need to create a new column with nice names for counties and nice names for places and assign
# them to column names that will match those we wish to merge with in the place information csv file
for ind, row in place_and_county.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the row's CountyName to a new column called County (matching the corresponding csv column)
    place_and_county.loc[ind, 'County'] = row['CountyName'].title()
for ind, row in place_and_county.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the row's PLACENAME to a new column called Name (matching the corresponding csv column)
    place_and_county.loc[ind, 'Name'] = row['PLACENAME'].title()

# load the input csv Place_information file from the data_files folder using pd.read_csv()
place_info = pd.read_csv('data_files/Place_information.csv')

# join the place_and_county GeoDataFrame with the place_info csv data by merging on the shared variables/columns
# (place and county) using gdf.merge()
places_wi = place_and_county.merge(place_info, on=["Name", "County"])
# we now have an updated places GeoDataFrame with info - places_wi

# now for the analysis
# first we need to define our search buffer distance
my_buffer_dist = 5000  # here we will use a buffer of 5km (we use 5000 here as our map units are in metres)
my_buff_column = '5km'  # here we will define our distance to be added to end of our column labels

# find the number of places within our search buffer distance of each site using number_objects_within_distance()
# function previously defined and assign results to a new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the count to a new column called 'NumPlaces(our buffer distance)'
    sites.loc[ind, ('NumPlaces' + my_buff_column)] = number_objects_within_distance(row, my_buffer_dist, places_wi)

# find the names of places within our search buffer distance of each site using
# get_unique_values_of_attribute_of_objects_within_distance() function previously defined and assign results to a new
# column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the names to a new column called 'Places(our buffer distance)'
    sites.at[ind, ('Places' + my_buff_column)] = get_unique_values_of_attribute_of_objects_within_distance(
                                                 row, my_buffer_dist, places_wi, 'Name')

# find the sum of the populations of the places within our search buffer distance of each site using
# sum_attribute_of_objects_within_distance() function previously defined and assign result to a new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the population sum to a new column called 'SumPop(our buffer distance)'
    sites.loc[ind, ('SumPop' + my_buff_column)] = sum_attribute_of_objects_within_distance(
                                                  row, my_buffer_dist, places_wi, 'Population')
# remove decimal places from results using the numpy's ceil() function
sites[('SumPop' + my_buff_column)] = ceil(sites[('SumPop' + my_buff_column)])

# find the number of different LGDs the places within our search buffer distance of each site are situated within
# using get_num_unique_values_of_attribute_of_objects_within_distance() function previously defined and assign results
# to a new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the count to a new column called 'NumPlaceLGDs(our buffer distance)'
    sites.at[ind, ('NumPlaceLGDs' + my_buff_column)] = get_num_unique_values_of_attribute_of_objects_within_distance(
                                                       row, my_buffer_dist, places_wi, 'LGD')

# find the LGD names the places within our search buffer distance of each site are situated within using
# get_unique_values_of_attribute_of_objects_within_distance() function previously defined and assign results to a new
# column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the names to a new column called 'PlaceLGDs(our buffer distance)'
    sites.at[ind, ('PlaceLGDs' + my_buff_column)] = get_unique_values_of_attribute_of_objects_within_distance(
                                                    row, my_buffer_dist, places_wi, 'LGD')

# find the length of rivers within our search buffer distance of each site using sum_length_of_lines_within_distance()
# function previously defined and assign results to a new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # divide the sum by 1000 to convert from metres to kilometres and assign the kilometre sum to a new column called
    # 'LengthRivers(our buffer distance)'
    sites.loc[ind, ('LengthRivers' + my_buff_column)] = sum_length_of_lines_within_distance(
                                                        row, my_buffer_dist, rivers)/1000
# round results to 2 decimal places (rounding to nearest 10m)
sites[('LengthRivers' + my_buff_column)] = sites[('LengthRivers' + my_buff_column)].round(2)

# find the total area of the Agricultural Critical Risk Areas (ACRAs) within our search buffer distance of each site
# using sum_area_of_polygons_within_distance() function previously defined and assign result to a new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # divide the area sum by 1000000 to convert from metres squared to kilometres squared and assign the kilometre area
    # sum to a new column called 'Areakm2ACRAs(our buffer distance)'
    sites.loc[ind, ('Areakm2ACRAs' + my_buff_column)] = sum_area_of_polygons_within_distance(
                                                        row, my_buffer_dist, ACRAs)/1000000
# round results to 2 decimal places
sites[('Areakm2ACRAs' + my_buff_column)] = sites[('Areakm2ACRAs' + my_buff_column)].round(2)

# find the type and total area of Priority Habitat Peatland within our search buffer distance of each site
# find the total area of Priority Habitat Peatland within our search buffer distance of each site within 5km of each
# site using sum_area_of_polygons_within_distance() function previously defined and assign result to a new column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # divide the area sum by 1000000 to convert from metres squared to kilometres squared and assign the kilometre area
    # sum to a new column called 'Areakm2Peatland(our buffer distance)'
    sites.loc[ind, ('Areakm2Peatland' + my_buff_column)] = sum_area_of_polygons_within_distance(
                                                           row, my_buffer_dist, peatland)/1000000
# round results to 2 decimal places
sites[('Areakm2Peatland' + my_buff_column)] = sites[('Areakm2Peatland' + my_buff_column)].round(2)

# find the type of Priority Habitat Peatland within our search buffer distance of each site using
# get_unique_values_of_attribute_of_objects_within_distance() function previously defined and assign results to a new
# column
for ind, row in sites.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the names to a new column called 'PeatlandType(our buffer distance)'
    sites.at[ind, ('PeatlandType' + my_buff_column)] = get_unique_values_of_attribute_of_objects_within_distance(
                                                       row, my_buffer_dist, peatland, 'DESCRIPTIO')

# for clarity and to allow the find_nearest_and_distance_to_nearest() function previously defined to run we need to
# ensure the identifier columns within all GeoDataFrames are unique (for example, we cannot pass two identifier
# columns called 'Name')
# for the sites we will rename the Name column to Site Name
sites = sites.rename(columns={'Name': 'Site Name'})

# for the Nature Reserves and National Nature Reserves we will rename the Name column to Reserve Name
nat_reserves = nat_reserves.rename(columns={'NAME': 'Reserve Name'})

# for the Historic Parks and Gardens we want to create names that are not in all capitals
for ind, row in parks_gardens.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the row's SITE (naming column) to a new column called Park/Garden Name
    parks_gardens.loc[ind, 'Park/Garden Name'] = row['SITE'].title()

# for the Ramsar Sites we want to create names that are not in all capitals
for ind, row in ramsars.iterrows():  # iterate over each row in the GeoDataFrame
    # assign the row's NAME to a new column called Ramsar Site Name
    ramsars.loc[ind, 'Ramsar Site Name'] = row['NAME'].title()

# for the remaining GeoDataFrame we create a new name column whilst also adding some character's to the object names
# using the add_characters_to_end_of_object_names() function previously defined
# for the Areas of Outstanding Natural Beauty (AONBs) we want to add the word AONB and call the new name column AONB
# Name
add_characters_to_end_of_object_names(AONBs, 'NAME', ' AONB', 'AONB Name')

# for the Areas of Special Scientific Interest (ASSIs) we want to add the word ASSI and call the new name column ASSI
# Name
add_characters_to_end_of_object_names(ASSIs, 'NAME', ' ASSI', 'ASSI Name')

# for the Special Areas of Conservation (SACs) we want to add the word SAC and call the new name column SAC Name
add_characters_to_end_of_object_names(SACs, 'NAME', ' SAC', 'SAC Name')

# for the Special Protection Areas (SPAs) we want to add the word SPA and call the new name column SPA Name
add_characters_to_end_of_object_names(SPAs, 'NAME', ' SPA', 'SPA Name')

# find the nearest nature reserve to each site using find_nearest_and_distance_to_nearest() function previously defined
# assigning distance to nearest nature reserve and the reserve's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  nat_reserves, 'Reserve Name',
                                             'DistNearestNatureReserve')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestNatureReserve'] = (sites['DistNearestNatureReserve']/1000).round(2)

# find the nearest park/garden to each site using find_nearest_and_distance_to_nearest() function previously defined
# assigning distance to nearest park/garden and the park/garden's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  parks_gardens, 'Park/Garden Name',
                                             'DistNearestParkGarden')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestParkGarden'] = (sites['DistNearestParkGarden']/1000).round(2)

# find the nearest ramsar site to each site using find_nearest_and_distance_to_nearest() function previously defined
# assigning distance to nearest ramsar site and the ramsar site's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  ramsars, 'Ramsar Site Name',
                                             'DistNearestRamsarSite')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestRamsarSite'] = (sites['DistNearestRamsarSite']/1000).round(2)

# find the nearest AONB to each site using find_nearest_and_distance_to_nearest() function previously defined assigning
# distance to nearest AONB and the AONB's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  AONBs, 'AONB Name', 'DistNearestAONB')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestAONB'] = (sites['DistNearestAONB']/1000).round(2)

# find the nearest ASSI to each site using find_nearest_and_distance_to_nearest() function previously defined assigning
# distance to nearest ASSI and the ASSI's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  ASSIs, 'ASSI Name', 'DistNearestASSI')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestASSI'] = (sites['DistNearestASSI']/1000).round(2)

# find the nearest SAC to each site using find_nearest_and_distance_to_nearest() function previously defined assigning
# distance to nearest SAC and the SAC's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  SACs, 'SAC Name', 'DistNearestSAC')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestSAC'] = (sites['DistNearestSAC']/1000).round(2)

# find the nearest SPA to each site using find_nearest_and_distance_to_nearest() function previously defined assigning
# distance to nearest SPA and the SPA's name to new columns
sites = find_nearest_and_distance_to_nearest(sites, 'Site Name',  SPAs, 'SPA Name', 'DistNearestSPA')
# convert the distance from metres to kilometres by diving by 1000 and round results to 2 decimal places (rounding to
# nearest 10m)
sites['DistNearestSPA'] = (sites['DistNearestSPA']/1000).round(2)

# create DataFrame to export to csv by converting the sites GeoDataFrame to a DataFrame by dropping the geometry column
site_results = sites.drop(columns=['geometry']).copy()

# save the results DataFrame as a csv called Proximity_analysis.csv to the output_files folder with the index removed
site_results.to_csv('output_files/Proximity_analysis.csv', index=False)