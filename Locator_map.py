import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from shapely.ops import unary_union

# function to generate matplotlib handles for polygon features on the map to add to a legend
def generate_handles(labels, colors, edge='k', linewidth=1, alpha=1):
    '''
    Generates matplotlib handles for polygons to add to a legend.
    
    - Get length of input lists
    - Create empty list of handles
    - Iterate through lists creating a handle for each item and adding it to the handles list

    inputs: labels = list of lables
            colors = list of colors
            edge = edge color to use default='k'
            linewidth = linewidth to use default=1
            alpha = alpha (transparency) to use default=1 (opaque)
    
    returns: list of matplotlib handles
    '''
    lc = len(colors) # get the length of the color list
    handles = [] # create an empty handles list
    for i in range(len(labels)): # iterate through the input lists to access each item with the help of its index
        handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=colors[i % lc], edgecolor=edge,
                                          linewidth=linewidth, alpha=alpha))
    return handles

# function to generate matplotlib handles for polygon features with hatches on the map to add to a legend
def generate_handles_with_hatch(labels, colors, hatch, edge='k', linewidth=1, alpha=1):
    '''
    Generates matplotlib handles for polygons with hatching to add to a legend.

    - Get length of input lists
    - Create empty list of handles
    - Iterate through lists creating a handle for each item and adding it to the handles list

    inputs: labels = list of lables
            colors = list of colors
            hatch = hatch type
            edge = edge color to use default='k'
            linewidth = linewidth to use default=1
            alpha = alpha (transparency) to use default=1 (opaque)

    returns: list of matplotlib handles
    '''
    lc = len(colors) # get the length of the color list
    handles = [] # create an empty handles list
    for i in range(len(labels)): # iterate through the input lists to access each item with the help of its index
        handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=colors[i % lc], hatch=hatch, edgecolor=edge,
                                          linewidth=linewidth, alpha=alpha))
    return handles

# function to generate matplotlib handles for point features with hatches on the map to add to a legend
def generate_handles_points(labels, markers, colors, marker_sizes):
    '''
    Generates matplotlib handles for points to add to a legend.

    - Get length of input lists
    - Create empty list of handles
    - Iterate through lists creating a handle for each item and adding it to the handles list creating list of lists
    - Flatten list of lists to single list of handles that can be used in legend

    inputs: labels = list of lables
            markers = list of marker types
            colors = list of marker colors
            marker_sizes = list of marker sizes

    returns: list of matplotlib handles
    '''
    lc = len(colors)  # get the length of the color list
    handles = [] # create an empty handles list
    for i in range(len(labels)): # iterate through the input lists to access each item with the help of its index
        handles.append(ax.plot([], [], markers[i % lc], color=colors[i % lc], ms=marker_sizes[i % lc]))
    # this gives a list of handle lists for the point types but we need to flatten it to get a single list
    # that can be used to add the handles to the legend as legend handles themselves cannot be lists
    handles = [item for sublist in handles for item in sublist] # flatten the list
    return handles

# function to create a scale bar of length 20 km in the upper right corner of the map
def scale_bar(ax, location=(0.92, 0.95)):
    '''
    Create and plot a scale bar of length 20km at a specified location

    - Turn specified location into coordinates in metres
    - Plot scale bar lines at specified location
    - Add scale bar text labels below scale bar lines

    inputs: ax = plotting axes
            location = desired scale bar location

    returns: scale bar of 20km length plotted at specified location
    '''
    x0, x1, y0, y1 = ax.get_extent() # get extent of plotted area (ax)
    # turn specified scale bar location into coordinates in metres giving scale bar x (sbx) and y (sby)
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    # plot a line of 20000 metres length (20km) - the length of the scale bar from the scale bar x coordinate
    # creating a base scale bar line that is black
    ax.plot([sbx, sbx - 20000], [sby, sby], color='k', linewidth=9, transform=ax.projection)
    # plot a thinner line from the end of the scale bar line to midway (10000 metres (10km)) and color it black
    ax.plot([sbx, sbx - 10000], [sby, sby], color='k', linewidth=6, transform=ax.projection)
    # plot a thinner line from midway (10000 metres (10km)) to the other end of the scale bar line (20000 metres (20km))
    # and color it white
    ax.plot([sbx-10000, sbx - 20000], [sby, sby], color='w', linewidth=6, transform=ax.projection)

    # add text labels of 20km, 10km and 0km
    ## add a text label of 20km to the end of the scale bar (at the scale bar x coordinate)
    ax.text(sbx, sby-4500, '20 km', transform=ax.projection, fontsize=8)
    # add a text label of 10km to the scale bar midway (with an offset to ensure label is centralised to midway)
    ax.text(sbx-12500, sby-4500, '10 km', transform=ax.projection, fontsize=8)
    # add a text label of 0km to the other end of the scale bar (with an offset to ensure label is centralised to the
    # end of the scale bar)
    ax.text(sbx-24500, sby-4500, '0 km', transform=ax.projection, fontsize=8)

# load the input shapefile datasets from the data_files folder using gpd.read_file(os.path.abspath('<file_path>'))
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))
lakes = gpd.read_file(os.path.abspath('data_files/Lakes.shp'))
places = gpd.read_file(os.path.abspath('data_files/NI_Places.shp'))
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))

# create a figure of size 14x11 (representing the page size in inches)
myFig = plt.figure(figsize=(14, 11))

# create a Universal Transverse Mercator reference system to transform our data
myCRS = ccrs.UTM(29)  # Northern Ireland (NI) is in UTM Zone 29, so we pass 29 to ccrs.UTM()

# transform data files to myCRS using gdf.to_crs() (UTM(29) has an epsg of 32629)
# transform all input data_files to ensure all data is on the same reference system using inplace=true as we want to
# transform the datasets here and not create new ones
counties.to_crs(epsg=32629, inplace=True)
lakes.to_crs(epsg=32629, inplace=True)
places.to_crs(epsg=32629, inplace=True)
sites.to_crs(epsg=32629, inplace=True)

# create a NI_outline
NI_Union = unary_union(counties.geometry) # create an outline by joining the geometries of the counties
NI_Outline = gpd.GeoDataFrame(geometry=gpd.GeoSeries(NI_Union)) # create GeoDataFrame based on GeoSeries output of union
NI_Outline.set_crs(epsg=32629, inplace=True) # set CRS of new outline GeoDataFrame to prevent naive geometry

# create an axes object in the figure, using a UTM projection where we can actually plot our data
ax = plt.axes(projection=myCRS)

# first, add the outline of NI using cartopy's ShapelyFeature()
outline_feature = ShapelyFeature(NI_Outline['geometry'], # first argument is the geometry
                                 myCRS, # second argument is the CRS
                                 edgecolor='k', # outline the feature in black
                                 facecolor='w') # set the face color to white
ax.add_feature(outline_feature)  # add the features we've created to the map using ax.add_feature

xmin, ymin, xmax, ymax = NI_Outline.total_bounds # find total bounds of NI outline (we will use this to zoom the map)

# using the boundary of NI outline, zoom the map to our area of interest using ax.set_extent()
ax.set_extent([xmin-17000, xmax+5000, ymin-5000, ymax+5000], crs=myCRS) # because total_bounds gives output as
# xmin, ymin, xmax, ymax, but set_extent takes xmin, xmax, ymin, ymax, we re-order the coordinates here
# additionally, due to the length of the legend we want to have an additional gap to the left of the map so we add a
# larger number of 17000 to xmin to create space for legend

# get a list of unique names for the county boundaries
county_names = list(counties.CountyName.unique())
county_names.sort()  # sort the counties alphabetically by name

# pick colors to use to display the county boundaries creating a list
county_colors = ['thistle', 'palegreen', 'paleturquoise', 'lightcoral', 'lemonchiffon', 'navajowhite']

# add the county outlines to the map using the colors that we've picked by iterating over the unique values in the
# county_names list/'CountyName' field and using cartopy's ShapelyFeature()
for ii, name in enumerate(county_names):
    feat = ShapelyFeature(counties.loc[counties['CountyName'] == name, 'geometry'],  # first argument is the geometry
                          myCRS,  # second argument is the CRS
                          edgecolor='k',  # outline the feature in black
                          facecolor=county_colors[ii],  # set the face color to the corresponding color from the list
                          linewidth=1,  # set the outline width to be 1 pt
                          alpha=0.75)  # set the alpha (transparency) to be 0.75 (out of 1)
    ax.add_feature(feat)  # once we have created the feature, we have to add it to the map using ax.add_feature()

# clip lakes to extent of Northern Ireland using gpd.clip()
NI_lakes = gpd.clip(lakes, NI_Outline, keep_geom_type=True)

# add the lakes to the map using cartopy's ShapelyFeature()
lakes_feat = ShapelyFeature(NI_lakes['geometry'],  # first argument is the geometry
                            myCRS,  # second argument is the CRS
                            edgecolor='cornflowerblue',  # set the edgecolor to be cornflowerblue
                            facecolor='cornflowerblue',  # set the facecolor to be cornflowerblue
                            linewidth=1)  # set the outline width to be 1 pt
ax.add_feature(lakes_feat)  # add the collection of features to the map using ax.add_feature()

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

# create a list of place types - we want these to be displayed in a specific order (size) so we will specify the order
# of the list here instead of getting unique values and sorting them
place_types = ['City', 'Town', 'Suburb', 'Large Village', 'Village', 'Hamlet', 'Townland', 'Location']

# pick colors to use to display the different place types creating a list
place_colors = ['navy', 'brown', 'limegreen', 'darkviolet', 'mediumorchid', 'palevioletred',
                'dimgrey', 'grey'] # create a list to use when plotting based on place_types order
# pick marker types to use to display the different place types creating a list
place_marker = ['o', '^', '^', 'd', 'd', '*', 'o', '*'] # create a list to use when plotting based on place_types order
# pick marker sizes to use to display the different place types creating a list
place_marker_size = [12, 9, 8, 8, 6, 7, 6, 6] # create a list to use when plotting based on place_types order

# add the place markers to the map using the marker types, colors and marker sizes that we've picked by iterating over
# the unique values in the place_types list/'Type' field and using ax.plot() (cartopy's ShapelyFeature() creates polygons
# so we use ax.plot() instead here to plot the point data)
for ii, type in enumerate(place_types):
    plot = places_wi.loc[places_wi['Type'] == type]
    ax.plot(plot.geometry.x, plot.geometry.y, # use geometry x and y of the point
            place_marker[ii], # set the marker type to the corresponding marker type from the list
            color=place_colors[ii], # set the marker color to the corresponding color from the list
            ms=place_marker_size[ii], # set the marker size to the corresponding marker size from the list
            transform=myCRS) # set the crs

# get a list of unique names for the sites
site_names = list(sites.Name.unique())
site_names.sort()  # sort the site names alphabetically by name

# pick colors to use to display the sites creating a list
site_colors = ['orangered', 'olivedrab', 'dodgerblue']

# add the sites to the map using the colors that we've picked by iterating over the unique values in the site_names
# list/'Name' field and using cartopy's ShapelyFeature
for ii, name in enumerate(site_names):
    feat = ShapelyFeature(sites.loc[sites['Name'] == name, 'geometry'],  # first argument is the geometry
                          myCRS,  # second argument is the CRS
                          edgecolor='k',  # outline the feature in black
                          facecolor=site_colors[ii],  # set the face color to the corresponding color from the list
                          hatch='++', # add a hatching to features
                          linewidth=1.5) # set the outline width to be 1.5 pt
    ax.add_feature(feat)  # once we have created the feature, we have to add it to the map using ax.add_feature()

# generate handles for legend
# generate a list of handles for the county datasets using the generate_handle() function created at the start of this
# script
county_handles = generate_handles(county_names, county_colors, linewidth=1, alpha=0.75)

# generate a handle for the lakes datasets using the generate_handles() function created at the start of this script
water_handle = generate_handles(['Lakes'], ['cornflowerblue'], edge='cornflowerblue')

# generate handles for the places datasets using the generate_handles_points() function created at the start of this
# script
place_type_handles = generate_handles_points(place_types, place_marker, place_colors, place_marker_size)

# generate a list of handles for the sites datasets using the generate_handles_with_hatch() function created at the
# start of this script
site_handles = generate_handles_with_hatch(site_names, site_colors, hatch='++', linewidth=1.5)

# update county_names to take it out of uppercase text
nice_names = ['County ' + name.title() for name in county_names]

# create a new list of place types for the legend
place_types_u = place_types
place_types_u[-1] = 'Other' # change the type 'Location' which is the last type in the list (hence index -1) to 'Other'

# ax.legend() takes a list of handles and a list of labels corresponding to the objects you want to add to the legend
handles = site_handles + county_handles + place_type_handles + water_handle # handles list
labels = site_names + nice_names + place_types_u + ['Lakes'] # labels list

leg = ax.legend(handles, # add the list of handles
                labels, # add the list of labels
                title='Map Legend', # add a title to the legend
                title_fontsize=12, # set the font size of the title
                fontsize=10, # set the font size of the labels
                loc='upper left', # set location of the legend on the figure
                frameon=True, # add a border to the legend (True adds border and False removes it)
                framealpha=1) # Set frame transparency (alpha 1 is fully opaque)

# add girdlines to the fulre using ax.gridlines()
gridlines = ax.gridlines(draw_labels=True,  # draw labels for the grid lines
                         xlocs=[-8, -7.5, -7, -6.5, -6, -5.5],  # add longitude lines at 0.5 deg intervals
                         ylocs=[54, 54.5, 55, 55.5])  # add latitude lines at 0.5 deg intervals
gridlines.right_labels = False  # turn off the right-side labels
gridlines.top_labels = False  # turn off the top labels

# add the scale bar to the axis
scale_bar(ax)

# save the figure as Locator_map.png tpo the output_files folder, cropped to the axis (bbox_inches='tight'), and a dpi of 300
myFig.savefig('output_files/Locator_map.png', bbox_inches='tight', dpi=300)