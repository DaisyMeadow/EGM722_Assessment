import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
from shapely.ops import unary_union

# generate matplotlib handles to create a legend of the features we put in our map.
def generate_handles(labels, colors, edge='k', alpha=1):
    lc = len(colors)  # get the length of the color list
    handles = []
    for i in range(len(labels)):
        handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=colors[i % lc], edgecolor=edge, alpha=alpha))
    return handles

def generate_handles_points(labels, markers, colors, marker_sizes):
    lc = len(colors)  # get the length of the color list
    handles = []
    for i in range(len(labels)):
        handles.append(ax.plot([], [], markers[i % lc], color=colors[i % lc], ms=marker_sizes[i % lc]))
    return handles

# create a scale bar of length 20 km in the upper right corner of the map
# adapted this question: https://stackoverflow.com/q/32333870
# answered by SO user Siyh: https://stackoverflow.com/a/35705477
def scale_bar(ax, location=(0.92, 0.95)):
    x0, x1, y0, y1 = ax.get_extent()
    sbx = x0 + (x1 - x0) * location[0]
    sby = y0 + (y1 - y0) * location[1]

    ax.plot([sbx, sbx - 20000], [sby, sby], color='k', linewidth=9, transform=ax.projection)
    ax.plot([sbx, sbx - 10000], [sby, sby], color='k', linewidth=6, transform=ax.projection)
    ax.plot([sbx-10000, sbx - 20000], [sby, sby], color='w', linewidth=6, transform=ax.projection)

    ax.text(sbx, sby-4500, '20 km', transform=ax.projection, fontsize=8)
    ax.text(sbx-12500, sby-4500, '10 km', transform=ax.projection, fontsize=8)
    ax.text(sbx-24500, sby-4500, '0 km', transform=ax.projection, fontsize=8)


# load the datasets
counties = gpd.read_file(os.path.abspath('data_files/NI_Counties.shp'))
lakes = gpd.read_file(os.path.abspath('data_files/Lakes.shp'))
places = gpd.read_file(os.path.abspath('data_files/NI_Places.shp'))
sites = gpd.read_file(os.path.abspath('data_files/Site_Locations.shp'))

# create a figure of size 10x10 (representing the page size in inches)
myFig = plt.figure(figsize=(10, 10))

myCRS = ccrs.UTM(29)  # create a Universal Transverse Mercator reference system to transform our data.
# NI is in UTM Zone 29, so we pass 29 to ccrs.UTM()

# transform data files to myCRS (UTM(29) has an epsg of 32629)
# transform all data_files to ensure all data is on the same reference system
counties.to_crs(epsg=32629, inplace=True)
lakes.to_crs(epsg=32629, inplace=True)
places.to_crs(epsg=32629, inplace=True)
sites.to_crs(epsg=32629, inplace=True)

# create NI_outline
NI_Union = unary_union(counties.geometry) # create NI outline by joining geometries of counties
NI_Outline = gpd.GeoDataFrame(geometry=gpd.GeoSeries(NI_Union)) # create GeoDataFrame based on GeoSeries output of union
NI_Outline.set_crs(epsg=32629, inplace=True) # set CRS of new outline GeoDataFrame

ax = plt.axes(projection=myCRS)  # finally, create an axes object in the figure, using a UTM projection,
# where we can actually plot our data.

# create outline
# first, we just add the outline of Northern Ireland using cartopy's ShapelyFeature
outline_feature = ShapelyFeature(NI_Outline['geometry'], myCRS, edgecolor='k', facecolor='w')

xmin, ymin, xmax, ymax = NI_Outline.total_bounds # find total bounds of NI outline (we will use this to zoom the map)
ax.add_feature(outline_feature)  # add the features we've created to the map.

# using the boundary of the shapefile features, zoom the map to our area of interest
ax.set_extent([xmin-10000, xmax+5000, ymin-5000, ymax+5000], crs=myCRS)  # because total_bounds
# gives output as xmin, ymin, xmax, ymax,
# but set_extent takes xmin, xmax, ymin, ymax, we re-order the coordinates here
# more of a gap after xmin to create space for legend

# pick colors, add features to the map
county_colors = ['violet', 'seagreen', 'cornflowerblue', 'firebrick', 'darkkhaki', 'orange']

# get a list of unique names for the county boundaries
county_names = list(counties.CountyName.unique())
county_names.sort()  # sort the counties alphabetically by name

# next, add the municipal outlines to the map using the colors that we've picked.
# here, we're iterating over the unique values in the 'CountyName' field.
# we're also setting the edge color to be black, with a line width of 0.5 pt.
# Feel free to experiment with different colors and line widths.
for ii, name in enumerate(county_names):
    feat = ShapelyFeature(counties.loc[counties['CountyName'] == name, 'geometry'],  # first argument is the geometry
                          myCRS,  # second argument is the CRS
                          edgecolor='k',  # outline the feature in black
                          facecolor=county_colors[ii],  # set the face color to the corresponding color from the list
                          linewidth=1,  # set the outline width to be 1 pt
                          alpha=0.25)  # set the alpha (transparency) to be 0.25 (out of 1)
    ax.add_feature(feat)  # once we have created the feature, we have to add it to the map using ax.add_feature()

# clip lakes to extent of Northern Ireland
NI_lakes = gpd.clip(lakes, NI_Outline, keep_geom_type=True)

# here, we're setting the edge color to be the same as the face color. Feel free to change this around,
# and experiment with different line widths.
lakes_feat = ShapelyFeature(NI_lakes['geometry'],  # first argument is the geometry
                            myCRS,  # second argument is the CRS
                            edgecolor='royalblue',  # set the edgecolor to be royalblue
                            facecolor='royalblue',  # set the facecolor to be royalblue
                            linewidth=1)  # set the outline width to be 1 pt
ax.add_feature(lakes_feat)  # add the collection of features to the map

# spatial join county to places
place_and_county = gpd.sjoin(places, counties, how='inner', lsuffix='left', rsuffix='right') # perform the spatial join

# create new column with nice names for county and nice names for places
for ind, row in place_and_county.iterrows():  # iterate over each row in the GeoDataFrame
    place_and_county.loc[ind, 'County'] = row['CountyName'].title()  # assign the row's CountyName to a new column
                               # that we will use to merge with additional place information - the values are changed
                               # to facilitate the merge as the place_and_county county names are all capitals

for ind, row in place_and_county.iterrows():  # iterate over each row in the GeoDataFrame
    place_and_county.loc[ind, 'Name'] = row['PLACENAME'].title()  # assign the row's PLACENAME to a new column
                               # that we will use to merge with additional place information - the values are changed
                               # to facilitate the merge as the place_and_county placenames are all capitals

# load csv
place_info = pd.read_csv('data_files/Place_information.csv')

# join with csv by merging on the shared variables (place and county)
places_wi = place_and_county.merge(place_info, on=["Name", "County"])
# we now have an updated places GeoDataFrame with info - places_wi

# pick colors, add features to the map
place_colors = ['indigo', 'darkmagenta', 'mediumorchid', 'mediumvioletred', 'orchid', 'palevioletred',
                'pink', 'lavenderblush']

# create a list of place types - we want these to be displayed in a specific order (size) so we will specify the order
# of the list here instead of getting unique values and sorting them
place_types = ['City', 'Town', 'Suburb', 'Large Village', 'Village', 'Hamlet', 'Townland', 'Location']
# we need to decide what shape symbol we want to use to display each type
place_symbol = ['H', 'H', 'o', 'o', 'h', 'h', 's', 's'] # create a list to use when plotting based on place_types order
# we need to decide the size of the symbol we want to use to display each type
place_symbol_size = [6, 6, 6, 6, 10, 10, 6, 6] # create a list to use when plotting based on place_types order

# next, add the municipal outlines to the map using the colors that we've picked.
# here, we're iterating over the unique values in the 'CountyName' field.
# we're also setting the edge color to be black, with a line width of 0.5 pt.
# Feel free to experiment with different colors and line widths.
# ShapelyFeature creates a polygon, so for point data we can just use ax.plot()
for ii, type in enumerate(place_types):
    plot = places_wi.loc[places_wi['Type'] == type]
    ax.plot(plot.geometry.x, plot.geometry.y, place_symbol[ii], color=place_colors[ii], ms=place_symbol_size[ii],
            transform=myCRS)

# sites
# pick colors, add features to the map
site_colors = ['orangered', 'olivedrab', 'dodgerblue']

# get a list of unique names for the county boundaries
site_names = list(sites.Name.unique())
site_names.sort()  # sort the counties alphabetically by name

# next, add the municipal outlines to the map using the colors that we've picked.
# here, we're iterating over the unique values in the 'CountyName' field.
# we're also setting the edge color to be black, with a line width of 0.5 pt.
# Feel free to experiment with different colors and line widths.
for ii, name in enumerate(site_names):
    feat = ShapelyFeature(sites.loc[sites['Name'] == name, 'geometry'],  # first argument is the geometry
                          myCRS,  # second argument is the CRS
                          edgecolor='k',  # outline the feature in black
                          facecolor=site_colors[ii],  # set the face color to the corresponding color from the list
                          linewidth=1) # set the outline width to be 1 pt
    ax.add_feature(feat)  # once we have created the feature, we have to add it to the map using ax.add_feature()

# generate handles for legend
# generate a list of handles for the county datasets
county_handles = generate_handles(counties.CountyName.unique(), county_colors, alpha=0.25)

# note: if you change the color you use to display lakes, you'll want to change it here, too
water_handle = generate_handles(['Lakes'], ['royalblue'])

# generate handles for the places datasets
place_type_handles = generate_handles_points(place_types, place_symbol, place_colors, place_symbol_size)
# we now have a list of lists for the point types on the map but we need to flatten it to get a list that can be used
# to add the handles to the legend as legend handles cannot be lists
place_type_handles_for_legend = [item for sublist in place_type_handles for item in sublist]

# update county_names to take it out of uppercase text
nice_names = [name.title() for name in county_names]

# create new list of place types for legend
place_types_u = place_types
place_types_u[-1] = 'Other' # change the type 'Location' which is the last type in the list (hence index -1) to 'Other'

# get centroids of polygons for site labels
sites["center"] = sites["geometry"].centroid
site_points = sites.copy()
site_points.set_geometry("center", inplace = True)

texts = []
for x, y, label in zip(site_points.geometry.x, site_points.geometry.y, site_points["Name"]):
    texts.append(plt.text(x, y, label, fontsize = 10))

# ax.legend() takes a list of handles and a list of labels corresponding to the objects you want to add to the legend
handles = county_handles + place_type_handles_for_legend + water_handle
labels = nice_names + place_types_u + ['Lakes']

leg = ax.legend(handles, labels, title='Map Legend', title_fontsize=12,
                fontsize=10, loc='upper left', frameon=True, framealpha=1)

gridlines = ax.gridlines(draw_labels=True,  # draw  labels for the grid lines
                         xlocs=[-8, -7.5, -7, -6.5, -6, -5.5],  # add longitude lines at 0.5 deg intervals
                         ylocs=[54, 54.5, 55, 55.5])  # add latitude lines at 0.5 deg intervals
gridlines.right_labels = False  # turn off the right-side labels
gridlines.top_labels = False  # turn off the top labels

# add the scale bar to the axis
scale_bar(ax)

# save the figure as map.png, cropped to the axis (bbox_inches='tight'), and a dpi of 300
myFig.savefig('Locator_map.png', bbox_inches='tight', dpi=300)