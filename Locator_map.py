import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from cartopy.feature import ShapelyFeature
import cartopy.crs as ccrs
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from shapely.ops import unary_union

# generate matplotlib handles to create a legend of the features we put in our map.
def generate_handles(labels, colors, edge='k', alpha=1):
    lc = len(colors)  # get the length of the color list
    handles = []
    for i in range(len(labels)):
        handles.append(mpatches.Rectangle((0, 0), 1, 1, facecolor=colors[i % lc], edgecolor=edge, alpha=alpha))
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

# create a figure of size 10x10 (representing the page size in inches)
myFig = plt.figure(figsize=(10, 10))

myCRS = ccrs.UTM(29)  # create a Universal Transverse Mercator reference system to transform our data.
# NI is in UTM Zone 29, so we pass 29 to ccrs.UTM()

# transform data files to myCRS (UTM(29) has an epsg of 32629)
# transform all data_files to ensure all data is on the same reference system
counties.to_crs(epsg=32629, inplace=True)
lakes.to_crs(epsg=32629, inplace=True)
places.to_crs(epsg=32629, inplace=True)

# create NI_outline
NI_Outline = unary_union(counties.geometry) # create NI outline by joining geometries of counties

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
county_colors = ['firebrick', 'seagreen', 'royalblue', 'coral', 'violet', 'cornsilk']

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
NI_Lakes = gpd.clip(lakes, NI_Outline, keep_geom_type=True)

# here, we're setting the edge color to be the same as the face color. Feel free to change this around,
# and experiment with different line widths.
lakes_feat = ShapelyFeature(NI_Lakes['geometry'],  # first argument is the geometry
                            myCRS,  # second argument is the CRS
                            edgecolor='mediumblue',  # set the edgecolor to be mediumblue
                            facecolor='mediumblue',  # set the facecolor to be mediumblue
                            linewidth=1)  # set the outline width to be 1 pt
ax.add_feature(lakes_feat)  # add the collection of features to the map


# spatial join county to places
place_and_county = gpd.sjoin(counties, places, how='inner', lsuffix='left', rsuffix='right') # perform the spatial join

# create new column with nice names for county and nice names for places




# load csv
place_info = pd.read_csv('data_files/Airports.csv')

# join with csv by merging on the shared variables (place and county)
places_wi = place_and_county.merge(place_info, on=["STATION", "DATE"])

# we now have an updated places GeoDataFrame with info places_wi




# ShapelyFeature creates a polygon, so for point data we can just use ax.plot()
# Use intermediate variable and .loc to select each location type to then add to ax.plot()
just_cities = places_wi.loc[places_wi['STATUS'] == 'City']
city_handle = ax.plot(just_cities.geometry.x, just_cities.geometry.y, 'D', color='r', ms=6, transform=myCRS)

just_towns = places_wi.loc[places_wi['STATUS'] == 'Town']
town_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)

just_suburbs = places_wi.loc[places_wi['STATUS'] == 'Suburb']
suburb_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)

just_large_village= places_wi.loc[places_wi['STATUS'] == 'Large Village']
large_village_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)

just_villages = places_wi.loc[places_wi['STATUS'] == 'Village']
village_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)

just_hamlets = places_wi.loc[places_wi['STATUS'] == 'Hamlet']
hamlet_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)

just_townlands = places_wi.loc[places_wi['STATUS'] == 'Townland']
townland_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)

just_others = places_wi.loc[places_wi['STATUS'] == 'Location']
other_handle = ax.plot(just_towns.geometry.x, just_towns.geometry.y, 's', color='g', ms=6, transform=myCRS)


# pick colors, add features to the map
place_colors = ['firebrick', 'seagreen', 'royalblue', 'coral', 'violet', 'cornsilk', 'cornsilk', 'cornsilk']

# get a list of unique names for the county boundaries
place_types = ['City', 'Town', 'Suburb', 'Large Village', 'Village', 'Hamlet', 'Townland', 'Location']

# next, add the municipal outlines to the map using the colors that we've picked.
# here, we're iterating over the unique values in the 'CountyName' field.
# we're also setting the edge color to be black, with a line width of 0.5 pt.
# Feel free to experiment with different colors and line widths.
for ii, type in enumerate(place_types):
    plot = places_wi.loc[places_wi['Type'] == type]
    ax.plot(plot.geometry.x, plot.geometry.y, 's', color=place_colors[ii], ms=6, transform=myCRS)





# generate handles for legend
# generate a list of handles for the county datasets
county_handles = generate_handles(counties.CountyName.unique(), county_colors, alpha=0.25)

# note: if you change the color you use to display lakes, you'll want to change it here, too
water_handle = generate_handles(['Lakes'], ['mediumblue'])

places_handles = generate_handles(places_wi.Type.unique(), place_colors, alpha=0.25)

# update county_names to take it out of uppercase text
nice_names = [name.title() for name in county_names]

# create new list just for legend
place_types_u = place_types
place_types_u[-1] = 'Other' # change the type 'Location' which is the last type in the list (hence index -1) to 'Other'




# ax.legend() takes a list of handles and a list of labels corresponding to the objects you want to add to the legend
handles = county_handles + water_handle + places_handles
labels = nice_names + ['Lakes'] + place_types_u

leg = ax.legend(handles, labels, title='Map Legend', title_fontsize=12,
                fontsize=10, loc='upper left', frameon=True, framealpha=1)

gridlines = ax.gridlines(draw_labels=True,  # draw  labels for the grid lines
                         xlocs=[-8, -7.5, -7, -6.5, -6, -5.5],  # add longitude lines at 0.5 deg intervals
                         ylocs=[54, 54.5, 55, 55.5])  # add latitude lines at 0.5 deg intervals
gridlines.right_labels = False  # turn off the right-side labels
gridlines.top_labels = False  # turn off the top labels

# add the text labels for the specific cities
for ind, row in just_cities.iterrows():  # just_cities.iterrows() returns the index and row
    x, y = row.geometry.x, row.geometry.y  # get the x,y location for each town
    ax.text(x, y, row['Name'].title(), fontsize=8, transform=myCRS)  # use plt.text to place a label at x,y

# add the scale bar to the axis
scale_bar(ax)

# save the figure as map.png, cropped to the axis (bbox_inches='tight'), and a dpi of 300
myFig.savefig('output_files/Locator_map.png', bbox_inches='tight', dpi=300)