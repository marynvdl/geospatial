import geopandas as gpd
import matplotlib
import matplotlib.pyplot as plt
import shapely.geometry
import shapely.ops
import pandas as pd
import folium
import earthpy
from ipyleaflet import Map, GeoData, basemaps, LayersControl
import json
from ipywidgets import HTML
import numpy as np
import multiprocessing
import time
from functools import partial
import os
from os import listdir
from os.path import isfile, join

mypath = 'data/river/output'
output_files = [f[:-8] for f in listdir(mypath)]


#-----------------------------
# Buffering and clipping rivers
#-----------------------------

river_results = gpd.GeoDataFrame()

parks = gpd.read_file("data/river/161pas.shp").to_crs(epsg=4269)

for index, park in parks.iterrows():
    park_name = park.loc['name']
    park_file_name = ''.join(e for e in park_name if e.isalnum())

    if park.loc['name'] not in output_files or park_file_name not in output_files:
        print(park_name)

        # Buffer
        park_buffer = parks.loc[parks['name'] == park_name].to_crs(
            '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ')

        park_buffer['geometry'] = park_buffer['geometry'].geometry.buffer(500000)
        park_buffer = park_buffer.to_crs(epsg=4269)

        # Park
        park = parks.loc[parks['name'] == park_name].to_crs(
            '+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs ')
        park['geometry'] = park['geometry'].geometry.buffer(500)

        park = park.to_crs(epsg=4269)

        # Rivers
        try:
            rwa_rivers
        except NameError:
            rwa_rivers = gpd.read_file("data/river/161clipped_rivers.shp").to_crs(epsg=4269)
            print('importing')

        # Buffer rivers
        buffer_border = gpd.GeoDataFrame(geometry=park_buffer.boundary)
        buffer_rivers = gpd.sjoin(rwa_rivers, park_buffer, how='inner', predicate='within')

        # Park rivers
        park_rivers = gpd.sjoin(rwa_rivers, park, how='inner', predicate='within')

        park_rivers = park_rivers.reset_index(drop=True)

        # Set property if in park
        buffer_rivers['inpark_par'] = 0

        for index1, river1 in buffer_rivers.iterrows():
            if river1['NOID'] in park_rivers['NOID'].values:
                buffer_rivers.at[index1, 'inpark_par'] = 1

        buffer_rivers = buffer_rivers.rename(columns={'index_left': 'i_l'})
        buffer_rivers = buffer_rivers.rename(columns={'index_right': 'i_r'})
        buffer_rivers = buffer_rivers.reset_index(drop=True)

        # Rivers entering the park from outside
        buffer_rivers['o_outside_dis'] = 0.0
        buffer_rivers['o_inflow_par'] = 0
        buffer_rivers['o_through_par'] = 0
        buffer_rivers['o_park_dis'] = 0.0

        # Rivers starting in park
        buffer_rivers['i_start_dis'] = 0.0
        buffer_rivers['i_park_dis'] = 0.0
        buffer_rivers['i_park_fra'] = 0.0

        # All rivers
        buffer_rivers['park_dis'] = 0.0
        buffer_rivers['intersect_par'] = 0
        buffer_rivers['park_fra'] = 0.0

        # Rivers on boundary of park
        park_rivers['intersect_par'] = 0
        park_rivers['o_outside_dis'] = 0.0
        buffer_rivers['o_park_fra'] = 0.0

        park_border = gpd.GeoDataFrame(geometry=park.boundary)
        all_inters = gpd.sjoin(buffer_rivers, park_border, op='intersects')

        for index1, river1 in buffer_rivers.iterrows():
            if river1['NOID'] in all_inters['NOID'].values:
                buffer_rivers.at[index1, 'intersect_par'] = 1
                buffer_rivers.at[index1, 'inpark_par'] = 1

        for index1, river1 in buffer_rivers.iterrows():
            if river1['NOID'] in park_rivers['NOID'].values:
                buffer_rivers.at[index1, 'inpark_par'] = 1

        for index1, river1 in all_inters.iterrows():
            if river1['NOID'] in park_rivers['NOID'].values:
                park_rivers.at[index1, 'intersect_par'] = 1

        all_inters = all_inters.reset_index(drop=True)
        inters = gpd.GeoDataFrame()
        for index1, r1 in all_inters.iterrows():
            already_intersected = 0
            for index2, r2 in all_inters.iterrows():
                if r1['NOID'] == r2['NDOID']:
                    already_intersected = 1
            if already_intersected != 1:
                inters = inters.append(all_inters.iloc[index1])

        # Make df of all rivers inside and intersecting park
        in_park = park_rivers.append(inters)


        #-----------------------------
        # Find upstream river outside park
        #-----------------------------

        in_inters = gpd.GeoDataFrame()
        inters = inters.reset_index(drop=True)

        # Loop through all rivers intersecting with park border
        for index1, r1 in inters.iterrows():
            # All rivers starting
            if r1['NUOID'] is None:
                pass
            # All rivers not starting
            else:
                for index2, r2 in buffer_rivers.loc[buffer_rivers['NDOID']==r1['NOID']].iterrows():
                    # If river (r2) flows into river on border (r1)
                    if r1['NOID'] == r2['NDOID']:
                        # to make sure upstream river (r2) is outside the park
                        if r2['inpark_par'] != 1:
                            # then upstream river (r2) is flowing into the park
                            inters.at[index1, 'o_inflow_par'] = 1
                            in_inters = in_inters.append(inters.iloc[index1])


        if inters.shape[0] > 0:

            # Set the o_inflow_par of buffer_rivers also to 1
            buffer_rivers.loc[buffer_rivers['NOID'].isin(inters.loc[inters['o_inflow_par']==1,'NOID']), 'o_inflow_par'] = 1

            # If there are rivers flowing into the park
            if in_inters.shape[0] > 0:

                # drop all duplicate in-flowing rivers
                in_inters = in_inters.drop_duplicates()
                print('Number of inflowing rivers: ' + str(len(in_inters)))

                # set property on buffer_rivers to know that river is entering the park from outside
            # may be REDUNDANT
                for index1, river1 in buffer_rivers.iterrows():
                    if river1['NOID'] in in_inters['NOID'].values:
                        buffer_rivers.at[index1, 'o_inflow_par'] = 1


            #-----------------------------
            # Calculating outside flow
            #-----------------------------


            if in_inters.shape[0] > 0:

                for index1, r1 in in_inters.iterrows():
                    print("\n")
                    print(r1['NOID'])
                    print('start_flow: ' + str(r1['DIS_AV_CMS']))

                    run = 0
                    # If the river has reached the park boundary and is now flowing out
                    out_flow = 0
                    new_inters = gpd.GeoDataFrame()

                    # While the river is still inside the park
                    while out_flow == 0 and run < 100000:
                        run += 1
                        print('run' + str(run))
                        # If in-flowing rivers on park boundary
                        if run == 1:
                            # For the next downstream rivers
                            for index2, r2 in buffer_rivers.loc[buffer_rivers['NOID'] == r1['NDOID']].iterrows():
                                # If inflowing river (r1) is flowing into downstream river (r2)
                                if r1['NDOID'] == r2['NOID']:

                                    # Set the downstream river (r2)'s outside flow to the inflowing river (r1)'s flow
                                    buffer_rivers.at[index2, 'o_outside_dis'] = r2['o_outside_dis'] + r1['DIS_AV_CMS']

                                    print('o_outside_dis: ' + str(buffer_rivers.at[index2, 'o_outside_dis']))

                                    # Add r2 to become the next round's upstream river
                                    new_inters = new_inters.append(buffer_rivers.iloc[index2])

                                    # Check if r2 is already crossing the park boundary (leaving the park)
                                    if r2['inpark_par'] == 0:
                                        out_flow = 1
                                    else:
                                        buffer_rivers.at[index2, 'o_through_par'] = 1


                        # If not the in-flowing river, but a next-in-line downstream river

                        else:
                            for index3, r3 in new_inters.iterrows():
                                new_inters = gpd.GeoDataFrame()
                                for index4, r4 in buffer_rivers.loc[buffer_rivers['NOID'] == r3['NDOID']].iterrows():
                                    # If upstream river (r3) is flowing into downstream river (r4)
                                    if r3['NDOID'] == r4['NOID']:

                                        # Set the downstream river (r4)'s outside flow to the inflowing river (r3)'s flow
                                        buffer_rivers.at[index4, 'o_outside_dis'] = r4['o_outside_dis'] + r1['DIS_AV_CMS']

                                        print('o_outside_dis: ' + str(buffer_rivers.at[index4, 'o_outside_dis']))

                                        # Add r4 to become the next round's upstream river
                                        new_inters = new_inters.append(buffer_rivers.iloc[index4])

                                        # Check if r4 is already crossing the park boundary (leaving the park)
                                        if r4['inpark_par'] == 0:
                                            out_flow = 1

                                        else:
                                            buffer_rivers.at[index4, 'o_through_par'] = 1



            #-----------------------------
            # Calculate park contribution
            #-----------------------------

            for index1, r1 in buffer_rivers.iterrows():
                if (r1['o_outside_dis'] != 0) and r1['inpark_par'] == 1:
                    buffer_rivers.at[index1, 'o_park_dis'] = r1['DIS_AV_CMS'] - r1['o_outside_dis']
                    buffer_rivers.at[index1, 'o_park_fra'] = buffer_rivers.at[index1, 'o_park_dis'] / r1['DIS_AV_CMS']


            #-----------------------------
            # Combine with rivers starting inside park
            #-----------------------------

            # Loop through all rivers intersecting with park border

            new_inters = inters
            print('No of inters:' + str(inters.shape[0]))
            runs = 500
            run = 0
            while run < runs:
                run += 1
                print(round(100 * run / runs))

                new_new_inters = new_inters
                new_inters = gpd.GeoDataFrame()

                for index1, r1 in new_new_inters.iterrows():

                    # RIVERS STARTING IN PARK

                    # If river originates in park, and not flowing from the outside in.
                    if run == 1 and r1['o_inflow_par'] == 0:
                        # Set its own i_start_dis to its dis
                        buffer_rivers.loc[buffer_rivers['NOID'] == r1['NOID'], 'i_start_dis'] = r1['DIS_AV_CMS'] + r1['i_start_dis']


                    # If the river is downstream from another
                    else:
                        # Subset all rivers upstream from r1
                        upstream_ids1 = map(int, r1['NUOID'].split('_'))
                        upstream1 = buffer_rivers.loc[buffer_rivers['NOID'].isin(upstream_ids1)].reset_index()

                        # r1's i_start_dis will be equal to the sum of all the other's i_start_dis
                        buffer_rivers.at[index1, 'i_start_dis'] = upstream1['i_start_dis'].sum()

                        if (r1['o_through_par'] == 0) and (r2['intersect_par'] == 0):
                            buffer_rivers.at[index1, 'o_park_dis'] = upstream1['o_park_dis'].sum()
                            buffer_rivers.at[index1, 'o_park_fra'] = upstream1['o_park_fra'].sum()

                        # POPULATE ALL RIVERS DOWNSTREAM

                    for index2, r2 in buffer_rivers.loc[buffer_rivers['NOID'] == r1['NDOID']].iterrows():

                        if (r1['NDOID'] == r2['NOID']):
                            new_inters = new_inters.append(buffer_rivers.iloc[index2])



            buffer_rivers['i_park_fra'] = buffer_rivers['i_start_dis']/buffer_rivers['DIS_AV_CMS']
            buffer_rivers['park_dis'] = buffer_rivers['o_park_dis'] + buffer_rivers['i_start_dis']
            buffer_rivers['park_fra'] = buffer_rivers['park_dis']/buffer_rivers['DIS_AV_CMS']



            for index1, r1 in buffer_rivers.iterrows():
                if buffer_rivers.at[index1, 'inpark_par'] == 1 and buffer_rivers.at[index1, 'o_through_par'] != 1 and buffer_rivers.at[index1, 'o_inflow_par'] == 0:
                    buffer_rivers.at[index1, 'park_fra'] = 1


            relevant_rivers = buffer_rivers.loc[(buffer_rivers['park_fra'] > 0) | (buffer_rivers['inpark_par'] == 1)]
            relevant_rivers['park'] = park_name

            park_file_name = ''.join(e for e in park_name if e.isalnum())
            file_name = 'data/river/output/' + park_file_name + '.geojson'



            relevant_rivers.to_file(file_name, driver='GeoJSON')
            print('Done')
