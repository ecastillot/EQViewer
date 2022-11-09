import pygmt
import types
import random
from scipy import interpolate
import os
import string
import pandas as pd
import numpy as np
from . import utils as ut
import datetime as dt
import geopandas as gpd
from obspy.geodetics.base import gps2dist_azimuth
pygmt.config(FORMAT_GEO_MAP="ddd.xx")

def map(region,
        catalogs=[],
        stations=[],
        profiles=[],
        fms = [],
        wells=None,
        cmap_args = {"cmap":'rainbow', 
                        "reverse":True,
                        "series":[0, 2500] },
        relief=False,
        shapes_before_catalog=[],
        shapes_after_catalog=[],
        map_scale_args = {"xloc":None,
                        "yloc":None,
                        "distance":5},
        legend_args = {"xloc":None, #xloc,yloc are the x-y map coordinates to locate the scale bar.
                        "yloc":None},
        land_map="gray",
        rivers="['2/blue','2/blue']",
        borders='1/1p,black',
        legend=True,
        fig = None):

    if fig == None:
        fig = pygmt.Figure()   
        projection="M4i" 
    else:
        projection='M?',


    if relief:
        fig.grdimage(
            '@earth_relief_01s',
            region=region,
            projection=projection,
            cmap=True,
            shading=True,
            frame=["afg","WNse"]
        )
        fig.coast(
            region=region,
            # projection='M4i',
            projection=projection,
            shorelines=True,
            rivers= rivers,
            borders=borders,
            # water='lightblue',
            # land='grey',
            frame=["afg","WNse"],
        )
        # pygmt.makecpt(**cpt_kwargs)
    else:
        fig.coast(
            region=region,
            # projection='M4i',
            projection=projection,
            shorelines=True,
            borders=borders,
            rivers= rivers,
            water='lightblue',
            land=land_map,
            frame=["afg","WNse"],
        )

    if shapes_before_catalog:
        for shape in shapes_before_catalog:
            shape = shape.to_dict()
            print(shape)
            fig.plot(**shape)

    # if isinstance(catalog,pd.DataFrame):
    if catalogs:
        new_catalogs = []
        for catalog in catalogs:
            catalog = catalog.to_dict()
            data = catalog["data"]
            data = data[(data["longitude"] >= region[0]) &\
                            (data["longitude"] <= region[1]) ]
            data = data[(data["latitude"] >= region[2]) &\
                            (data["latitude"] <= region[3]) ]
            catalog["data"] = data
            new_catalogs.append(catalog.copy())


            color = catalog["color"]
            size = catalog["size"]
            style = catalog["style"]
            label = catalog["label"]
            cmap = catalog["cmap"]
            pen = catalog["pen"]

            if type(size) is types.LambdaType:
                _size = data.magnitude.apply(size)
            elif size == None:
                _size=size
            else:
                raise Exception("size parameter must be a lambda function")


            if cmap:
                # cmap_args = cmap_args.to_dict()
                color = data[cmap_args["color_target"]]
                cmap_label = cmap_args["label"]
                _cmap_args = cmap_args.copy() 

                if not cmap_args["series"]:
                    min_val = data[_cmap_args["color_target"]].min()
                    max_val = data[_cmap_args["color_target"]].max()
                    _cmap_args["series"] = [min_val,max_val]
                _cmap_args.pop('color_target', None)
                _cmap_args.pop('label', None)
                pygmt.makecpt(**_cmap_args)
                # print(data)
                fig.plot(
                        x=data["longitude"],
                        y=data["latitude"],
                        size=_size,
                        color=color,
                        cmap=True,
                        style=style,
                        pen=pen,
                    )
                    
                fig.colorbar(frame=f'af+l"{cmap_label}"',
                            position="JBC+e")

            else:
                fig.plot(
                    x=data["longitude"],
                    y=data["latitude"],
                    size=_size,
                    label=label,
                    color=color,
                    style=style,
                    pen=pen,
                )
    else:
        new_catalogs = []

    if stations:
        for station in stations:
            station = station.to_dict()
            fig.plot(
                x=station["data"]["longitude"],
                y=station["data"]["latitude"],
                color=station["color"],
                label=station["label"],
                style=station["style"],
                pen=station["pen"],
            )
            if station["name_in_map"]:
                fig.text(x=station["data"]["longitude"], 
                        y=station["data"]["latitude"], 
                        text=station["data"]["station"],
                        font=f"10p,Helvetica,black",
                        fill=None,offset="-0.05c/0.15c")
                        
    if wells:
        
        injection_cmap = wells.injection_cmap
        if injection_cmap != None:
            injection_cmap = injection_cmap.to_dict()
            if isinstance(injection_cmap,dict):
                inj_cmap_label = injection_cmap["label"]
                _injection_cmap = injection_cmap.copy() 

                if not injection_cmap["series"]:
                    injections = [ x.injection for x in wells.wells]
                    injections = pd.concat(injections)
                    if not injections.empty:
                        min_inj = injections[wells.injection_cmap.color_target].min()
                        max_inj = injections[wells.injection_cmap.color_target].max()
                        step = (max_inj-min_inj)/10
                        min_val = min_inj
                        max_val = max_inj
                        step = step
                        _injection_cmap["series"] = [min_val,max_val,step]
        else:
            _injection_cmap = None

            
        for _well in wells.wells:
            w = _well.to_dict()
            well = w["name"]
            well_df = w["data"]
            color_mech_sta = w["color"]
            cmap_wells = w["cmap"]
            injection = w["injection"]
            injection_cmap = w["injection_cmap"]

            if cmap_wells:
                cmap_args = cmap_args.to_dict()
                _cmap_args = cmap_args.copy() 
                cmap_label = _cmap_args["label"]
                _cmap_args.pop('color_target', None)
                _cmap_args.pop('label', None)
                pygmt.makecpt(**_cmap_args)

                fig.plot(
                        x=well_df["lon"],
                        y=well_df["lat"],
                        cmap=True,
                        color=well_df["z"],
                        size=None,
                        style="s0.05",
                        pen=f"+0.0001p+i"
                        )

                if not catalogs:
                    fig.colorbar(frame=f'af+l"{cmap_label}"',
                                position="JBC+e")
            else:

                if isinstance(color_mech_sta,dict):
                    possible_colors = list(color_mech_sta.keys())

                    if well not in possible_colors:
                        color="black"
                        label=None
                    else:
                        color = color_mech_sta[well]
                        label=well
                else: 
                    color = color_mech_sta
                    label=None


                fig.plot(
                        x=well_df["lon"],
                        y=well_df["lat"],
                        label=label,
                        pen=f"{color}"
                        )
            data = well_df[["lon","lat","z", "MD","TVD"]]
            
            data = data[(data["lon"] >= region[0]) &\
                                         (data["lon"] <= region[1]) ]
            data = data[(data["lat"] >= region[2]) &\
                            (data["lat"] <= region[3]) ]
            if data.empty:
                continue
            data = data.drop_duplicates(subset=["MD","z"])
            data = data.drop_duplicates(subset=["TVD","z"])
            x_MD = data["MD"].to_numpy()
            x_TVD = data["TVD"].to_numpy()
            lat = data["lat"].to_numpy()
            lon = data["lon"].to_numpy()

            f_MD2lat= interpolate.interp1d(x_MD, lat,
                                            kind="linear",
                                            fill_value="extrapolate")
            f_MD2lon= interpolate.interp1d(x_MD, lon,
                                            kind="linear",
                                            fill_value="extrapolate")
            f_TVD2lat= interpolate.interp1d(x_TVD, lat,
                                            kind="linear",
                                            fill_value="extrapolate")
            f_TVD2lon= interpolate.interp1d(x_TVD, lon,
                                            kind="linear",
                                            fill_value="extrapolate")
           
            
            if not injection.empty:
                _injection = injection
                for i,row in _injection.iterrows():
                    depth_type = row["depth_type"]
                    min_depth = row["min"]
                    max_depth = row["max"]
                    if wells.injection_cmap != None:
                        inj_size = [row[wells.injection_cmap.color_target]]*30
                    else:
                        inj_size = [row["water_flow"]]*30

                    if depth_type.upper() == "MD": 
                        inj_md = np.linspace(min_depth,max_depth,num=30,endpoint=True)
                        inj_md = np.clip(inj_md, min(x_MD), max(x_MD))
                        inj_lat = f_MD2lat(inj_md)
                        inj_lon = f_MD2lon(inj_md)
                    elif depth_type.upper() == "TVD": 
                        inj_tvd = np.linspace(min_depth,max_depth,num=30,endpoint=True)
                        inj_tvd = np.clip(inj_tvd, min(x_TVD), max(x_TVD))
                        inj_lat = f_TVD2lat(inj_tvd)
                        inj_lon = f_TVD2lon(inj_tvd)

                    if isinstance(_injection_cmap,dict):
                        _injection_cmap.pop('color_target', None)
                        _injection_cmap.pop('label', None)

                        pygmt.makecpt(**_injection_cmap)

                        fig.plot(
                            x=inj_lon,
                            y=inj_lat,
                            cmap=True,
                            color=inj_size,
                            size=None,
                            style="g0.3",
                            )
                        with pygmt.config(FORMAT_FLOAT_MAP="%.1e"):
                            fig.colorbar(frame=["af",f'y+l{inj_cmap_label}'],
                                        position='JMR+o1c/0c+e',
                                        transparency=80)
                    else:
                        fig.plot(
                            x=inj_lon,
                            y=inj_lat,
                            cmap=False,
                            color="blue",
                            size=None,
                            style="g0.3",
                            )
    
    if fms:
        for fm in fms:
            fm = fm.to_dict()

            fm_scale = str(fm["scale_for_m5"])
            data = fm["data"]
            color = fm["color"]
            cmap = fm["cmap"]
            main_n = fm["main_n"]

            data = data.rename(columns={f"dip_n{main_n}":"dip",
                                       f"strike_n{main_n}":"strike",
                                       f"rake_n{main_n}":"rake" })

            df_cols = data.columns.values.tolist()
            mf_cols = ["longitude","latitude","depth","strike","dip","rake",
                        "magnitude","plot_longitude","plot_latitude"]
            cols = list(set(df_cols) & set(mf_cols))
            data = data[cols]
            fm_offset = '1p'

            if cmap:
                color = data[cmap_args["color_target"]]
                _cmap_args = cmap_args.copy() 
                cmap_label = _cmap_args["label"]
                # _cmap_args["series"] 
                _cmap_args.pop('color_target', None)
                _cmap_args.pop('label', None)
                pygmt.makecpt(**_cmap_args)

                fig.meca(
                        spec=data,
                        convention="aki",
                        scale = fm_scale,
                        C=True,
                        offset = fm_offset
                    )
                if not catalogs:
                    fig.colorbar(frame=f'af+l"{cmap_label}"',
                                position="JBC+e")

            else:
                fig.meca(
                    spec=data,
                    convention="aki",
                    scale = fm_scale,
                    G=color
                    # plot_longitude=False,
                    # plot_latitude=False
                )

    if shapes_after_catalog:
        for shape in shapes_after_catalog:
            shape = shape.to_dict()
            fig.plot(**shape)

    if profiles:
        for k,profile in enumerate(profiles):
            profile = profile.to_dict()
            line = profile["coords"]
            width = profile["width"]
            colorline = profile["colorline"]
            C,D = line

            
            # fig.plot(x=[C[0], D[0]], y=[C[1], D[1]],
            #         projection="M", pen=f"2p,{colorline}")

            ul = ut.points_parallel_to_line(line,abs(width[0]))
            cul,dul = ul
            bl = ut.points_parallel_to_line(line,abs(width[1]),upper_line=False)
            cbl,dbl = bl

            # fig.plot(x=[cul[0], dul[0]], y=[cul[1], dul[1]],
            #         projection="M", pen=f"1.5p,{colorline},.")
            # fig.plot(x=[cbl[0], dbl[0]], y=[cbl[1], dbl[1]],
            #         projection="M", pen=f"1.5p,{colorline},.")

            fig.plot(x=[cbl[0], dbl[0],dul[0],cul[0],cbl[0]], 
                    y=[cbl[1], dbl[1],dul[1],cul[1],cbl[1]],
                    projection="M", pen=f"1.5p,{colorline},4_2:2p")


            if profile["name"] != None:
                ln,rn = profile["name"]
                fig.text(x=C[0], y=C[1], text=f"{ln}",
                        font=f"10p,Helvetica,black",
                        fill="white",offset="-0.05c/0.05c")
                        # fill=None,offset="-0.05c/0.05c")
                fig.text(x=D[0], y=D[1], text=f"{rn}",
                        font=f"10p,Helvetica,black",
                        fill="white",offset="0.05c/0.05c")
                        # fill=None,offset="0.05c/0.05c")

    if legend:
        xs = legend_args["xloc"]
        ys = legend_args["yloc"]
        if (xs != None) or (ys!= None):
            fig.legend(position=f"g{xs}/{ys}+o0.1c",box='+gwhite+p1p')
        else:
            fig.legend()
    
    xs = map_scale_args["xloc"]
    ys = map_scale_args["yloc"]
    d = map_scale_args["distance"]
    if (xs == None) or (ys== None):
        fig.basemap(map_scale=f"jBL+w{d}k+o0.5c/0.5c+f+lkm+at")
    else:
        fig.basemap(map_scale=f"g{xs}/{ys}+c{xs}/{ys}+w{d}k+o0.5c/0.5c+f+lkm+at")
    

    return fig,new_catalogs


def profile_plots(region,catalogs,profiles,
                fms=[],
                wells=None,
                depth=[0,3e3],
                subsize = ("12c", "12c"),
                figsize = None,
                cmap_args = {"cmap":'rainbow', 
                        "reverse":True,
                        "series":[0, 2500] },
                depth_unit="m",
                color=None,
                reverse_xy=False,
                margins=["1c","1c"],
                fig=None,
                legend=True):
    if fig == None:
        fig = pygmt.Figure()

    n = len(profiles)
    n_square = np.sqrt(n)
    c = int(n_square )
    r = int(n_square )
    if n ==2:
        c = 1
        r = 2
    elif n%n_square != 0:
        c += 1
        r += 1
        

    if depth_unit == "m":
        divisor = 1
    elif depth_unit == "km":
        divisor = 1e3

    cmaps = [True for x in profiles if x.to_dict()["cmap"]==True]

    if reverse_xy:
        y_label = "Distance"
        x_label = "Depth"
        f_label= "ESnw"
        sharex = "t"
        sharey = False
    else:
        x_label = "Distance"
        y_label = "Depth"
        f_label= "WSrt"
        # f_label= "WSen"
        sharex = "b"
        sharey="l"

    # if cmaps:
    #     if y_label=="Distance":
    #         margins = ["1.6c","0.1c"]
    #         # margins = ["0.1c","0.1c"]
    #     else:
    #         margins = ["0.1c","1.4c"]
    #         # margins = ["0.1c","0.1c"]
    # else:
    #     margins = ["0.1c","0.1c"]

    

    with pygmt.clib.Session() as session:
        session.call_module('gmtset', 'FONT 10p')
        subplot =   fig.subplot(
                        nrows=r, ncols=c, 
                        # figsize=figsize, 
                        subsize=subsize, 
                        # frame=[f'xafg200+l"{x_label} ({depth_unit})"', 
                        #         f'yafg200+l"{y_label} ({depth_unit})"',
                        #         # f'yafg{depth[1]}+l"{y_label} ({depth_unit})"',
                        #         f_label
                        #         ],
                        # frame=[
                        #         f'xafg{grid[0]}+l"{x_label} ({depth_unit})"', 
                        #         f'yafg{grid[1]}+l"{y_label} ({depth_unit})"',
                        #         # f'yafg{depth[1]}+l"{y_label} ({depth_unit})"',
                        #         f_label
                        #         ],
                        frame=f_label,
                        # sharex=sharex,  # shared x-axis on the bottom side
                        # sharey=sharey,  # shared y-axis on the left side
                        margins=margins
                    ) 
    personal_keys = ["color","label","size","style","pen"]

    with subplot:
        for f,profile in enumerate(profiles):
            profile = profile.to_dict()
            print(f+1,"->",profile["name"])
            C,D = profile["coords"]
            width = profile["width"]
            name = profile["name"]
            grid = profile["grid"]
            legend = profile["legend"]

            r,a,ba = gps2dist_azimuth(C[1],C[0],D[1],D[0])

            if reverse_xy:
                corte_region= depth + [0, r/divisor] 
                projection="x?/?"
                # projection="X?/?"
            else:
                corte_region= [0, r/divisor] + depth
                projection="x?/-?"
                # projection="X?/-?"

            if grid == None:
                grid1 = corte_region[:2]
                grid1 = round((grid1[1] -grid1[0])/5,2)
                grid2 = corte_region[2:]
                grid2 = round((grid2[1] -grid2[0])/5,2)
                grid = [grid1,grid2]
                # print(grid)

            # corte_region= [0, r/divisor] + depth
            fig.basemap(
                region=corte_region, 
                projection=projection, 
                frame=[
                    f'xafg{grid[0]}+l"{x_label} ({depth_unit})"', 
                                f'yafg{grid[1]}+l"{y_label} ({depth_unit})"',
                                # f'yafg{depth[1]}+l"{y_label} ({depth_unit})"',
                                f_label
                                ],
                panel=f
                        )
            with fig.set_panel(panel=f):

                fig.text(
                            position="BR",
                            text=f"{name[0]}-{name[1]}",
                        )

                for _catalog in catalogs:
                    _catalog = _catalog.to_dict()

                    catalog = _catalog.copy()

                    for key,value in profile.items():
                        if key in personal_keys:
                            if value != None:
                                catalog[key] = value

                    data = catalog["data"]
                    data = data[(data["longitude"] >= region[0]) &\
                            (data["longitude"] <= region[1]) ]
                    data = data[(data["latitude"] >= region[2]) &\
                                    (data["latitude"] <= region[3]) ]
                    catalog["data"] = data


                    color = catalog["color"]
                    cat_cmap = catalog["cmap"]
                    size = catalog["size"]
                    style = catalog["style"]
                    pen = catalog["pen"]

                    data = data.dropna(subset=["magnitude"])
                    data = data[["longitude","latitude","depth",
                                    "origin_time","magnitude"]]

                    try:
                        df = pygmt.project(
                                    data=data,
                                    unit=True,
                                    center=C,
                                    endpoint=D,
                                    convention="pz",
                                    width=width,
                                    verbose=True,
                                )
                    except:
                        continue

                    df = df.rename(columns={0:"distance",1:"depth",2:"origin_time",
                                                3:"magnitude"})
                    df["magnitude"] = df["magnitude"].astype(float)
                    df["distance"] = df["distance"].astype(float)*1e3/divisor
                    df["origin_time"] = pd.to_datetime(df["origin_time"])
                    df = df[["distance","depth","origin_time","magnitude"]]
                    df = df.sort_values("magnitude",ascending=False,ignore_index=True)

                    if type(size) is types.LambdaType:
                        _size = df.magnitude.apply(size)
                    elif size == None:
                        _size = size
                    else:
                        raise Exception("size parameter must be a lambda function")

                    if reverse_xy:
                        y=df.distance
                        x=df.depth
                        colorbar_pos = 'JMR+o2c/0c+e'

                    else:
                        colorbar_pos = 'JBC+e'
                        x=df.distance
                        y=df.depth

                    if profile["cmap"]:

                        cbar_profile_args = profile["cbar_profile_args"]
                        cbar_label = cbar_profile_args["label"]
                        color = df[cbar_profile_args["color_target"]]
                        _cbar_profile_args = cbar_profile_args.copy() 
                        
                        if not cbar_profile_args["series"]:
                            min_val = df[_cbar_profile_args["color_target"]].min()
                            max_val = df[_cbar_profile_args["color_target"]].max()
                            _cbar_profile_args["series"] = [min_val,max_val]
                       
                        _cbar_profile_args.pop('color_target', None)
                        _cbar_profile_args.pop('label', None)


                        pygmt.makecpt(**_cbar_profile_args)


                        fig.plot(
                                x=x,
                                y=y,
                                # projection="X?/-?", 
                                size=_size,
                                color=color,
                                cmap=True,
                                style=style,
                                pen=pen,
                            )
                        

                        # fig.colorbar(frame=f'af+l"{cbar_label}"')
                        # with pygmt.config(FORMAT_DATE_MAP="yyyy"):
                        fig.colorbar(frame=["af",f'y+l"{cbar_label}"',
                                        "sxa1y","pxa1o"],
                                    position=colorbar_pos)

                    else:
                        if cat_cmap:
                            color = df[cmap_args["color_target"]]
                            cmap_label = cmap_args["label"]
                            _cmap_args = cmap_args.copy() 

                            if not cmap_args["series"]:
                                min_val = data[_cmap_args["color_target"]].min()
                                max_val = data[_cmap_args["color_target"]].max()
                                _cmap_args["series"] = [min_val,max_val]
                            _cmap_args.pop('color_target', None)
                            _cmap_args.pop('label', None)

                            pygmt.makecpt(**_cmap_args)
                            fig.plot(
                                x=x,
                                y=y,
                                # projection="X?/-?", 
                                size=_size,
                                color=color,
                                cmap=True,
                                style=style,
                                pen=pen,
                                # pen="0.11",
                            )

                        else:
                            fig.plot(
                                x=x,
                                y=y,
                                # projection="X?/-?", 
                                size=_size,
                                color=color,
                                style=style,
                                pen=pen,
                                # pen="0.11",
                            )


                labeled = False
                if isinstance(wells,dict):
                    color_mech_sta = wells["color"]
                    mech_state = wells["mech_sta_path"]
                    injection = wells["injection_path"] 
                    injection_cmap = wells["injection_cmap"] 

                    if injection != None:
                        injection = pd.read_csv(injection)
                        injection["name"] = injection["name"].apply(lambda x: x.strip())
                        min_inj = injection["water_flow"].min()
                        max_inj = injection["water_flow"].max()
                        # step = (max_inj-min_inj)/10
                        step = 1e4
                    else:
                        injection = pd.DataFrame()


                    xl = pd.read_excel(mech_state, sheet_name = None)
                    for well, well_df in xl.items():
                        data = well_df[["lon","lat","z", "MD","TVD"]]

                        data = data[(data["lon"] >= region[0]) &\
                                         (data["lon"] <= region[1]) ]
                        data = data[(data["lat"] >= region[2]) &\
                                        (data["lat"] <= region[3]) ]
                        if data.empty:
                            continue

                        data = data.drop_duplicates(subset=["MD","z"])
                        data = data.drop_duplicates(subset=["TVD","z"])

                        x_MD = data["MD"].to_numpy()
                        x_TVD = data["TVD"].to_numpy()
                        y_depth = data["z"].to_numpy()

                        f_MD2depth= interpolate.interp1d(x_MD, y_depth,
                                                        kind="linear",
                                                        fill_value="extrapolate")
                        f_TVD2depth = interpolate.interp1d(x_TVD, y_depth, 
                                                        kind="linear",
                                                        fill_value="extrapolate")

                        try:
                            df = pygmt.project(
                                        data=data,
                                        unit=True,
                                        center=C,
                                        endpoint=D,
                                        convention="pz",
                                        width=width,
                                        verbose=True,
                                    )
                        except:
                            continue
                        


                        df = df.rename(columns={0:"distance",1:"depth",
                                                2:"md",3:"tvd"})
                        df["distance"] = df["distance"].astype(float)*1e3/divisor
                        df["md"] = df["md"].astype(float)
                        df["tvd"] = df["tvd"].astype(float)

                        # df = df[(df["distance"] >= width[0]*1e3) &\
                        #                  (df["distance"] <= width[1]*1e3) ]
                        
                        if df.empty:
                            continue
                        x_depth = df["depth"].to_numpy()
                        y_distance = df["distance"].to_numpy()
                        if len(x_depth) <2:
                            continue
                        f_depth2distance= interpolate.interp1d(x_depth, y_distance,
                                                            kind="linear",
                                                            fill_value="extrapolate")
                        # f_depth2distance= interpolate.splrep( x_depth,y_distance,k=5)


                        ##injection
                        if not injection.empty:

                            _injection = injection[injection["name"] == well.strip()]
                            if not _injection.empty:
                                for i,row in _injection.iterrows():
                                    depth_type = row["depth_type"]
                                    min_depth = row["min"]
                                    max_depth = row["max"]
                                    
                                    if depth_type.upper() == "MD": 
                                        inj_md = np.linspace(min_depth,max_depth,num=30,endpoint=True)
                                        # step = abs(max_depth-min_depth)/30
                                        # inj_md = np.arange(min_depth,max_depth+step,step)
                                        inj_md = np.clip(inj_md, min(x_MD), max(x_MD))
                                        inj_depth = f_MD2depth(inj_md)

                                        inj_depth = np.clip(inj_depth, min(x_depth), max(x_depth))
                                        inj_distance = f_depth2distance(inj_depth)
                                        # inj_distance = interpolate.splev(inj_depth, f_depth2distance)

                                    elif depth_type.upper() == "TVD": 
                                        inj_tvd = np.linspace(min_depth,max_depth,num=30,endpoint=True)
                                        # step = abs(max_depth-min_depth)/30
                                        # inj_tvd = np.arange(min_depth,max_depth+step,step)
                                        inj_tvd = np.clip(inj_tvd, min(x_TVD), max(x_TVD))
                                        inj_depth = f_TVD2depth(inj_tvd)

                                        inj_depth = np.clip(inj_depth, min(x_depth), max(x_depth))
                                        inj_distance = f_depth2distance(inj_depth)
                                        # inj_distance = interpolate.splev(inj_depth, f_depth2distance)


                                    if (inj_depth.size != 0) and (inj_distance.size != 0):

                                        if reverse_xy:
                                            x=inj_depth
                                            y=inj_distance
                                        else:
                                            x=inj_distance
                                            y=inj_depth                        

                                        x = list(dict.fromkeys(x))
                                        y = list(dict.fromkeys(y))
                                        inj_size = [row["water_flow"]]*len(x)

                                        if len(x) ==1:
                                            continue
                                        # fig.plot(
                                        # x=x,
                                        # y=y,
                                        # pen=f"8p,blue",
                                        # transparency=40
                                        # )  


                                        # color = df[injection_cmap["color_target"]]
                                        # inj_cmap_label = injection_cmap["label"]
                                        if isinstance(injection_cmap,dict):
                                            _injection_cmap = injection_cmap.copy() 

                                            if not injection_cmap["series"]:
                                                min_val = min_inj
                                                max_val = max_inj
                                                step = step
                                                _injection_cmap["series"] = [min_val,max_val,step]
                                            
                                            _injection_cmap.pop('label', None)

                                            pygmt.makecpt(**_injection_cmap)
                                            # pygmt.makecpt(cmap='cool',
                                            #             reverse=True,
                                            #             transparency=80,
                                            #             series=[min_inj,max_inj,step])
                                            fig.plot(
                                                x=x,
                                                y=y,
                                                cmap=True,
                                                color=inj_size,
                                                size=None,
                                                style="g0.3",
                                                )   
                                        else:
                                            pass

                            else:
                                pass


                        ##mech state
                        if reverse_xy:
                            y=df["distance"]
                            x=df["depth"]
                        else:
                            x=df["distance"]
                            y=df["depth"]

                        if isinstance(color_mech_sta,dict):
                            possible_colors = list(color_mech_sta.keys())

                            if well not in possible_colors:
                                color="black"
                                proflabel=None
                            else:
                                color = color_mech_sta[well]
                                if legend:
                                    proflabel=well
                                    labeled = True
                                else:
                                    proflabel=None
                                    labeled = False

                        else:
                            color=color_mech_sta
                            proflabel=None


                        fig.plot(
                                x=x,
                                y=y,
                                label=proflabel,
                                pen=f"{color}"
                                )    
                        

                if fms:
                    for fm in fms:
                        fm_scale = str(fm["scale_for_m5"])
                        data = fm["data"]
                        cmap = fm["cmap"]
                        main_n = fm["main_n"]
                        color = fm["color"]

                        data = data.rename(columns={f"dip_n{main_n}":"dip",
                                                f"strike_n{main_n}":"strike",
                                                f"rake_n{main_n}":"rake" })

                        df_cols = data.columns.values.tolist()
                        mf_cols = ["longitude","latitude","depth","strike","dip","rake",
                                    "magnitude","plot_longitude","plot_latitude","origin_time"]
                        cols = list(set(df_cols) & set(mf_cols))
                        data = data[cols]
                        fm_offset = '1p'

                        data = data[["longitude","latitude","depth",
                                    "origin_time","magnitude",
                                    "strike","dip","rake"]]
                        try:
                            df = pygmt.project(
                                        data=data,
                                        unit=True,
                                        center=C,
                                        endpoint=D,
                                        convention="pz",
                                        width=width,
                                        verbose=True,
                                    )
                        except:
                            continue

                        df = df.rename(columns={0:"longitude",1:"latitude",2:"origin_time",
                                                    3:"magnitude",4:"strike",5:"dip",6:"rake"})
                        
                        df["longitude"] = df["longitude"].astype(float)*1e3/divisor
                        df["depth"]= df["latitude"]
                        # print(df)
                        if reverse_xy:
                            df["_latitude"]=df["latitude"]
                            df["latitude"]=df["longitude"]
                            df["longitude"]= df["_latitude"]
                            df = df.drop('_latitude', 1)

                            colorbar_pos = 'JMR+o2c/0c+e'
                        else:
                            colorbar_pos = 'JBC+e'

                        if cmap:
                            cbar_profile_args = profile["cbar_profile_args"]
                            cbar_label = cbar_profile_args["label"]
                            color = df[cbar_profile_args["color_target"]]
                            _cbar_profile_args = cbar_profile_args.copy() 
                            
                            if not cbar_profile_args["series"]:
                                min_val = df[_cbar_profile_args["color_target"]].min()
                                max_val = df[_cbar_profile_args["color_target"]].max()
                                _cbar_profile_args["series"] = [min_val,max_val]
                        
                            _cbar_profile_args.pop('color_target', None)
                            _cbar_profile_args.pop('label', None)


                            pygmt.makecpt(**_cbar_profile_args)
                            
                            df = df[["longitude","latitude","depth","magnitude",
                                    "strike","dip","rake"]]

                            fig.meca(
                                    spec=df,
                                    convention="aki",
                                    scale = fm_scale,
                                    C=True,
                                    offset = fm_offset
                                )
                            # fig.colorbar(frame='af+l"Depth (m)"')
                            fig.colorbar(frame=["af",f'y+l"{cbar_label}"'],
                                    position=colorbar_pos)

                        else:
                            df = df[["longitude","latitude","depth","magnitude",
                                    "strike","dip","rake"]]

                            fig.meca(
                                spec=df,
                                convention="aki",
                                scale = fm_scale,
                                G=color
                                # plot_longitude=False,
                                # plot_latitude=False
                            )

                if legend:
                    if labeled == True:
                        fig.legend( position='jTL')
    fig.shift_origin(yshift="h-1c")
    return fig

if __name__ == "__main__":
    from shape import ShapeObject
    lats = 3.78;latn = 3.94;lonw = -73.72;lone = -73.57 #castilla
    reg = [lonw , lone, lats, latn ]

    field_path = "/home/emmanuel/Ecopetrol/Sismologia_Ecopetrol/data/shapes_before_catalog/castilla/castilla_chichimene.shp"
    castilla_field = ShapeObject(path=field_path,
                        attrs={"FIELD_NAME":["Castilla",
                                            "Castilla Norte",
                                            "Castilla Este"]},
                        label='"Castilla fields"',
                        pen="1p,red"                                        
                        ) 

    block_path = "/home/emmanuel/G-Ecopetrol/ecastillo/Avances/2022/Asociacion_castilla/qgis/Poligono campos/cubarral.shp"
    castilla_block = ShapeObject(path=block_path,
                        attrs={"CONTRATO_N":["CUBARRAL"]},
                        label='"Cubarral block"',
                        pen="1p,black"                                        
                        ) 

    shapes_before_catalog = [castilla_field,castilla_block]



    
    shapes_after_catalog = [castilla_field,castilla_block]
    print(shapes_before_catalog)
    map(reg,shapes_before_catalog=shapes_before_catalog)
    # map(reg,shapes=shapes,relief=True)