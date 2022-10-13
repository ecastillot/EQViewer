from turtle import width
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


def map_single_plot(region,
        catalogs=[],
        stations=[],
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
        legend=True,
        fig = None,
        save=True):

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
            rivers= ['2/blue','2/blue'],
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
            rivers= ['2/blue','2/blue'],
            water='lightblue',
            land=land_map,
            frame=["afg","WNse"],
        )

    if shapes_before_catalog:
        for shape in shapes_before_catalog:
            fig.plot(**shape.__dict__)
            

    # if isinstance(catalog,pd.DataFrame):
    if catalogs:
        new_catalogs = []
        for catalog in catalogs:
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
                
                fig.plot(
                        x=data["longitude"],
                        y=data["latitude"],
                        sizes=_size,
                        color=color,
                        cmap=True,
                        style=style,
                        pen=pen,
                    )

        
                with pygmt.config(FONT=6):
                    fig.colorbar(frame=["af",f'y+l"{cmap_label}"'],
                            position="JML+o1c/0c+w6c/0.25c+mc")

            else:
                fig.plot(
                    x=data["longitude"],
                    y=data["latitude"],
                    sizes=_size,
                    label=label,
                    color=color,
                    style=style,
                    pen=pen,
                )
    else:
        new_catalogs = []
    if stations:
        for station in stations:
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

    if fms:
        for fm in fms:

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
                _cmap_args.pop('color_target', None)
                pygmt.makecpt(**_cmap_args)

                fig.meca(
                        spec=data,
                        convention="aki",
                        scale = fm_scale,
                        C=True,
                        offset = fm_offset
                    )
                fig.colorbar(frame='af+l"Depth (m)"')

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
            fig.plot(**shape.__dict__)

    if isinstance(wells,dict):
        mech_state = wells["mech_sta_path"]
        color_mech_sta = wells["color_mech_sta"]
        cmap_wells = wells["cmap"]
        injection = wells["mech_sta_path"]
        xl = pd.read_excel(mech_state, sheet_name = None)
        for well, well_df in xl.items():

            if cmap_wells:
                _cmap_args = cmap_args.copy() 
                _cmap_args.pop('color_target', None)
                pygmt.makecpt(**_cmap_args)

                fig.plot(
                        x=well_df["Lon (°)"],
                        y=well_df["Lat (°)"],
                        cmap=True,
                        color=well_df["Z (m)"],
                        size=None,
                        style="s0.05",
                        pen=f"+0.0001p+i"
                        )
            else:

                if isinstance(color_mech_sta,dict):
                    possible_colors = list(color_mech_sta.keys())

                    if well not in possible_colors:
                        color="black"
                        proflabel=None
                    else:
                        color = color_mech_sta[well]
                        if legend:
                            proflabel=well
                        else:
                            proflabel=None
                        labeled = True
                else:
                    color=color_mech_sta
                    proflabel=None


                fig.plot(
                        x=well_df["Lon (°)"],
                        y=well_df["Lat (°)"],
                        label=proflabel,
                        pen=f"{color}"
                        )


    # if legend:
    #     # print(legend)
    #     # fig.legend(position='JML')
    #     fig.legend()

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
    
        
    # fig.show()
    return fig,new_catalogs

def profile_single_plot(region,catalogs,profile,
                fms=[],
                wells=None,
                depth=[0,3e3],
                figsize = ("12c", "12c"),
                depth_unit="m",
                color=None,
                reverse_xy=False,
                fig=None,
                save=True):

    if depth_unit == "m":
        divisor = 1
    elif depth_unit == "km":
        divisor = 1e3

    if reverse_xy:
        y_label = "Distance"
        x_label = "Depth"
        f_label= "ESnw"
        sharex = "t"
        sharey = False
    else:
        x_label = "Distance"
        y_label = "Depth"
        f_label= "WSen"
        sharex = False
        sharey="l"

    C,D = profile["coords"]
    width = profile["width"]
    name = profile["name"]

    r,a,ba = gps2dist_azimuth(C[1],C[0],D[1],D[0])

    if fig == None:
        fig = pygmt.Figure()   

    if reverse_xy:
        corte_region= depth + [0, r/divisor] 
        projection="X?/?"
    else:
        corte_region= [0, r/divisor] + depth
        projection="X?/-?"

    fig.basemap(region=corte_region, 
                projection=projection,
                frame=[f_label,f'xafg+l"{x_label} ({depth_unit})"', 
                        f'yafg{depth[1]}+l"{y_label} ({depth_unit})"'
                                ])

    fig.text(
                            position="BR",
                            text=f"{name[0]}-{name[1]}",
                        )

    personal_keys = ["color","label","size","style","pen"]
    for _catalog in catalogs:

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
            fig.colorbar(frame=["af",f'y+l"{cbar_label}"'],
                        position=colorbar_pos)

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
        color_mech_sta = wells["color_mech_sta"]
        mech_state = wells["mech_sta_path"]
        injection = wells["injection_path"]
        # well_label = wells["well_label"]

        injection = pd.read_csv(injection)
        injection["name"] = injection["name"].apply(lambda x: x.strip())

        xl = pd.read_excel(mech_state, sheet_name = None)
        for well, well_df in xl.items():
            data = well_df[["Lon (°)","Lat (°)","Z (m)", "MD (ft)","TVD (ft)"]]

            data = data[(data["Lon (°)"] >= region[0]) &\
                            (data["Lon (°)"] <= region[1]) ]
            data = data[(data["Lat (°)"] >= region[2]) &\
                            (data["Lat (°)"] <= region[3]) ]
            if data.empty:
                continue

            x_MD = data["MD (ft)"].to_numpy()
            x_TVD = data["TVD (ft)"].to_numpy()
            y_depth = data["Z (m)"].to_numpy()

            

            f_MD2depth= interpolate.interp1d(x_MD, y_depth, fill_value="extrapolate")
            f_TVD2depth = interpolate.interp1d(x_TVD, y_depth, fill_value="extrapolate")

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


            # print(df)
            df = df.rename(columns={0:"distance",1:"depth",
                                    2:"md",3:"tvd"})
            df["distance"] = df["distance"].astype(float)*1e3/divisor
            df["md"] = df["md"].astype(float)
            df["tvd"] = df["tvd"].astype(float)
            x_depth = df["depth"].to_numpy()
            y_distance = df["distance"].to_numpy()

            f_depth2distance= interpolate.interp1d(x_depth, y_distance, fill_value="extrapolate")

            ##injection
            _injection = injection[injection["name"] == well.strip()]
            if not _injection.empty:
                for i,row in _injection.iterrows():
                    depth_type = row["depth_type"]
                    min_depth = row["min"]
                    max_depth = row["max"]


                    if depth_type.upper() == "MD": 
                        inj_md = np.linspace(min_depth,max_depth,num=30,endpoint=True)
                        # inj_md = np.clip(inj_md, min(x_MD), max(x_MD))
                        inj_depth = f_MD2depth(inj_md)

                        # inj_depth = np.clip(inj_depth, min(x_depth), max(x_depth))
                        inj_distance = f_depth2distance(inj_depth)

                    elif depth_type.upper() == "TVD": 
                        inj_tvd = np.linspace(min_depth,max_depth,num=30,endpoint=True)
                        # inj_tvd = np.clip(inj_tvd, min(x_TVD), max(x_TVD))
                        inj_depth = f_TVD2depth(inj_tvd)

                        # inj_depth = np.clip(inj_depth, min(x_depth), max(x_depth))
                        inj_distance = f_depth2distance(inj_depth)


                    if (inj_depth.size != 0) and (inj_distance.size != 0):

                        if reverse_xy:
                            x=inj_depth
                            y=inj_distance
                        else:
                            x=inj_distance
                            y=inj_depth                        

                        fig.plot(
                        x=x,
                        y=y,
                        pen=f"8p,blue",
                        transparency=40
                        )     

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
                    proflabel=well
                    labeled = True
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

    if labeled:
        fig.legend( position='jTL')
    return fig

def xy_map_and_profile(region,
        catalogs=[],
        profile={},
        stations=[],
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
        depth=[0,3e3],
        depth_unit="m",
        figsize = ("12c", "12c"),
        legend_map=False,
        save=True):

    # if isinstance(wells,dict):
    #     color_mech_sta = wells["color_mech_sta"]
    #     mech_state = wells["mech_sta_path"]
    #     xl = pd.read_excel(mech_state, sheet_name = None)

    #     if isinstance(color_mech_sta,list):
    #         color_wells = {}
    #         for well, well_df in xl.items():
    #             r= np.random.randint(color_mech_sta[0])
    #             g = np.random.randint(color_mech_sta[1])
    #             b = np.random.randint(color_mech_sta[2])
    #             color = f"{r}/{g}/{b}"
    #             color_wells[well] = color
    #         wells["color_mech_sta"] = color_wells
        

    fig = pygmt.Figure()  
    with pygmt.clib.Session() as session:
        session.call_module('gmtset', 'FONT 8p')
        subplot =   fig.subplot(
                        nrows=2, ncols=2, 
                        # figsize=figsize, 
                        subsize=(8,8), 
                        sharex=False,  # shared x-axis on the bottom side
                        sharey=False,  # shared y-axis on the left side
                        margins=["0.1c","0.1c"]
                    ) 



    x_coords = ((region[0],(region[2]+region[3])/2),(region[1],(region[2]+region[3])/2))

    rx,a,ba = gps2dist_azimuth(x_coords[0][1],x_coords[0][0],
                              x_coords[1][1],x_coords[1][0])
    rx = rx/1e3

    
    y_coords = (((region[0]+region[1])/2,region[2]),
                ((region[0]+region[1])/2,region[3]))
    ry,a,ba = gps2dist_azimuth(y_coords[0][1],y_coords[0][0],
                              y_coords[1][1],y_coords[1][0])
    ry = ry/1e3

    if not profile:
        profile = {
                "colorline":"magenta",
                "color": "gray",
                "label":None,
                "size":None,
                "style" :None,
                "pen":"black",
                "cmap":False,
                "cbar_profile_args" : {"cmap":'viridis', 
                                        "color_target":"depth",
                                        "label":"Depth (m)",
                                        "reverse":False,
                                        "series":[0, 3e3] }
                }
    profile_x = profile.copy()
    profile_y = profile.copy()

    _profile_x = {"name":("x","x'"),      
                "coords":x_coords,
                "width":(-rx/2,rx/2),
                "colorline":"magenta",}
                
    _profile_y = {"name":("y","y'"),      
                "coords":y_coords,
                "width":(-ry/2,ry/2),
                "colorline":"magenta",
                }

    profile_x.update(_profile_x)
    profile_y.update(_profile_y)

    profiles = [profile_x,profile_y]

    with subplot:
        with fig.set_panel(panel=0):
            fig,catalogs = map_single_plot(region = region,
                    catalogs=catalogs,
                    stations=stations,
                    fms = fms,
                    wells=wells,
                    cmap_args = cmap_args,
                    relief=relief,
                    shapes_before_catalog=shapes_before_catalog,
                    shapes_after_catalog=shapes_after_catalog,
                    map_scale_args = map_scale_args,
                    legend=legend_map,
                    legend_args = legend_args,
                    land_map=land_map,
                    fig = fig,
                    save=False)
        with fig.set_panel(panel=1):
            fig = profile_single_plot(region = region,
                catalogs=catalogs,
                profile=profile_y,
                fms=fms,
                wells=wells,
                depth=depth,
                figsize = (("1c", "1c")),
                depth_unit=depth_unit,
                # color=color,
                reverse_xy=True,
                fig=fig,
                save=False)
        with fig.set_panel(panel=2):
            fig = profile_single_plot(region = region,
                catalogs=catalogs,
                profile=profile_x,
                fms=fms,
                wells=wells,
                depth=depth,
                figsize = (("1c", "1c")),
                depth_unit=depth_unit,
                # color=color,
                reverse_xy=False,
                fig=fig,
                save=False)
    
    # fig.show()
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