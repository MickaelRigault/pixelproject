#! /usr/bin/env python
#

""" Tools to take an astrobject insturment and create a slice of it """

import numpy as np

from shapely import geometry, vectorized
#
# Grid Project into another using shapely
#

def project_grid_into_grid(verts_1, verts_2, value_1):
    """ """
    # Base on geopandas
    import geopandas
    from shapely.geometry import Polygon
    # GeoSeries
    grid_gs_1 = geopandas.GeoSeries([Polygon(v) for v in verts_1])
    grid_gs_2 = geopandas.GeoSeries([Polygon(v) for v in verts_2])
    # GeoDataFrame
    df1 = geopandas.GeoDataFrame({'geometry': grid_gs_1, 'data':value_1, "id":np.arange(len(grid_gs_1))})
    df2 = geopandas.GeoDataFrame({'geometry': grid_gs_2, 'data':0,       "id":np.arange(len(grid_gs_2))})
    # Interact
    res_interact = geopandas.overlay(df1, df2, how='intersection')
    # Measure Overlap
    def get_area(g):
        return g.area
    res_union["area"] = res_union["geometry"].apply(get_area)
    res_union["wdata"] = res_union["area"]*res_union["data_1"]
    
    # The projected
    value_2_serie = res_union.groupby("id_2")["wdata"].sum()
    # Actual value2
    values_2 = np.zeros(len(verts_2))
    values_2[value_2_serie.index.values] = value_2_serie.values
    return values_2, res_union


# =============================== #
#  Slow Backup if not GeoPandas   #
# =============================== #
def project_to_grid(grid_1, val_1, grid_2 ):
    """ """
    weight = [overlap_to_verts(grid_1, v_ ) for v_ in grid_2]
    return np.sum(values*weights,axis=1)
    
def overlap_to_verts(grid_1, vert_, buffer = 2):
    """ """
    poly_ = geometry.Polygon(vert_)
    lgrid, sgrid, _should_be_2 = np.shape(grid_1)
        
    weight_map = np.zeros(lgrid)
    # Are the grid corners inside the buffered poly
    corners_in = vectorized.contains(poly_.buffer(buffer), *grid_1.reshape(lgrid*sgrid, 2).T
                                    ).reshape(lgrid,sgrid)
    
    corner_sum = np.sum(corners_in, axis=1)
    for i in np.argwhere(np.asarray(corner_sum, dtype="bool")).flatten():
        # 300microsec
        poly_vert_ = geometry.Polygon(verts_flat[i])
        if polygon.contains(poly_vert_):
            weight_map[i]= 1 
        else:
            weight_map[i]= poly_.intersection(poly_vert_).area
    
    return weight_map




