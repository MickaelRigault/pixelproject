#! /usr/bin/env python
#

""" Tools to take an astrobject insturment and create a slice of it """

import numpy as np

from shapely import geometry, vectorized
#
# Grid Project into another using shapely
#

def restride(arr, binfactor):
    """
    Rebin arr by binfactor.

    Let `arr.shape = (s1, s2, ...)` and `binfactor = (b1, b2, ...)` (same
    length), new shape will be `(s1/b1, s2/b2, ... b1, b2, ...)` (squeezed).

    * If binfactor is an iterable of length < arr.ndim, it is prepended with 1's.
    * If binfactor is an integer, it is considered as the bin factor for all axes.

    Bin 2D-array by a factor 2:
    >>> restride(np.ones((6, 8)), 2).shape
    (3, 4, 2, 2)

    Bin 2D-array by uneven factor (3, 2):
    >>> restride(np.ones((6, 8)), (3, 2)).shape
    (2, 4, 3, 2)

    Bin 3D-array by factor 2 over the last 2 axes, and take bin average:
    >>> q = np.arange(2*2*3*3*2).reshape(2, 2*3, 3*2)
    >>> restride(q, (2, 2)).mean(axis=(-1, -2))
    array([[[ 3.5,  5.5,  7.5],
            [15.5, 17.5, 19.5],
            [27.5, 29.5, 31.5]],

           [[39.5, 41.5, 43.5],
            [51.5, 53.5, 55.5],
            [63.5, 65.5, 67.5]]])

    Bin 3D-array by factor 2, and take bin average:
    >>> restride(q, 2).mean(axis=(-1, -2, -3))
    array([[21.5, 23.5, 25.5],
           [33.5, 35.5, 37.5],
           [45.5, 47.5, 49.5]])

    .. Note:: for a 2D-array, `restride(arr, (3, 2))` is equivalent to::

         np.moveaxis(arr.ravel().reshape(arr.shape[1]/3, arr.shape[0]/2, 3, 2), 1, 2)
    """
    # From Yannick Copin | y.copin@ipnl.in2p3.fr
    
    try:                        # binfactor is list-like
        # Convert binfactor to [1, ...] + binfactor
        binshape = [1] * (arr.ndim - len(binfactor)) + list(binfactor)
    except TypeError:           # binfactor is not list-like
        binshape = [binfactor] * arr.ndim
        
    assert len(binshape) == arr.ndim, "Invalid bin factor (shape)."
    assert (~np.mod(arr.shape, binshape).astype('bool')).all(), \
        "Invalid bin factor (modulo)."

    # New shape
    rshape = [ d // b for d, b in zip(arr.shape, binshape) ] + binshape
    # New stride
    rstride = [ d * b for d, b in zip(arr.strides, binshape) ] + list(arr.strides)

    rarr = np.lib.stride_tricks.as_strided(arr, rshape, rstride)

    return rarr.squeeze()       # Remove length-1 axes

    

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
    def localdef_get_area(l):
        return l.geometry.area/self.gridin.geodataframe.iloc[l.id_1].geometry.area

    res_interact["area"] = res_interact.apply(localdef_get_area, axis=1)
    res_interact["wdata"] = res_interact["area"]*res_interact["data_1"]
    
    # The projected
    value_2_serie = res_interact.groupby("id_2")["wdata"].sum()
    # Actual value2
    values_2 = np.zeros(len(verts_2))
    values_2[value_2_serie.index.values] = value_2_serie.values
    return values_2, res_interact


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




