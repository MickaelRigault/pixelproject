#! /usr/bin/env python
#

import numpy as np
UNIT_SQUARE = np.asarray([[0,0],[0,1],[1,1],[1,0]])-0.5

from propobject import BaseObject
from shapely import geometry
import pandas
import geopandas

class GridProjector( BaseObject ):

    PROPERTIES = ["gridin", "gridout"]
    DERIVED_PROPERTIES = ["gridinterest"]
    
    def __init__(self, grid_in=None, grid_out=None):
        """ """
        if grid_in is not None:
            self.set_grid(grid_in, "in")
        if grid_out is not None:
            self.set_grid(grid_out, "out")

    # =================== #
    #   Methods           #
    # =================== #
    # --------- #
    #  SETTER   #
    # --------- #
    def set_grid(self, grid, which="in"):
        """ """
        if which not in ["in","out"]:
            raise ValueError("Which should either be 'in' our 'out'")
        self._properties["grid%s"%which] = grid
        self._derived_properties["gridinterest"] = None

    def _measure_gridinterest_(self):
        """ """
        # -- internal -- #        
        def localdef_get_area(g):
            return g.area
        # -------------- #
        
        if self.gridin is not None and self.gridout is not None:
            self._derived_properties["gridinterest"] = geopandas.overlay(self.gridin.geodataframe, self.gridout.geodataframe, how='intersection')
            self.gridinterest["area"] = self.gridinterest["geometry"].apply(localdef_get_area)
        else:
            warnings.warn("Cannot measure gridinterest, because gridin and/or gridout is/are None")

    # -------------- #
    #  Measurement   #
    # -------------- #
    def project_data(self, data, as_serie=True):
        """ Use gridinteresect 

        Parameters
        ----------
        data: [ndarray or string or pandas.Serie]
            data associated to gridin that should be projected in gridout.
            could be:
            - ndarray: must have the same length as gridin
            - string: name of a gridin column (pandas)
            - pandas.Serie: serie that will be matched with gridin
        """
        # Calcul itself
        projected_data = self._project_data_( self._parse_data_(data) )

        if as_serie:            
            return projected_data

        projected_data_array = np.zeros( len(self.gridout.geodataframe) )
        projected_data_array[projected_data.index.values] = projected_data.values
        return projected_data_array
        
    def _project_data_(self, data):
        """ """
        self.gridinterest["_tmp"] = data[ self.gridin.geodataframe.loc[ self.gridinterest["id_1"]].index ] * self.gridinterest["area"]
        return self.gridinterest.groupby("id_2")["_tmp"].sum()

    def _parse_data_(self,data):
        """ 
        Parameters
        ----------
        data: [ndarray or string or pandas.Serie]
            data associated to gridin that should be projected in gridout.
            could be:
            - ndarray: must have the same length as gridin
            - string: name of a gridin column (pandas)
            - pandas.Serie: serie that will be matched with gridin
            
        Returns
        -------
        ndarray
        """
        if type(data) == str:
            if data not in self.gridin.geodataframe.columns:
                raise ValueError("Unknown gridin column '%s'"%data)
            return self.gridin.geodataframe[data].values
            
        elif type(data) == pandas.Series:
            return data.values
            
        elif len(data) != len(self.gridin.geodataframe):
            raise ValueError("data given as ndarray but lengthes do not match")
        
        return data
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def gridin(self):
        """ """
        return self._properties["gridin"]
    
    @property
    def gridout(self):
        """ """
        return self._properties["gridout"]

    @property
    def gridinterest(self):
        """ """
        if self._derived_properties["gridinterest"] is None:
            self._measure_gridinterest_()
        return self._derived_properties["gridinterest"]


    
class Grid( BaseObject ):
    """ """
    PROPERTIES = ["pixels", "shape"]
    SIDE_PROPERTIES = ["indexes"]
    DERIVED_PROPERTIES = ["vertices","geodataframe"]
    
    def __init__(self, pixels=None, shape=UNIT_SQUARE, indexes=None):
        """ """
        if pixels is not None:
            self.set_pixels(pixels,shape=shape)

        if indexes is not None:
            self.set_indexes(indexes)
            
    # =================== #
    #   Methods           #
    # =================== #
    @staticmethod
    def set_from(datainput):
        """ Create a new Grid objects from the given input data:
        
        Parameters
        ----------
        datainput: [geopandas.geodataframe.GeoDataFrame or ndarray]
            this could either be a:
            - geodataframe (and this calls self.set_geodataframe)
            - ndarray: if 3-shaped, this calls set_vertices ;
                       if 2-shaped, this calls set_pixels.
        
        Returns
        -------
        Grid
        """
        this = Grid()
        if type(datainput) == geopandas.geodataframe.GeoDataFrame:
            this.set_geodataframe(datainput)
            return this
        if type(datainput) == np.ndarray:
            if len(np.shape( datainput) ) == 3: # vertices
                this.set_vertices(datainput)
            elif len(np.shape( datainput) ) == 3: # pixels
                this.set_pixels(datainput)
            else:
                raise TypeError("cannot parse the shape of the given datainput")
            return this
        raise TypeError("cannot parse the format of the given input")
    
            
    # --------- #
    #  SETTER   #
    # --------- #
    def set_indexes(self, indexes, update=True):
        """ """
        if self.pixels is not None and len(indexes)[0] != self.npixels:
            raise AssertionError("not the same number of indexes as the number of pixels")
        self._side_properties["indexes"]  = indexes
        if update:
            self._update_geodataframe_()
        
    def set_pixels(self, pixels, shape=None, update=True):
        """ """
        # Setting the pixels
        if np.shape(pixels)[-1] != 2:
            raise ValueError("pixels must be [N,2] arrays")
        self._properties["pixels"] = np.asarray(pixels)
        
        if shape is not None:
            self.set_pixelshapes(shape, update=False)
            
        if update:
            self._update_geodataframe_()
            
    def set_pixelshapes(self, shape, update=True):
        """ """
        # Setting the pixel shape.s
        if len(np.shape(shape))==2:
            self._properties["shape"] = np.asarray(shape)
        elif len(np.shape(shape))==3:
            if self.pixels is not None and np.shape(shape)[0] != self.npixels:
                raise AssertionError("`shape` must be unique or have the same lenth as pixels")
            self._properties["shape"] = np.asarray(shape)
        else:
            raise ValueError("Cannot parse the given shape, must be [M,2] or [N,M,2] when N is the number of pixel and M the number of vertices")
        if update:
            self._update_geodataframe_()

    def set_vertices(self, vertices, overwrite=False, **kwargs):
        """ """
        if not overwrite and (self.pixels is not None and self.shape is not None):
            raise ValueError("Pixels and shape already defined. set the overwrite option to true, to update vertices")

        self._derived_properties["vertices"] = np.asarray(vertices)
        # Redefine pixels and shape
        pixels = np.mean(self.vertices,axis=1)
        shape = self.vertices - pixels[:,None]
        shape_unique = np.unique(shape, axis=0)
        if len(shape_unique)==1:
            shape = shape_unique[0]

        self.set_pixels(pixels, shape, **kwargs)

    def set_geodataframe(self, geodataframe, overwrite=False):
        """ """
        if not overwrite and (self.pixels is not None and self.shape is not None):
            raise ValueError("Pixels and shape already defined. set the overwrite option to true, to update geodataframe")

        if "geometry" not in geodataframe.columns:
            raise TypeError("The given geodataframe does not have 'geometry' column. It is required")
        
        self._derived_properties["geodataframe"] = geodataframe
        
        if "id" not in geodataframe.columns:
            self.geodataframe["id"] = self.indexes if self.pixels is not None else np.arange( len(geodataframe) )

        # - get the vertices:
        def get_verts(poly_):
            return np.stack(poly_.exterior.xy).T[:-1]
        
        vertices = np.stack( geodataframe["geometry"].apply(get_verts).values )
        self.set_vertices(vertices, update=False) # don't update the geodataframe
        
        
    # --------- #
    #  UPDATE   #
    # --------- #
    def _update_geodataframe_(self):
        """ """
        dataseries = self.get_geoseries()
        self._derived_properties["geodataframe"] = \
          geopandas.GeoDataFrame({'geometry': dataseries,
                                        'id':self.indexes})
        
    def add_data(self, data, name, indexes=None, inplace=True):
        """ """
        if indexes is None:
            indexes = self.indexes
        s_ = pandas.Series(data, name=name, index=indexes)
        if not inplace:
            return self.geodataframe.join(s_)
        self._derived_properties["geodataframe"] = self.geodataframe.join(s_)
        
    # --------- #
    #  GETTER   #
    # --------- #
    def get_geoseries(self):
        """ """
        import geopandas
        return geopandas.GeoSeries([geometry.Polygon(v) for v in self.vertices])
    
    # --------- #
    #  PLOTTER  #
    # --------- #
    def show(self, column=None, ax=None, edgecolor="0.7", facecolor="None", **kwargs):
        """ """
        if column is not None:
            facecolor=None
        return self.geodataframe.plot(column, ax=ax,facecolor=facecolor,edgecolor=edgecolor, **kwargs)


    
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def pixels(self):
        """ """
        return self._properties["pixels"]
    
    @property
    def npixels(self):
        """ """
        return len(self.pixels)
    
    @property
    def shape(self):
        """ """
        if self._properties["shape"] is None:
            self._properties["shape"] = UNIT_SQUARE
        return self._properties["shape"]

    # -- Side
    @property
    def indexes(self):
        """ """
        if self._side_properties["indexes"] is None:
            self._side_properties["indexes"] = np.arange(self.npixels)
        return self._side_properties["indexes"]
    
    # -- Derived
    
    @property
    def vertices(self):
        """ """
        if self._derived_properties["vertices"] is None and (self.pixels is not None and self.shape is not None):
            self._derived_properties["vertices"] = self.pixels[:,None]+self.shape
        return self._derived_properties["vertices"]
        
    @property
    def is_shape_unique(self):
        """ """
        return len(np.shape(self.shape))==2

    @property
    def geodataframe(self):
        """ """
        if self._derived_properties["geodataframe"] is None:
            self._update_geodataframe_()
        return self._derived_properties["geodataframe"]
