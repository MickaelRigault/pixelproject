#! /usr/bin/env python
#

import warnings
import numpy as np
UNIT_SQUARE = np.asarray([[0,0],[0,1],[1,1],[1,0]])-0.5

from propobject import BaseObject
from shapely import geometry
import pandas
import geopandas

# ======================= #
#                         #
#    Functions            #
#                         #
# ======================= #

def get_simple_grid(xbounds, ybounds, shift_origin=None):
    """ """
    xbounds = np.atleast_1d(xbounds)
    if len(xbounds)==1:
        xmin,xmax = 0,xbounds[0]
    else:
        xmin,xmax = xbounds
        
    ybounds = np.atleast_1d(ybounds)
    if len(ybounds)==1:
        ymin,ymax = 0,ybounds[0]
    else:
        ymin,ymax = ybounds

    pixels  = np.mgrid[xmin:xmax,ymin:ymax]
    pixels2_flat = np.concatenate(pixels.T, axis=0)
    if shift_origin is not None:
        # not += because conflict between int and float array
        pixels2_flat = pixels2_flat+ shift_origin 
    return Grid(pixels2_flat, UNIT_SQUARE)

# ======================= #
#                         #
#    Classes              #
#                         #
# ======================= #

class GridProjector( BaseObject ):
    """ """
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
        def localdef_get_area(l):
            return l.geometry.area/self.gridin.geodataframe.iloc[l.id_1].geometry.area
        # -------------- #
        
        if self.gridin is not None and self.gridout is not None:
            #
            # Most likely there is a faster method if is_shape_unique
            #
            self._derived_properties["gridinterest"] = geopandas.overlay(self.gridin.geodataframe,
                                                                         self.gridout.geodataframe,
                                                                         how='intersection')
            self.gridinterest["area"] = self.gridinterest.apply(localdef_get_area, axis=1)
        else:
            warnings.warn("Cannot measure gridinterest, because gridin and/or gridout is/are None")

    # -------------- #
    #  Measurement   #
    # -------------- #
    def project_data(self, data, as_serie=True, use="sum"):
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
        projected_data = self._project_data_(self._parse_data_(data), use=use)

        if as_serie:            
            return projected_data

        projected_data_array = np.zeros( len(self.gridout.geodataframe) )
        projected_data_array[projected_data.index.values] = projected_data.values
        return projected_data_array
        
    def _project_data_(self, data, use="sum"):
        """ """
        self.gridinterest["_tmp"] = data[ self.gridin.geodataframe.loc[ self.gridinterest["id_1"]].index
                                        ] * self.gridinterest["area"]
        
        return getattr(self.gridinterest.groupby("id_2")["_tmp"],use)()

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
    
    PROPERTIES         = ["pixels", "shape"]
    SIDE_PROPERTIES    = ["indexes"]
    DERIVED_PROPERTIES = ["vertices","geodataframe", "triangulation"]
    
    def __init__(self, pixels=None, shape=UNIT_SQUARE, indexes=None):
        """ """
        if pixels is not None:
            self.set_pixels(pixels,shape=shape)

        if indexes is not None:
            self.set_indexes(indexes)
            
    # =================== #
    #   Methods           #
    # =================== #
    @classmethod
    def from_stamps(cls, stamp, origin=[0,0]):
        """ stamps are 2d array, something you could to ax.imshow(stamps) 
        data will be stored as 'data' in the grid's dataframe
        """
        this = get_simple_grid(*np.shape(stamp), shift_origin=origin)
        this.add_data(np.ravel(stamp), "data")
        return this

    @classmethod
    def from_vertices(cls, vertices, indexes=None):
        """ directly provide the vertices

        Parameters:
        -----------
        vertices: [list of array or dictionary]
            The vertices of all the grid entries. 
            Could have two format:
            - list of array: [[vert_1],[vert_2],....], then you may want to provide indexes
            - dictionary: {id_1:vert_1,id_2: vert_2, ...} 
            if a dictionary is provided, the indexes will be set by the vertices.

        indexes: [list or None] -optional-
            (Ignored if vertices is a dict)
            If you provide vertices as a list of vertices, you can provide the
            indexes of each of the vertices. 
            -> if None, then indexes = np.arange(len(vertices))

        Returns
        -------
        Grid
        """
        this = cls()
        if type(vertices) is dict:
            indexes, vertices = list(vertices.keys()), list(vertices.values())
            
        this.set_vertices(vertices)
        if indexes is not None:
            this.set_indexes(indexes)
            
        return this
        
    @classmethod
    def set_from(cls, datainput):
        """ Creates a new Grid objects from the given input data:
        
        Parameters
        ----------
        datainput: [geopandas.geodataframe.GeoDataFrame or ndarray]
            this could either be a:
            - geodataframe (and this calls self.set_geodataframe)
            - geoSeries
            - ndarray: if 3-shaped, this calls set_vertices ;
                       if 2-shaped, this calls set_pixels.
        
        Returns
        -------
        Grid
        """
        this = cls()
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
        """ provide the indexes associated to each pixels

        Parameters
        ----------
        indexes: [ndarray]
           indexes associated to the pixels. 
           This should have the length equal to th number of pixels (if any).

        update: [bool] -optional-
           should the geodataframe be updated ?
           [use True if you are not sure]

        Returns
        -------
        Void
        """
        if self.pixels is not None and len(indexes) != self.npixels:
            raise AssertionError("not the same number of indexes as the number of pixels")
        self._side_properties["indexes"]  = indexes
        if update:
            self._update_geodataframe_()
        
    def set_pixels(self, pixels, shape=None, update=True):
        """ provide the pixels. 
        
        Pixels define the position up on which the geometries are defined.
        NB: vertices = pixels+shape
        """
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

        try:
            pixels = np.mean(vertices, axis=1)
        except:
            # Means vertices have different size.
            self._derived_properties["vertices"] = vertices
            pixels = np.asarray([np.mean(v_, axis=0) for v_ in vertices])
            self.set_pixels(pixels, None, **kwargs)
            return

        self._derived_properties["vertices"] = np.asarray(vertices)
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
        
        vertices = geodataframe["geometry"].apply(get_verts).values
        self.set_vertices(vertices, update=False) # don't update the geodataframe

    
        
        
    # --------- #
    #  UPDATE   #
    # --------- #
    def _update_geodataframe_(self):
        """ """
        dataseries = self.get_geoseries()
        x,y     = self.pixels.T
        self._derived_properties["geodataframe"] = \
          geopandas.GeoDataFrame({'geometry': dataseries,
                                        'id':self.indexes,
                                  'x':x,'y':y})
        
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
        """ build a new geodataframe and returns it. """
        import geopandas
        return geopandas.GeoSeries([geometry.Polygon(v) for v in self.vertices])

    def get_triangulation_grid(self):
        """ Returns a grid of triangulation. """
        return Grid.set_from( np.concatenate(self.triangulation, axis=0) )

    def get_pixels_in(self, polygon, invert=False):
        """ checks if the centroid of the pixel is in or out the given shapely polygon.
        
        Parameters
        ----------
        polygon: [shapely.geometry]
            reference polygon 

        invert: [bool] -optional-
            Get the pixel inside the polygon [invert=False] or outsite [invert=True]

        Returns
        -------
        list of pixels and boolean mask
        """
        from shapely import vectorized
        flagin = vectorized.contains(polygon, *self.pixels.T)
        if invert:
            flagin = ~flagin
        return self.pixels[flagin], flagin
    
        
    # --------- #
    # Project   #
    # --------- #
    def project_to(self, othergrid, column="*", asgrid=True, use="sum"):
        """ project data in the given grid

        Parameters
        ----------
        othergrid: [Grid]
            New grid where data should be projected to
            
        column: [str/None/list of] -optional-
            Which data should be projected ? 
            If None or '*' all the non-structural columns will be
            (structural columns are 'geometry', 'id', 'x', 'y')

        asgrid: [bool] -optional-
            Should this return a new Grid (actually same object as othergrid)
            or a dict [asgrid=False]?

        Returns
        -------
        Grid or dict (see asgrid)
        """
        gproj = GridProjector(self, othergrid)
        if column is None or column in ["*","all"]:
            column = [k for k in self.geodataframe if k not in ['geometry', 'id', 'x', 'y']]
        datas = {k:gproj.project_data(k, use=use) for k in column}
        if not asgrid:
            return datas
        # building and setting the new grid
        gout = othergrid.__class__.set_from(othergrid.geodataframe)
        for k in column:
            gout.add_data(datas[k],k)
        return gout
    
    def project_to_wcs(self, wcs_, asgrid=True, **kwargs):
        """ provide an astropy.wcs.WCS and this will project 
        the current grid into it (assuming grid's vertices coordinates are in pixels) 

        Parameters
        ----------
        wcs_: [astropy.wcs.WCS]
            The world coordinate solution

        asgrid: [bool] -optional-
            Should this return a load Grid object or an array of vertices (in degree)
            
        **kwargs goes to wcs_.all_pix2world

        Returns
        -------
        Grid or array (see asgrid)

        """
        verts             = self.vertices
        verts_shape       = np.shape(verts)
        flatten_verts     = np.concatenate(verts, axis=0)
        # 
        flatten_verts_wcs = np.asarray(wcs_.all_pix2world(flatten_verts[:,0],
                                                          flatten_verts[:,1], 0,
                                                          **kwargs)).T
        
        #
        verts_wcs = flatten_verts_wcs.reshape(verts_shape)
        if not asgrid:
            return verts_wcs

        g_wcs = Grid.set_from(verts_wcs)
        g_wcs.geodataframe["x_pix"],g_wcs.geodataframe["y_pix"] = self.pixels.T
        return g_wcs


    def evaluate(self, func, vectorized=True):
        """ Evaluate the given function throughout the grid.
        This evulation is using polynome triangulation to integrate the 
        given function inside the polyname using triangle integration.
        
        -> dependency: the integration is made using quadpy.

        Examples:
        # Remark the np.stack(x, axis=-1). 
        # This is mandatory since integration is going to send 
        #  x = [ [[....],[...]], [[....],[...]], ... ] for triangles
        ```python
        def get_2dgauss(x, mu=[4,4], cov=[[1,0],[0,2]]):
            """ """
            return stats.multivariate_normal.pdf(np.stack(x, axis=-1), mean=mu, cov=cov)
        ```

        """
        try:
            import quadpy
        except ImportError:
            raise ImportError("Integration is made using quadpy. pip install quadpy")

        # Is Triangulation made ?
        if self._derived_properties["triangulation"] is None:
            warnings.warn("triangles not defined: deriving triangulation.")
            self.derive_triangulation()
            
        # Let's get the triangles
        trs = np.stack(self.triangulation)
        shape_trs = np.shape(trs)
        
        if len(shape_trs)==4 and vectorized: # All Polygon have the same topology (same amount of vertices)
            tr_flat = np.stack(np.concatenate(trs, axis=0), axis=-2)
            val = quadpy.triangle.strang_fix_cowper_09().integrate(func,tr_flat).reshape(shape_trs[:2])
        else:
            val = np.asarray([quadpy.triangle.strang_fix_cowper_09().integrate(func,np.stack(t_, axis=-2)) for t_ in trs])

        return np.sum(val, axis=1)
        
    def derive_triangulation(self, fast_unique=True):
        """ """
        
        def triangulate(geom):
            """ Return triangulate format that quadpy likes """
            from shapely import ops
            triangles = ops.triangulate(geom)
            return np.stack([np.asarray(t.exterior.coords.xy).T[:-1] for t in triangles])

        if not self.is_shape_unique or not fast_unique:
            self._derived_properties["triangulation"] = self.geodataframe["geometry"].apply(triangulate)
        else:
            self._derived_properties["triangulation"] = self.pixels[:,None,None] + triangulate(geometry.Polygon(self.shape))

    # --------- #
    #  PLOTTER  #
    # --------- #
    def show(self, column=None, ax=None, edgecolor="0.7", facecolor="None", **kwargs):
        """ """
        if column is not None:
            facecolor=None
        return self.geodataframe.plot(column, ax=ax,facecolor=facecolor,
                                          edgecolor=edgecolor, **kwargs)
    
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
    
    @property
    def triangulation(self):
        """ Triangulation of the vertices. Based on Delaunay tesselation, see shapely.ops.triangulate """
        if self._derived_properties["triangulation"] is None:
            self.derive_triangulation()
        return self._derived_properties["triangulation"]
