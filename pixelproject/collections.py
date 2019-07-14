#! /usr/bin/env python
#
import pandas
import numpy as np

from propobject import BaseObject


class GridCollection( BaseObject ):

    PROPERTIES = ["grids"]
    DERIVED_PROPERTIES = ["geodataframe"]
    
    def __init__(self):
        """ """
        

    def _update_geodataframe_(self):
        """ """
        def new_name(l):
            return "%s_%s"%(l["cid"],l["id"])
        
        self._derived_properties["geodataframe"] = pandas.concat([g.geodataframe for g in self.grids.values()], ignore_index=True)
        self.geodataframe["sid"] = self.geodataframe.apply(new_name, axis=1)
        self.geodataframe["id"]  = np.arange( len(self.geodataframe) )
        
    # =================== #
    #   Methods           #
    # =================== #
    # --------- #
    #  SETTER   #
    # --------- #
    def add_grid(self, grid, update=True, noclean_up=False):
        """ """
        id_ = self.ngrids
        grid.geodataframe["cid"] = id_
        self.grids[id_] = grid
        # clean up
        if not noclean_up:
            self._derived_properties["geodataframe"] = None
            
        if update:
            self._update_geodataframe_()
            
    # =================== #
    #   Properties        #
    # =================== #
    @property
    def grids(self):
        """ """
        if self._properties["grids"] is None:
            self._properties["grids"] = {}
        return self._properties["grids"]

    @property
    def ngrids(self):
        """ """
        return len(self.grids)
    
    @property
    def geodataframe(self):
        """ """
        if self._derived_properties["geodataframe"] is None:
            self._update_geodataframe_()
        return self._derived_properties["geodataframe"]
