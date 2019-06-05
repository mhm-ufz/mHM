#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Author
------
David Schaefer

Purpose
-------
A sanitizing layer for the netCDF4 library. Adds a number of convenince methods
and aims for a cleaner user interface. All classes avaliable are children of their
netCDF4 counterparts.
"""

import uuid
from netCDF4 import Dataset, Group, Dimension, Variable
from netCDF4 import (chartostring, date2index, date2num, getlibversion,
                     num2date, stringtoarr, stringtochar)
from collections import OrderedDict


def _tupelize(arg):
    if isinstance(arg, str):
        return (arg,)
    try:
        return tuple(arg)
    except TypeError:
        return (arg,)


def copyGroup(ncin, group, skipdims=None, skipgroups=None, skipvars=None,
              skipattrs=None, fixdims=False, vardata=False, varparams=None):
    """
    Arguments
    ---------
    ncin                  : Instance of an object with a createGroup method
                            (i.e. NcDataset, NcGroup)
    group                 : Instance of an object with dimensions/variables/attributes/groups attributes
                            (i.e. NcDataset, NcGroup)
    skipdims (optional)   : string or list/tuple of strings
                            Name(s) of dimension(s) to skip
    skipgroups (optional) : string or list/tuple of strings
                            Name(s) of group(s) to skip
    skipvars (optional)   : string or list/tuple of strings
                            Name(s) of variable(s) to skip
    skipattrs (optional)  : string or list/tuple of strings
                            Name(s) of attribute(s) to skip
    fixdims (optional)    : convert unlimited dimensions to fixed size dimensions
    vardata (optional)    : boolean, copy variable data
    varparams(optional)   : dict, variable paramaters will be passed to
                            createVariable (i.e. zlib, complevel, chunksizes, ...)

    Return
    ------
    NcGroup

    Purpose
    -------
    Copy the given group to ncin
    """
    out = ncin.createGroup(group.name)
    out.set_fill_off()
    out.copyDimensions(group.dimensions, skip=skipdims, fix=fixdims)
    out.copyVariables(
        group.variables, skip=skipvars, data=vardata, varparams=varparams
    )
    out.copyAttributes(group.attributes, skipattrs)
    out.copyGroups(group.groups, skipgroups)
    return out


def copyDataset(ncin, group, skipdims=None, skipgroups=None, skipvars=None,
                skipattrs=None, fixdims=False, vardata=False, varparams=None):
    """
    Arguments
    ---------
    ncin                  : Instance of an object with a createGroup method
                            (i.e. NcDataset, NcGroup)
    group                 : Instance of an object with dimensions/variables/attributes/groups attributes
                            (i.e. NcDataset, NcGroup)
    skipdims (optional)   : string or list/tuple of strings
                            Name(s) of dimension(s) to skip
    skipgroups (optional) : string or list/tuple of strings
                            Name(s) of group(s) to skip
    skipvars (optinal)    : string or list/tuple of strings
                            Name(s) of variable(s) to skip
    skipattrs (optinal)   : string or list/tuple of strings
                            Name(s) of attribute(s) to skip
    fixdims (optional)    : convert unlimited dimensions to fixed size dimensions
    vardata (optional)    : boolean, copy variable data
    varparams(optional)   : dict, variable paramaters will be passed to
                            createVariable (i.e. zlib, complevel, chunksizes, ...)


    Return
    ------
    NcDataset/NcGroup

    Purpose
    -------
    Copy the content of given group to ncin
    """
    ncin.set_fill_off()
    ncin.copyDimensions(group.dimensions, skip=skipdims, fix=fixdims)
    ncin.copyVariables(
        group.variables, skip=skipvars, data=vardata, varparams=varparams
    )
    ncin.copyAttributes(group.attributes, skipattrs)
    ncin.copyGroups(group.groups, skipgroups)
    return ncin


def copyGroups(ncin, groups, skip=None):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createGroup method
                      (i.e. NcDataset, NcGroup)
    groups          : Dictionary
                      key   : group name (string)
                      value : instance of an object with dimensions/variables/attributes/groups attributes
    skip (optional) : string or list/tuple of strings
                      Name(s) of group(s) to skip

    Return
    ------
    None

    Purpose
    -------
    Copy the given groups to ncin
    """
    for g in groups.values():
        if g.name not in _tupelize(skip):
            ncin.copyGroup(g)


def copyDimension(ncin, dim, fix=False, fail=True):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createDimension method
                      (i.e. NcDataset, NcGroup)
    group           : Instance of NcDimension
    fix (optional)  : convert an unlimited dimension to a fixed size dimension
    fail (Optional[bool]): raise an exception if dimension already exists

    Return
    ------
    netCDF4.Dimension

    Purpose
    -------
    Copy the given dimension to ncin
    """
    length = None if dim.isunlimited() and not fix else len(dim)
    try:
        return ncin.createDimension(dim.name, length)
    except Exception:
        if fail:
            raise
        return ncin.dimensions[dim.name]


def copyDimensions(ncin, dimensions, skip=None, fix=False, fail=True):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createDimension method
                      (i.e. NcDataset, NcGroup)
    dimension       : Dictionary
                      key   : dimension name (string)
                      value : instance of NcDimension
    skip (optional) : string or list/tuple of strings
                      Name(s) of dimension(s) to skip
    fix (optional)  : convert unlimited dimensions to fixed size dimensions
    fail (Optional[bool]): raise exception if dimension already exists

    Return
    ------
    None

    Purpose
    -------
    Copy the given dimensions to ncin
    """
    for d in dimensions.values():
        if d.name not in _tupelize(skip):
            ncin.copyDimension(d, fix=fix, fail=fail)


def copyAttributes(ncin, attributes, skip=None):
    """
    Arguments
    ---------
    ncin              : Instance of an object with a createAttribute method
                        (i.e. NcDataset, NcGroup, NcVariable)
    attributes        : Dictionary
                        key   : string
                        value : string/any numeric type
    skip (optional)   : string or list/tuple of strings
                        Name(s) of attribute(s) to skip

    Return
    ------
    None

    Purpose
    -------
    Copy the given attributes to ncin
    """
    for k, v in attributes.items():
        if k not in _tupelize(skip):
            if k == "missing_value":
                try:
                    v = ncin.dtype.type(v)
                except Exception:
                    pass
            ncin.createAttribute(k, v)


def copyVariable(ncin, var, data=True, dims=False, fail=True, **kwargs):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a createCopy method
                      (i.e. NcDataset, NcGroup, NcVariable)
    var             : Instance of NcVariable
    data (optional) : boolean, copy variable data
    dims (Optional[bool]): copy missing dimensions
    fail (Optional[bool]): raise an exception if variable exists, ignored at the moment
    kwargs          : will be passed to createVariable. Allows to set
                      parameters like chunksizes, deflate_level, ...

    Return
    ------
    NcVariable

    Purpose
    -------
    Copy the given variables to ncin. Copy the data if data=True
    """
    invardef = var.definition
    
    if data is not True:
        invardef["chunksizes"] = None
    invardef.update(kwargs)

    if dims:
        nc = var.parent
        try:
            shape = data.shape
        except AttributeError:
            shape = var.shape
        for name, length in zip(var.dimensions, shape):
            l = None if nc.dimensions[name].isunlimited() else length
            ncin.createDimension(name, l, fail=fail)

    vname = invardef.pop("name")
    try:
        invar = ncin.createVariable(vname, invardef.pop("dtype"), **invardef)
    except RuntimeError:
        if fail:
            raise
        invar = ncin.variables[vname]

    invar.copyAttributes(var.attributes)
    if data is True and var.shape:
        invar[:] = var[:]
    elif data is not False:
        # i.e. if an array is given
        invar[:] = data
    return invar


def copyVariables(ncin, variables, skip=None, data=True, dims=False, fail=True, varparams=None):
    """
    Arguments
    ---------
    ncin                    : Instance of an object with a createCopy method
                              (i.e. NcDataset, NcGroup, NcVariable)
    variables               : Dictionary
                              key   : variables name (string)
                              value : instance of NcVariable
    skip (optional)         : string or list/tuple of strings
                              Name(s) of variable(s) to skip
    data (optional)         : boolean
    dims (Optional[bool])   : copy nonexisting dimensions
    fail (Optional[bool])   : raise an exception if a variable already exists
    varparams(optional)     : dict, variable paramaters will be passed to
                              createVariable (i.e. zlib, complevel, chunksizes, ...)

    Return
    ------
    NcVariable

    Purpose
    -------
    Copy the given variables to ncin. Copy the data if data=True
    """

    if varparams is None:
        varparams = dict()
    for v in variables.values():
        if v.name not in _tupelize(skip):
            ncin.copyVariable(v, data, dims, fail, **varparams)


def createDimensions(ncin, dim_dict, fail=True):
    for name, length in dim_dict.items():
        ncin.createDimension(name, length, fail=fail)


def getVariableDefinition(ncvar):
    out = ncvar.filters() if ncvar.filters() else {}
    out.update({
        "name"       : ncvar.name,
        "dtype"      : ncvar.dtype,
        "dimensions" : ncvar.dimensions,
        "chunksizes" : ncvar.chunking() if not isinstance(ncvar.chunking(), str) else None,
        "fill_value" : getattr(ncvar, "_FillValue", None),
    })
    return out


def getDates(ncin, timesteps=None, timevar="time", units=None, calendar=None):
    """
    Arguments
    ---------
    ncin                 : Instance of an object holding variables (NcDataset/NcGroup)
    timesteps (optional) : list/tuple/nd.array of Numerical values.
                           The time_steps to return dates for. If not given the content of the
                           entire time variable will be returned.
    units (optional)     : string
                           time units following the CF Conventions. Needs to be given, if not
                           available as an attribute of the time variable.
    calendar (optional)  : string
                           calendar name following the CF conventions. Needs to be given, if not
                           available as an attribute of the time variable.

    Return
    ------
    List of datetime objects

    Purpose
    -------
    Return datetime objects associated to the time variable of ncin
    """
    var = ncin.variables[timevar]
    if not units:
        try:
            units = var.units
        except AttributeError:
            raise AttributeError(
                "Time variable does not specify an units attribute! Pass as argument."
            )

    if not calendar:
        try:
            calendar = var.calendar
        except AttributeError:
            calendar = "standard"

    if not timesteps:
        timesteps = var[:]

    dates = num2date(timesteps, units, calendar)

    try:
        return [d.date() for d in dates]
    except AttributeError:
        return dates


def setFillValue(ncin, value):
    """
    Arguments
    ---------
    ncin  : Instance of an object with a _FillValue attribute
            (i.e. NcVariable)

    Return
    ------
    Numeric

    Purpose
    -------
    Return the value of the attribute _FillValue
    """
    ncin.setncattr("_FillValue", value)


def getFillValue(ncin):
    """
    Arguments
    ---------
    ncin  : Instance of an object with a _FillValue attribute
            (i.e. NcVariable)

    Return
    ------
    Numeric

    Purpose
    -------
    Return the value of the attribute _FillValue
    """
    try:
        return ncin.getncattr("_FillValue")
    except AttributeError:
        return None


def setAttribute(ncin, name, value):
    """
    Arguments
    ---------
    ncin  : Instance of an object with a setncatts method
            (i.e. NcDataset/NcGroup/NcVariable)
    name  : string
    value : string or any numeric type

    Return
    ------
    None

    Purpose
    -------
    Set/Write the attribute given as name, value
    """
    ncin.setncattr(name, value)


def setAttributes(ncin, attdict):
    """
    Arguments
    ---------
    ncin    : Instance of an object with a setncatts method
              (i.e. NcDataset/NcGroup/NcVariable)
    attdict : dictionary
              key: attribute name (string)
              value: attribute value (string or any numeric type)

    Return
    ------
    None

    Purpose
    -------
    Set/Write the attributes given in attdict
    """
    ncin.setncatts(attdict)


def filterVariables(ncin, dims=None, ndim=None):
    """
    Arguments
    ---------
    ncin            : Instance of an object with a variables attribute
                      (i.e. NcDataset/NcGroup)
    dims (optional) : tuple/list of dimension strings
    ndim (optional) : int number of dimensions

    Return
    ------
    OrderedDict:
        key   : variable name (string)
        value : NcVariable instance

    Purpose
    -------
    Return all Variables that are based on the dimension(s) given in dims
    and/or have ndims dimensions.
    """
    out = OrderedDict()
    dims = set(dims or {})

    for v in ncin.variables.values():
        if dims.issubset(set(v.dimensions)):
            if ndim:
                if ndim == len(v.dimensions):
                    out[v.name] = v
            else:
                out[v.name] = v

    return out


def filterDimensions(ncin, lengths):
    """
    Arguments
    ---------
    lengths : integer, tuple/list of integers

    Return
    ------
    OrderedDict:
        key   : dimension name (string)
        value : NcDimension instance

    Purpose
    -------
    Return all Dimensions with a length given in the argument lengths.
    """
    try:
        lengths[0]
    except TypeError:
        lengths = (lengths,)

    return OrderedDict(
        [(d.name, d) for d in ncin.dimensions.values() if len(d) in lengths]
    )


def getGroups(ncin):
    out = OrderedDict()
    for g in getattr(ncin, "groups").values():
        out[g.name] = NcGroup(ncin, g.name, id=g._grpid)
    return out


def getVariables(ncin):
    out = OrderedDict()
    for v in getattr(ncin, "variables").values():
        out[v.name] = NcVariable(
            ncin, v.name, v.dtype, v.dimensions, id=v._varid
        )
    return out


def getDimensions(ncin):
    
    out = OrderedDict()
    for dim in getattr(ncin, "dimensions").values():
        out[dim.name] = NcDimension(ncin, dim.name, id=dim._dimid)
    return out


def getAttributes(ncin):
    out = OrderedDict()
    for k in ncin.ncattrs():
        if not k.startswith("_"):
            out[k] = ncin.getncattr(k)
    return out


def getParent(ncin):
    return ncin._grp


def attributeSetter(ncin, name, value):
    ncin.__dict__[name] = value


def attributeGetter(ncin, name):
    try:
        return ncin.__dict__[name]
    except KeyError:
        try:
            return getattr(super(ncin.__class__, ncin), name)
        except KeyError:
            raise AttributeError(
                "'{:}' object has no attribute '{:}'".format(ncin.__class__, name)
            )

def createGroup(ncin, name):
    grp = NcGroup(ncin, name)
    ncin.groups[name] = grp
    return grp

def createVariable(ncin, *args, **kwargs):
    var = NcVariable(ncin, *args, **kwargs)
    ncin.variables[var.name] = var
    return var

def createDimension(ncin, name, length, fail=True):
    try:
        dim = NcDimension(ncin, name, length)
        ncin.dimensions[dim.name] = dim
    except Exception:
        if fail:
            raise
        dim = ncin.dimensions[name]
    return dim


class NcDataset(Dataset):
    def __init__(
            self,
            filename = None,
            mode     = "r",
            clobber  = True,
            diskless = False,
            persist  = False,
            weakref  = False,
            format   = "NETCDF4",
    ):

        if filename is None:
            # in memory dataset
            filename = str(uuid.uuid4())
            mode = "w"
            diskless = True

        super(NcDataset, self).__init__(
            filename = filename,
            mode     = mode,
            clobber  = clobber,
            diskless = diskless,
            persist  = persist,
            weakref  = weakref,
            format   = format,
        )

        self.fname = filename
        for k, v in zip(self.groups, getGroups(self).values()):
            self.groups[k] = v
        for k, v in zip(self.dimensions, getDimensions(self).values()):
            self.dimensions[k] = v
        for k, v in zip(self.variables, getVariables(self).values()):
            self.variables[k] = v

    def tofile(self, fname):
        # preserve dataset options
        with NcDataset(fname, "w") as out:
            out.copyDataset(self, vardata=True)

    # def __enter__(self):
    #     return self

    # def __exit__(self, *args, **kwargs):
    #     self.close()

    copyDataset      = copyDataset
    copyDimension    = copyDimension
    copyDimensions   = copyDimensions
    copyAttributes   = copyAttributes
    copyVariable     = copyVariable
    copyVariables    = copyVariables
    copyGroup        = copyGroup
    copyGroups       = copyGroups
    createAttribute  = setAttribute
    createAttributes = setAttributes
    createDimensions = createDimensions
    createVariable   = createVariable
    createGroup      = createGroup
    createDimension  = createDimension
    filterVariables  = filterVariables
    filterDimensions = filterDimensions
    getDates         = getDates
    attributes       = property(fget=getAttributes)
    # restore a "normal" attribute access behaviour
    __setattr__      = attributeSetter
    __getattr__      = attributeGetter


class NcGroup(Group):
    def __init__(self, *args, **kwargs):
        super(NcGroup,self).__init__(*args, **kwargs)
        for k, v in zip(self.groups, getGroups(self).values()):
            self.groups[k] = v
        for k, v in zip(self.dimensions, getDimensions(self).values()):
            self.dimensions[k] = v
        for k,v in zip(self.variables, getVariables(self).values()):
            self.variables[k] = v

    copyDimension    = copyDimension
    copyDimensions   = copyDimensions
    copyAttributes   = copyAttributes
    copyVariable     = copyVariable
    copyVariables    = copyVariables
    copyGroup        = copyGroup
    copyGroups       = copyGroups
    createAttribute  = setAttribute
    createAttributes = setAttributes
    createDimensions = createDimensions
    createVariable   = createVariable
    createGroup      = createGroup
    createDimension  = createDimension
    filterVariables  = filterVariables
    filterDimensions = filterDimensions
    getDates         = getDates
    attributes       = property(fget=getAttributes)
    parent           = property(fget=getParent)
    # restore a "normal" attribute access behaviour
    __setattr__      = attributeSetter
    __getattr__      = attributeGetter


class NcVariable(Variable):
    def __init__(self, *args, **kwargs):
        super(NcVariable, self).__init__(*args, **kwargs)

    copyAttributes   = copyAttributes
    createAttribute  = setAttribute
    createAttributes = setAttributes
    attributes       = property(fget=getAttributes)
    definition       = property(fget=getVariableDefinition)
    fill_value       = property(fget=getFillValue, fset=setFillValue)
    parent           = property(fget=getParent)
    # restore a "normal" attribute access behaviour
    __setattr__      = attributeSetter
    __getattr__      = attributeGetter


# Just to be consistent...
class NcDimension(Dimension):
    def __init__(self, *args, **kwargs):
        super(NcDimension, self).__init__(*args, **kwargs)
    
    parent = property(fget=getParent)
