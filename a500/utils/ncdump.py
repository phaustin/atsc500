"""
  recursively dump netcdf file metadata -- taken from:

  http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
"""

import argparse
from netCDF4 import Dataset
import textwrap

def ncdump(nc_fid, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''

    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)

    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    
    if verb:
        groups=list(nc_fid.groups.items())
        if not groups:
            group_name="root"
            groups = [(group_name,nc_fid)]
        for group_name, group in groups:
            print(f"NetCDF variable information for group {group_name}:")
            nc_vars = [var for var in group.variables]  # list of nc variables
            for var in nc_vars:
                if var not in nc_dims:
                    print('\tName:', var)
                    print("\t\tdimensions:", group.variables[var].dimensions)
                    print("\t\tsize:", group.variables[var].size)
                    print_ncattr(var)
    return nc_attrs, nc_dims, nc_vars

def make_parser():
    linebreaks = argparse.RawTextHelpFormatter
    descrip = textwrap.dedent(__doc__)
    parser = argparse.ArgumentParser(formatter_class=linebreaks,
                                     description=descrip)
    parser.add_argument('ncfile', type=str, help='path to netcdf file to dump')
    return parser

def main(args=None):
    parser = make_parser()
    args=parser.parse_args(args)
    with Dataset(args.ncfile) as nc_in:
        ncdump(nc_in)
    

if __name__ == "__main__":
    main()
    
