#!/usr/local/other/python/GEOSpyD/2019.03_py3.7/2019-04-22/bin/python
'''
Script to create uniform CH4 offset file that can be read by GEOS to
add a universal offset to CH4 boundary conditions.
The offset value can be either passed as input argument or calculated
online as the difference between the target year and the reference year.
In case of the latter, the values from the following source are used to
determine the CH4 increase: 
- https://gml.noaa.gov/webdata/ccgg/trends/ch4/ch4_mm_gl.txt

EXAMPLES:
1. Create CH4 growth files relative to 2021 passing explicit values:

Create y/y offset file for January 2022 relative to 2021, using prescribed offset
python create_ch4_offset_file.py -y 2022 -m 1 -v 17.97e-9 -o 'ch4_offset_to_2021_%Y%m.nc'

Create y/y offset file for January 2023 relative to 2022, using prescribed offset
python create_ch4_offset_file.py -y 2023 -m 1 -v 14.62e-9 -o 'ch4_offset_to_2022_%Y%m.nc'

Create y/y offset file for January 2023 relative to 2021 by adding y/y offset to 2022/2021 file
python create_ch4_offset_file.py -y 2023 -m 1 -v 14.62e-9 -i 'ch4_offset_to_2021_2022%m.nc' -o 'ch4_offset_to_2021_%Y%m.nc'

2. Create CH4 growth file relative to 2021 using online computed offset:
python create_ch4_offset_file.py -y 2023 -m 1 -r 2021 -o 'ch4_offset_to_2021_%Y%m.nc'

HISTORY: 
20240118 - christoph.a.keller@nasa.gov - initial version
'''
import sys
import argparse
import logging
import datetime as dt
import time
import numpy as np
import xarray as xr
import pandas as pd

# NOAA CH4 data
URL="https://gml.noaa.gov/webdata/ccgg/trends/ch4/ch4_mm_gl.txt"

def main(args):
    '''
    Create a netCDF file with the given offset values. 
    '''
    log = logging.getLogger(__name__)
    # output time
    otime = dt.datetime(args.year,args.month,1)
    # offset value: read from online table if not provided
    offset = np.nan 
    if args.value is None:
       log.info('Reading CH4 trends from {}'.format(URL))
       dat = pd.read_csv(URL,comment="#",header=None,sep='\s+')
       # get reference values in ppb
       refconc = dat.loc[(dat[0]==args.refyear)&(dat[1]==args.month),3].values[0]
       targetconc = dat.loc[(dat[0]==args.year)&(dat[1]==args.month),3].values[0]
       offset = (targetconc-refconc)*1.0e-9
       log.info("Calculated offset from URL: {}".format(URL))
       log.info("reference concentration ({}-{}): {}".format(args.refyear,args.month,refconc))
       log.info("target concentration ({}-{}): {}".format(args.year,args.month,targetconc))
       log.info("offset (mol/mol): {}".format(offset))
    else:
       offset = args.value

    # define lat/lon coordinates and output array with offsets. 
    if args.ifile is not None:
       # Read input file and inherit coordinates from it. Add offset to values in that file
       ifile = otime.strftime(args.ifile) 
       log.info('reading {}'.format(ifile))
       dsi = xr.open_dataset(ifile)
       lons = dsi.lon.values
       lats = dsi.lat.values
       assert args.varname in dsi,"variable {} not found in {}".format(args.varname,args.ifile)
       assert len(dsi[args.varname].shape)==3,"variable {} must have 3 dimensions (time,lat,lon)".format(args.varname)
       arr  = dsi[args.varname].values[:,:,:] + offset 
    else:
       # Create from scratch 
       lons = np.arange(-180.,180.,2.5)
       lats = np.arange(-90.,90.1,2.) 
       arr = np.zeros((1,len(lats),len(lons)))
       arr[:] = offset

    # Create output dataset
    ds = xr.Dataset()
    ds[args.varname] = (('time','lat','lon'),arr)
    ds[args.varname].attrs =  {'standard_name':'CH4_offset','long_name':'CH4_offset','units':"mol/mol"}
    ds.coords['lat'] = (('lat'),lats)
    ds['lat'].attrs =  {'standard_name':'latitude','long_name':'latitude','units':'degrees_north'}
    ds.coords['lon'] = (('lon'),lons)
    ds['lon'].attrs =  {'standard_name':'longitude','long_name':'longitude','units':'degrees_east'}
    tunit = otime.strftime('days since %Y-%m-%d %H:%M:%S')
    ds.coords['time'] = (('time'),np.zeros((1,)))
    ds['time'].attrs =  {'standard_name':'time','long_name':'time','units':tunit,'calendar':'standard'}
    ofile = otime.strftime(args.ofile)
    ds.to_netcdf(ofile)
    log.info('file written to {}'.format(ofile))
    ds.close()
    return


def parse_args():
    p = argparse.ArgumentParser(description='Undef certain variables')
    p.add_argument('-y', '--year',type=int,help='output year',default=2023)
    p.add_argument('-m', '--month',type=int,help='output month',default=1)
    p.add_argument('-v', '--value',type=float,help='CH4 offset, in mol/mol dry',default=None)
    p.add_argument('-r', '--refyear',type=int,help='reference year',default=2021)
    p.add_argument('-i', '--ifile',type=str,help='input file, will add offset to values in that file',default=None)
    p.add_argument('-o', '--ofile',type=str,help='output file',default="ch4_offset_to_2021_x144_y91_%Y%m.nc")
    p.add_argument('-n', '--varname',type=str,help='variable name for offset',default="CH4_offset")
    return p.parse_args()


if __name__ == '__main__':
    log = logging.getLogger()
    log.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    log.addHandler(handler)
    main(parse_args())

