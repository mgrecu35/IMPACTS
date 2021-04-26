fname='https://goldsmr5.gesdisc.eosdis.nasa.gov/opendap/MERRA2/M2I3NPASM.5.12.4/2020/02/MERRA2_400.inst3_3d_asm_Np.20200227.nc4'

from netCDF4 import Dataset

fh=Dataset(fname)
