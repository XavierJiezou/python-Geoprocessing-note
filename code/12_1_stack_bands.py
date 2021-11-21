from osgeo import gdal
import numpy as np


def stack_bands(filename):
    ds = gdal.Open(filename)
    return np.dstack([ds.GetRasterBand(i).ReadAsArray() for i in range(1, ds.RasterCount+1)])


if __name__ == '__main__':
    print(stack_bands('data/Landsat/Washington/nat_color.tif').shape) # (8023, 8849, 3)
