import os
from osgeo import gdal


def three_bands_stacked(band1_fn: str, band2_fn: str, band3_fn: str, num_band: int, out_path: str):
    """Three bands stacked into an RGB image.
    Args:
        band1_fn (str): File name of red band.
        band2_fn (str): File name of green band.
        band3_fn (str): File name of blue band.
        num_band (int): Number of bands.
        out_path (str): Path to output file.
    Refer:
        https://github.com/cgarrard/osgeopy-code/blob/master/Chapter9/listing9_1.py
    """
    # Open band 1.
    in_ds = gdal.Open(band1_fn)
    in_band = in_ds.GetRasterBand(1)

    # Create a 3-band GeoTIFF with the same dimensions, data type, projection,
    # and georeferencing info as band 1. This will be overwritten if it exists.
    gtiff_driver = gdal.GetDriverByName('GTiff')
    dst_ds = gtiff_driver.Create(  # Refer: https://gdal.org/tutorials/raster_api_tut.html#using-create
        out_path,
        xsize=in_band.XSize,
        ysize=in_band.YSize,
        bands=num_band,
        eType=in_band.DataType
    )
    dst_ds.SetProjection(in_ds.GetProjection())
    dst_ds.SetGeoTransform(in_ds.GetGeoTransform())

    # Copy data from band 1, 2, 3 into the output image.
    for index, name in enumerate([band1_fn, band2_fn, band3_fn]):
        dst_ds.GetRasterBand(index+1).WriteArray(gdal.Open(name).ReadAsArray())

    # Compute statistics on each output band.
    dst_ds.FlushCache()
    for i in range(num_band):
        dst_ds.GetRasterBand(i+1).ComputeStatistics(False)

    # Build overview layers for faster display.
    dst_ds.BuildOverviews('average', [2, 4, 8, 16, 32])

    # This will effectively close the file and flush data to disk.
    del dst_ds


if __name__ == '__main__':
    three_bands_stacked(
        band1_fn='./data/Landsat/Washington/p047r027_7t20000730_z10_nn30.tif',  # Red
        band2_fn='./data/Landsat/Washington/p047r027_7t20000730_z10_nn20.tif',  # Green
        band3_fn='./data/Landsat/Washington/p047r027_7t20000730_z10_nn10.tif',  # Blue
        num_band=3,
        out_path='./data/Landsat/Washington/nat_color.tif'
    )
