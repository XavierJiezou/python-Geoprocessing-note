# Function to retrieve overview data from a raster.

from osgeo import gdal


def get_overview_data(fn: str, band_index: int = 1, level: int = -1):
    """Returns an array containing data from an overview.

    Args:
        fn (str): path to raster file.
        band_index (int, optional): band number to get overview for.
        level (int, optional): overview level, where 1 is the highest resolution; the coarsest can be retrieved with -1.
    """
    ds = gdal.Open(fn)
    band = ds.GetRasterBand(band_index)
    if level > 0:
        ov_band = band.GetOverview(level)
    else:
        num_ov = band.GetOverviewCount()
        ov_band = band.GetOverview(num_ov + level)
    return ov_band.ReadAsArray()


if __name__ == '__main__':
    get_overview_data('data/Landsat/Washington/nat_color.tif', 1, -1)
