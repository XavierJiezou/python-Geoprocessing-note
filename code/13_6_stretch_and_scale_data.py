# Function to stretch and scale data.

import matplotlib.pyplot as plt
from osgeo import gdal
import numpy as np
import os


def stretch_data(data: np.array, num_stddev: int = 2):
    """Returns the data with a standard deviation stretch applied.

    Args:
        data (np.array): array containing data to stretch.
        num_stddev (int, optional): number of standard deviations to use. Defaults to 2.

    Returns:
        np.array: stretched data.
    """
    mean = np.mean(data)
    std_range = np.std(data) * 2
    new_min = max(mean - std_range, np.min(data))
    new_max = min(mean + std_range, np.max(data))
    clipped_data = np.clip(data, new_min, new_max)
    return clipped_data / (new_max - new_min)


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

if __name__=='__main__':
    red_fn = 'data/Landsat/Washington/p047r027_7t20000730_z10_nn30.tif'
    green_fn = 'data/Landsat/Washington/p047r027_7t20000730_z10_nn20.tif'
    blue_fn = 'data/Landsat/Washington/p047r027_7t20000730_z10_nn10.tif'

    # Try plotting 3 bands.
    red_data = get_overview_data(red_fn)
    green_data = get_overview_data(green_fn)
    blue_data = get_overview_data(blue_fn)
    data = np.dstack((red_data, green_data, blue_data))
    plt.imshow(data)
    plt.axis('off')
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
    plt.savefig('img/13.13.a.png')
    plt.show()


    # Plot 3 stretched bands.
    red_data = stretch_data(get_overview_data(red_fn), 2)
    green_data = stretch_data(get_overview_data(green_fn), 2)
    blue_data = stretch_data(get_overview_data(blue_fn), 2)
    alpha = np.where(red_data + green_data + blue_data > 0, 1, 0)
    data = np.dstack((red_data, green_data, blue_data, alpha))
    plt.imshow(data)
    plt.axis('off')
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
    plt.savefig('img/13.13.b.png')
    plt.show()
