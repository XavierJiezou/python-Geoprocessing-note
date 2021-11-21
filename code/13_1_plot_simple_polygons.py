# Plot world landmasses as simple polygons.
import matplotlib.pyplot as plt
from osgeo import ogr


def plot_simple_polygons(fn):
    ds = ogr.Open(fn)
    lyr = ds.GetLayer(0)
    for row in lyr:
        geom = row.geometry()
        ring = geom.GetGeometryRef(0)
        coords = ring.GetPoints()
        x, y = zip(*coords)
        plt.plot(x, y, 'k')  # 'k' means black
    # Equalize the axis units so things aren't warped. Comment out
    # this line and see what happens.
    plt.axis('equal')
    plt.show()


if __name__ == '__main__':
    plot_simple_polygons('data/global/ne_110m_land.shp')
