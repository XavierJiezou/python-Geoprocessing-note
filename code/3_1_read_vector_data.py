# Script to read attribute from a shapefile
from osgeo import ogr  # Don't forget to import ogr


def read_vector_data(fn: str, num: int = 10):
    """Printing data from the first ten features in a shapefile.

    Args:
        fn (str): File name of the shapefile.
        num (int, optional): Number to printing features. Defaults to 10.
    """
    # Open the data source and get the layer
    ds = ogr.Open(fn, 0) # 0: read only 1: read & write
    if ds is None:
        raise FileNotFoundError(f"[Errno 2] No such file or directory: '{fn}'")
    else:
        ly = ds.GetLayer(0)

    i = 0
    for ft in ly:
        # Get the x,y coordinates
        pt = ft.geometry()
        x = pt.GetX()
        y = pt.GetY()

        # Get the attribute values
        name = ft.GetField('Name')
        pop = ft.GetField('POP_MAX')
        print(name, pop, x, y)
        i += 1
        if i == num:
            break
        else:
            continue
    del ds


if __name__ == '__main__':
    read_vector_data('data/global/ne_50m_populated_places.shp')


"""Expected Output
Bombo 75000 32.533299524864844 0.5832991056146284
Fort Portal 42670 30.27500161597942 0.671004121125236       
Potenza 69060 15.798996495640267 40.642002130098206
Campobasso 50762 14.655996558921856 41.56299911864397       
Aosta 34062 7.315002595706176 45.7370010670723
Mariehamn 10682 19.949004471869102 60.09699618489543        
Ramallah 24599 35.206209378189556 31.90294475142406
Vatican City 832 12.453386544971766 41.903282179960115      
Poitier 85960 0.3332765285345545 46.58329225573658
Clermont-Ferrand 233050 3.080008095928406 45.779982115759424
"""
