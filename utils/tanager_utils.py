"""
This module contains accessory functions used in the TOA smoke mask and SR burn severity analyses.
"""
import h5py
import rasterio
import re
from affine import Affine


def georeference_h5(img, northern=True):
    """
    Build georeferencing info (transform, CRS, extent) for a Tanager HDFEOS-gridded HDF5 product. 

    Returns
    -------
    dict with: transform, crs, shape, grid_name, xmin, xmax, ymin, ymax
    (xmin/xmax/ymin/ymax are ready to drop into imshow's extent=[...]).
    """
    def _parse_struct_metadata(img):
        """
        Parse an HDFEOS StructMetadata.0 blob for grid geometry.
        """

        raw = img["HDFEOS INFORMATION/StructMetadata.0"][()]
        meta = raw.decode() if isinstance(raw, bytes) else raw

        grid_name = re.search(r'GridName="([^"]+)"', meta).group(1)
        ul_x, ul_y = re.search(r"UpperLeftPointMtrs=\(([-\d.]+),([-\d.]+)\)", meta).groups()
        lr_x, lr_y = re.search(r"LowerRightMtrs=\(([-\d.]+),([-\d.]+)\)", meta).groups()
        zone_code = int(re.search(r"ZoneCode=(\d+)", meta).group(1))
        xdim = int(re.search(r"XDim=(\d+)", meta).group(1))
        ydim = int(re.search(r"YDim=(\d+)", meta).group(1))

        return {
            "grid_name": grid_name,
            "ul_x": float(ul_x), "ul_y": float(ul_y),
            "lr_x": float(lr_x), "lr_y": float(lr_y),
            "zone_code": zone_code, "xdim": xdim, "ydim": ydim,
        }
    
    info = _parse_struct_metadata(img)

    pixel_x = (info["lr_x"] - info["ul_x"]) / info["xdim"]
    pixel_y = (info["lr_y"] - info["ul_y"]) / info["ydim"]
    transform = Affine.translation(info["ul_x"], info["ul_y"]) * Affine.scale(pixel_x, pixel_y)

    epsg = (32600 if northern else 32700) + info["zone_code"]
    crs = rasterio.crs.CRS.from_epsg(epsg)

    shape = (info["ydim"], info["xdim"])
    return {
        "transform": transform,
        "crs": crs,
        "shape": shape,
        "grid_name": info["grid_name"],
        "xmin": info["ul_x"], "xmax": info["lr_x"],
        "ymin": info["lr_y"], "ymax": info["ul_y"],
    }
