"""
This module contains accessory functions used in the TOA smoke mask and SR burn severity analyses.
"""
import h5py
import rasterio
import re
import tqdm
import numpy as np
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view
from affine import Affine


def build_resampled_library(library_root, class_dirs, target_wavelengths, target_fwhm):
    """
    Walk each class subdirectory under library_root, load every .txt spectrum and
    resample each onto the tanager grid and return:

      library        : (n_bands, n_spectra) array for MesmaCore
      class_list     : (n_spectra,) array of class-name strings
      spectrum_names : (n_spectra,) array of source file stems

    Returns
    ----
    spectral library ready for MESMA, group class names, and source spectrum names
    """

    def _gaussian_response_resample(ref_wavelengths, ref_reflectance, target_wavelengths, target_fwhm, oversample_step=1.0):
        """
        Resample a reference spectrum onto tanager band centers by convolving
        with each band's Gaussian spectral response function (defined by
        center wavelength + FWHM).
        """
        ref_wavelengths = np.asarray(ref_wavelengths, dtype=float)
        ref_reflectance = np.asarray(ref_reflectance, dtype=float)

        # Remove wavelength gaps from library spectrum
        valid = ~np.isnan(ref_wavelengths) & ~np.isnan(ref_reflectance)
        ref_wavelengths = ref_wavelengths[valid]
        ref_reflectance = ref_reflectance[valid]

        # Sort the wavelengths
        order = np.argsort(ref_wavelengths)
        ref_wavelengths = ref_wavelengths[order]
        ref_reflectance = ref_reflectance[order]

        # Create a consistent sampling grid for weighting
        fine_grid = np.arange(ref_wavelengths.min(), ref_wavelengths.max(), oversample_step)
        fine_reflectance = np.interp(fine_grid, ref_wavelengths, ref_reflectance)

        resampled = np.full(len(target_wavelengths), np.nan)
        for i, (center, fwhm) in enumerate(zip(target_wavelengths, target_fwhm)):
            # Skip bands that fall outside the reference spectrum
            if not (ref_wavelengths.min() <= center <= ref_wavelengths.max()):
                continue
            # Convert fwhm from tanager metadata into standard deviation for the gaussian resample
            sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
            # Calculate spectral response using Gaussian weight of the fine grid point relative to this band center
            weights = np.exp(-0.5 * ((fine_grid - center) / sigma) ** 2)
            weight_sum = weights.sum()
            if weight_sum > 0:
                # Return the reference spectra aligned with tanager's band
                resampled[i] = np.sum(weights * fine_reflectance) / weight_sum

        return resampled

    
    library_spectra, class_list, spectrum_names = [], [], []

    for class_name in class_dirs:
        class_path = library_root / class_name
        n_loaded = 0
        for txt_file in sorted(class_path.glob("*.txt")):
            df = pd.read_csv(txt_file, sep="\t")
            df.columns = df.columns.str.lower()
            if "mean" not in df.columns:
                print(f"  SKIP {txt_file.name}: columns = {list(df.columns)}")
                continue

            resampled = _gaussian_response_resample(df["wavelength"].values, df["mean"].values, target_wavelengths, target_fwhm)
            library_spectra.append(resampled)
            class_list.append(class_name)
            spectrum_names.append(txt_file.stem)
            n_loaded += 1

        print(f"{class_name}: {n_loaded} spectra loaded")

    library = np.array(library_spectra).T   # (n_bands, n_spectra)
    return library, np.array(class_list), np.array(spectrum_names)


def generate_new_background(back, window_size=(3, 3), n_bands=None): # function to generate new background
    H, W = back.shape[:2] # pull out the height and width from the shape
    if n_bands is None: # calculate the number of bands 
        n_bands = back.shape[2]

    wh, ww = window_size # assign window variables 
    pad_h, pad_w = wh // 2, ww // 2 # pad variables  by half their value rounded down
    rng = np.random.default_rng() # random number generator 

    new_bands = [] # make empty list 
    for b in tqdm.tqdm(range(n_bands)): # cycle through the bands, tqdm gives a progress bar
        band = back[:, :, b] # select individual band

        # pad so the sliding-window output matches the original H, W
        padded = np.pad(band, ((pad_h, pad_h), (pad_w, pad_w)), mode='reflect') 
        window = sliding_window_view(padded, window_shape=window_size)  # (H, W, wh, ww) make the sliding window
        out_h, out_w = window.shape[:2] # create the out height and width

        flat = window.reshape(out_h, out_w, wh * ww) # make the window flat in the last dimesion (i.e. band)
        idx = rng.integers(0, wh * ww, size=(out_h, out_w)) # make random values in the range of the flattened band
        moving_sample = np.take_along_axis(flat, idx[..., None], axis=-1)[..., 0]  # randomly select a value from the background values in the flattened array

        new_bands.append(moving_sample) # build array 

    return np.stack(new_bands, axis=-1) # return the new array 


def georeference_h5(img, northern=True):
    """
    Build georeferencing info (transform, CRS, extent) for a Tanager HDFEOS-gridded HDF5 product.

    img         : an open hdf file using h5py and with
    northern    : is the image in the northern hemisphere (for projection codes)

    Returns
    -------
    dict with: transform, crs, shape, grid_name, xmin, xmax, ymin, ymax
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


def print_structure(name, obj):
    """Callback function to print the name and type of each object."""
    indent = "  " * name.count('/')
    obj_type = "Group" if isinstance(obj, h5py.Group) else "Dataset"
    print(f"{indent}{name} ({obj_type})")