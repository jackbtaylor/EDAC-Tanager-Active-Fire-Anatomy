"""
This module contains accessory functions used in the TOA smoke mask and SR burn severity analyses.
"""
import h5py
import rasterio
import re
import numpy as np
from affine import Affine


def build_resampled_library(library_root, class_dirs, target_wavelengths, target_fwhm):
    """
    Walk each class subdirectory under library_root, load every .txt spectrum and
    resample each onto the tanager grid and return:

      library        : (n_bands, n_spectra) array for MesmaCore
      class_list     : (n_spectra,) array of class-name strings
      spectrum_names : (n_spectra,) array of source file stems
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



def georeference_h5(img, northern=True):
    """
    Build georeferencing info (transform, CRS, extent) for a Tanager HDFEOS-gridded HDF5 product. 

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
