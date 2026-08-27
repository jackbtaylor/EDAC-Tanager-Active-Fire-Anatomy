# Anatomy of an Active Fire 🔥
A submission for Planet's *Tanager Open Data Competition* from the University of New Mexico's [Earth Data Analysis Center](https://edac.unm.edu/) Team 

**Lisa Sinclair** (Remote Sensing Scientist, Project Lead)<br> **Jack Taylor** (Programmer-Analyst, Author - burn severity analysis)<br> **Chris Girlamo** (Systems Administrator, Author - plume detection)

Contact: [Lisa Sinclair](mailto:llsinclair@unm.edu)

## Introduction
Active fire scenes present a complex spectral target for remote sensing. A single acquisition captures a smoke- and gas-laden atmosphere, fire-emitted radiance, and a highly heterogeneous surface. At fine spatial scales, flame, char, exposed soil, scorched vegetation, and unburned patches can occur within meters of one another. Conventional multispectral sensors represent that complexity with a limited number of broad spectral bands, and reduce it, in practice, to a single burn severity index. However, Tanager-1's 426 contiguous narrow bands spanning the visible through 2500 nm SWIR carry enough spectral detail to pull those unique atmospheric and surface components apart at fine spectral and spatial scales, all from a single acquisition.

The two notebooks in this repository leverage that spectral detail in sequence. The first identifies what the atmosphere is doing over the fire; the second sets those pixels aside and asks what the fire has done to the ground.

Notebook 1 (`1_toa_plume_detection.ipynb`) works in the *radiance* domain, before atmospheric correction. Atmospheric correction strips out the atmosphere's influence to isolate reflectance on the ground. The correction algorithm works by assuming the air is roughly uniform across the scene. A dense combustion plume is the opposite of uniform, so the plume is better studied in raw radiance than in a corrected product. Following the background modeling and suppression approach of [Ning et al. (2025)](https://doi.org/10.3390/agriculture15171835), the notebook synthesizes several mean background bands from the image, differences the full 426-band cube against that background, and normalizes the result band by band. It then averages across the bands near 2.0 µm, where smoke and combustion gases, particularly CO2, leave a distinct imprint on at-sensor radiance. Finally, it thresholds the result to separate plume-affected pixels from clear ones. Resolving that part of the spectrum is what makes the step possible: Tanager places 60 bands between 2000 and 2300 nm, where a broadband sensor places one. The notebook exports the extracted plume as both a GeoTIFF for use in a GIS and a NumPy array that feeds directly into Notebook 2.

Notebook 2 (`2_sr_burn_severity_analysis.ipynb`) works with Tanager's *surface reflectance* product. First, a combined mask removes three sets of pixels from the Area of Interest (AOI). Tanager's out of the box cloud and cirrus masks remove dense clouds, a reflectance threshold above 1.0 removes active flame, where a pixel is emitting energy rather than reflecting it. Finally, the plume mask from Notebook 1 removes the smoke plume and residual haze and leaves only clear ground. The final AOI is clipped to the July 24th fire perimeter digitized from the incident's Public Information Map published that day. The pixels within the AOI are spectrally unmixed using Multiple Endmember Spectral Mixture Analysis (MESMA), introduced by [Roberts et al. (1998)](https://doi.org/10.1016/S0034-4257(98)00037-6) with scripting tools developed by [Crabbé et al. (2020)](https://mesma.readthedocs.io/en/latest/). Reference spectra used for unmixing come from the Joint Fire Science Program's [spectral library](https://www.frames.gov/assessing-burn-severity/spectral-library/overview). The outcome is a set of continuous fractional cover maps, with individual products mapping the distribution of ash/char, green vegetation, non-photosynthetic vegetation, and soil/rock within the fire perimeter. These products are more physically interpretable than a unitless burn index, especially for non-remote sensing expert stakeholders who want to understand what is happening on the landscape. A land manager can estimate consumption from the ash fraction, locate surviving vegetation in the green vegetation fraction, flag the exposed-soil and charred areas that drive post-fire debris flow, and see where fuel remains available to burn.

Taken together, these notebooks demonstrate how valuable a single Tanager acquisition can be for both atmospheric and surface characterization of an event, even while an incident is active. 

## Image Background & Study Area:
> Hotshot crews are securing edges of unburned ground in Monroe Meadows and Monroe Canyon, eliminating isolated fuel pockets to reinforce contingency lines. Along Road 4088 to Hunts Lake, crews are using mastication and dozer work to support the development of strategic contingency lines, while water tenders and graders are in action maintaining access roads near Koosharem Canyon Road.
>
> Around the areas west of Magleby Reservoir, firefighters continue firing operations to create contingency opportunities to help prevent fire spread eastward. Hotshot crews in Manning Meadows are conducting additional direct firing operations to slow potential fire movement amid shifting winds and evolving weather patterns. These firing activities will likely cause increased smoke. 
>
> -- July 24th Morning Update, Monroe Canyon Fire

The Monroe Canyon Fire began on July 13th, 2025 from unknown causes. It ignited approximately 3 miles east of Monroe, in south-central Utah, USA. By the time it was declared contained on September 5th, it had burned over 73 thousand acres. When Tanager-1 flew over the fire on July 24th at 1pm, the fire was mapped at approximately 10 thousand acres and had already exhibited rapid growth and extreme fire behavior. 

The Tanager collection used for this workflow is available for free as a sample in the [Tanager Core Imagery Catalog](https://www.planet.com/data/stac/browser/tanager-core-imagery/fire/20250724_190927_83_4001/20250724_190927_83_4001.json). 

## Case Study Analysis:
This repository contains the tools necessary to gain a better understanding of the active fire behavior and the affected landscape surrounding the Monroe Canyon Fire on July 24th, 2025.

### Contents:
* 1_toa_plume_detection.ipynb
  * Use this notebook to analyze Tanager's *radiance* data, detecting a plume based on the spectral signature of CO2.
* 2_sr_burn_severity_analysis.ipynb
  * Use this notebook to analyze Tanager's *surface reflectance* data, using spectral unmixing to calculate portions of ash, vegetation, and soil per pixel.
 
### Installation:

1. Download the repository using `git clone https://github.com/jackbtaylor/EDAC-Tanager-Active-Fire-Anatomy`
2. Create a .venv in your preferred `python 3.12.x` development environment.
3. Install the required dependencies with `pip install -r requirements.txt`
4. Work through Notebook 1 to create required data for Notebook 2.
5. Explore results!
