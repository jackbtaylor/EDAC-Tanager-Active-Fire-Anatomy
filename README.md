# Background:
> Hotshot crews are securing edges of unburned ground in Monroe Meadows and Monroe Canyon, eliminating isolated fuel pockets to reinforce contingency lines. Along Road 4088 to Hunts Lake, crews are using mastication and dozer work to support the development of strategic contingency lines, while water tenders and graders are in action maintaining access roads near Koosharem Canyon Road.
>
> Around the areas west of Magleby Reservoir, firefighters continue firing operations to create contingency opportunities to help prevent fire spread eastward. Hotshot crews in Manning Meadows are conducting additional direct firing operations to slow potential fire movement amid shifting winds and evolving weather patterns. These firing activities will likely cause increased smoke. 
>
> -- July 24th Morning Update, Monroe Canyon Fire

The Monroe Canyon Fire began on July 13th, 2025 from unknown causes. It ignited approximately 3 miles east of Monroe, in south-central Utah, USA. By the time it was declared contained on September 5th, it had burned over 73 thousand acres. When Tanager-1 flew over the fire on July 24th at 1pm, the fire was mapped at approximately 10 thousand acres and had already exhibited rapid growth and extreme fire behavior. 

The Tanager collection used for this workflow is available for free in the [Tanager Core Imagery Catalog](https://www.planet.com/data/stac/browser/tanager-core-imagery/fire/20250724_190927_83_4001/20250724_190927_83_4001.json). 

# Analysis:
This repository contains the tools necessary to gain a better understanding of the active fire behavior and the affected landscape.

**Notebook Contents:**
* 1_toa_plume_detection.ipynb
  * Use this notebook to analyze Tanager's *radiance* data, detecting a plume based on the spectral signature of CO2.
* 2_sr_burn_severity_analysis.ipynb
  * Use this notebook to analyze Tanager's *surface reflectance* data, using spectral unmixing to calculate portions of ash, vegetation, and soil per pixel.
 
**Installation:**

1. Download the repository using `git clone https://github.com/jackbtaylor/EDAC-Tanager-Active-Fire-Anatomy`
2. Create a .venv in your preferred python development environment.
3. Install the required dependencies with `pip install -r requirements.txt`
4. Work through Notebook 1 to create required data for Notebook 2.
5. Explore results!
