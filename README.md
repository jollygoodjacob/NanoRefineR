# NanoRefineR: A Shiny app for postprocessing and visualizing noisy Nano-Hyperspec HSI data cubes

NanoRefineR is a small project I developed to enable users to apply certain postprocessing techniques to noisy Headwall Nano-Hyperspec HSI cubes and visualize the spectral differences in real time. This tool allows users to rapidly determine which post-processing corrections they may want to make to HSI cubes and can also be used to post-process data in batch.

NanoRefineR provides two main functionalities:
- a post-processing component that allows for spatial downsampling, spectral binning, and Savitzky-Golay smoothing of input datasets
- a visualizer with control over band combinations and the ability to click on areas of an image and see the spectra. 

<img src="https://github.com/jollygoodjacob/NanoRefineR/raw/master/imgs/img.png" style="width:100%; border:1px solid #ccc;"/>


## Processing
The processing panel allows users to batch process multiple hyperspectral cubes from a specified folder using the following selectable operations:

1. **Spatial Downscaling**
Reprojects and resamples the input cube to a user-defined coordinate reference system (CRS) and spatial resolution.


2. **Spectral Binning**
Aggregates spectral bands into bins of user-defined width (in nanometers). This is done by averaging reflectance values across bands falling into each bin. Binned bands are renamed using the bin center wavelength to preserve spectral traceability.

3. **Savitzky-Golay Smoothing**
Applies a Savitzky-Golay filter along the spectral dimension to reduce noise while preserving the shape of spectral features. Users can define the polynomial order and filter length (must be odd).

All processing steps are optional and can be enabled/disabled via checkboxes. Output files are saved in GeoTIFF format with a custom suffix, allowing easy distinction from the original input files.

Batch processing is parallelized using `doParallel` for faster execution across multiple CPU cores.

## Visualization
The visualization panel provides real-time inspection of processed results:

**File Selector**: Choose any processed output to visualize.

**RGB Plot Controls**: Dynamically select red, green, and blue wavelengths for generating RGB composites. The closest matching band is automatically selected.

**Interactive RGB Display**: View side-by-side RGB composites of the input and processed files.

**Spectral Explorer**: Click on any pixel in either image to view its spectral profile.

Both input and output spectra for the selected location are shown in parallel.

This layout is designed for intuitive quality control and to help identify the most effective postprocessing settings for noisy hyperspectral data.











