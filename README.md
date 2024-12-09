# VeCell

* **Developed for:** Naomie
* **Team:** Cohen-Salmon
* **Date:** October 2024
* **Software:** Fiji


### Images description

3D images taken on an AxioZoom microscope.

2 channels:
  1. *DsRed:* Astrocytes (mandatory)
  2. *Cy5:* Vessels (optional)
  
With each image, a *.roi* or *.zip* file containing one or multiple ROIs should be provided; otherwise, the image is not analyzed.

### Plugin description

In each ROI:
* Detect astrocytes with Cellpose
* Compute distance between each cell and its nearest neighbors
* Compute G-function related spatial distribution index of the population of cells
* If vessels channel provided:
  * Detect vessels with median filtering + DoG filtering + thresholding + closing filtering + median filtering
  * Compute vessels skeleton and provide vessels diameter, length, branches number, and junctions number
  * Compute distance between each cell and its nearest vessel


### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin
* **Cellpose** conda environment + *cyto* pretrained model or another fine-tuned model

### Version history

Version 5 released on October 23, 2024.

Improvements compared to version 4:
* Plugin renamed
* Code cleaned
* Dialog box changed
* Vessels segmentation improved:
  * Channel quantile-based normalization
  * 2 DoG filters can be applied if thin and thick vessels appear in the image
* Vessels skeleton branches filtered out by length
* Vessels skeleton length, branches number, and junctions number provided
