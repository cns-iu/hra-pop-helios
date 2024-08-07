# hra-pop-helios

Creates a k-nearest neighbors network for biological data and prepare it to be exported to Helios-Web (heliosweb.io).



# Documentation

We first apply a typical single-cell analysis preprocessing on the data. In particular, first we find and filter data by  the highly variable genes (using Seurat approach), apply log1p transformation, and then apply PCA or TruncatedSVD to the data.

We are currently using 10000 highly variable genes and 100 PCA dimensions. [These values can be tuned for better cell type separation or better visualizations]

For highly variable genes, we are using Seurat approach (implemented in ScanPy):
https://doi.org/10.1016/j.cell.2019.05.031

Next, we use approximate nearest neighbors according to cosine similarity to obtain the edges of the network. The method is described here:
https://doi.org/10.1145/1963405.1963487
We are using the implementation from the PyNNDescent package.

We use Helios-Web (http://heliosweb.io) to visualize the obtained network. Currently we have a javascript implementation of a force-directed layout algorithm, which uses a single thread. The current implementation is based on the forceAtlas2 algorithm (https://doi.org/10.1371/journal.pone.0098679).


An alternative implementation of that in C is available and can be used for larger networks (The layout won't be interactive).

Code is available from:
https://github.com/filipinascimento/bioknn


# Usage

First download the data. Currently we are using the GTEx 8 tissues snRNA-Seq Data, which can be downloaded from:
https://gtexportal.org/home/downloads/adult-gtex/single_cell (`GTEx_8_tissues_snRNAseq_atlas_071421.public_obs.h5ad`)

Create a folder named `Data`  Move the data file to that folder. Data should be in `h5ad` format.
Packages from `requirements.txt` should be installed with `pip install -r requirements.txt`

The current version uses a single script: `generateKNN.py`
After setting the parameters inside the script, run it with `python Scripts/generateKNN.py`
Parameters can be changed inside the script.

The script will generate a .xnet file in the `Networks` folder.


To run helios-web for visualization:
Preferentially using Chrome or Safari. Firefox should work but it remains untested.
Open heliosweb.io and choose a configuration (light, dark, density, etc)
Add `?use2d` as argument to the heliosweb url to use the 2d mode. E.g. `http://heliosweb.io/docs/example/?advanced&use2d`

Drag and drop the generated .xnet file to the window.
It should be replaced by the new network.

If the layout didn't start by itself, press space to run the algorithm. You can press space again to stop or resume.

# TODO
 - Convert into a package
 - Will update Helios-Web with beta version which supports ForceAtlas2 without the degree normalization
 - Include helios-web server (ailens), which should launch the browser with the network visualization
 - Use the Helios-Core (https://github.com/filipinascimento/helios) library for layouting pre visualization

