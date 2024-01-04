# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2022 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Transforms the global ship density data from the `World Bank Data Catalogue.

<https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density>`_
to the size of the considered cutout. The global ship density raster is later
used for the exclusion when calculating the offshore potentials.

Relevant Settings
-----------------

.. code:: yaml

    renewable:
        {technology}:
            cutout:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`renewable_cf`

Inputs
------

- ``data/bundle/shipdensity/shipdensity_global.zip``: Global shipping traffic
  density from `World Bank Data Catalogue
  <https://datacatalog.worldbank.org/search/dataset/0037580/>`_.

Outputs
-------

- ``resources/europe_shipdensity_raster.nc``: Reduced version of global shipping
  traffic density from `World Bank Data Catalogue
  <https://datacatalog.worldbank.org/search/dataset/0037580/>`_ to reduce
  computation time.

Description
-----------
"""

import logging
import os
import zipfile
import rioxarray
from _helpers import configure_logging, mock_snakemake
from build_natura_raster import determine_cutout_xXyY

logger = logging.getLogger(__name__)

def extract_and_process_ship_density(tif_path, snakemake):
    cutouts = snakemake.input.cutouts
    xs, Xs, ys, Ys = zip(*(determine_cutout_xXyY(cutout) for cutout in cutouts))

    with rioxarray.open_rasterio(tif_path) as ship_density:
        ship_density = ship_density.drop(["band"]).sel(
            x=slice(min(xs), max(Xs)), y=slice(max(Ys), min(ys))
        )
        ship_density.rio.to_raster(snakemake.output[0])

if __name__ == "__main__":
    if "snakemake" not in globals():
        snakemake = mock_snakemake("build_ship_raster", weather_year=1980)
    
    configure_logging(snakemake)
    
    tif_path = "shipdensity_global.tif"
    
    if not os.path.exists(tif_path):
        try:
            with zipfile.ZipFile(snakemake.input.ship_density) as zip_f:
                zip_f.extract(tif_path)
                extract_and_process_ship_density(tif_path, snakemake)
        except Exception as e:
            logger.error(f"Error processing ship density data: {e}")
    elif os.path.exists(tif_path):
        extract_and_process_ship_density(tif_path, snakemake)
