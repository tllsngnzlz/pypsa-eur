#!/usr/bin/env python
# -*- coding: utf-8 -*-

# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""
Build hydroelectric inflow time-series for each country.

Relevant Settings
-----------------

.. code:: yaml

    countries:

    renewable:
        hydro:
            cutout:
            clip_min_inflow:

.. seealso::
    Documentation of the configuration file ``config/config.yaml`` at
    :ref:`toplevel_cf`, :ref:`renewable_cf`

Inputs
------

- ``data/bundle/EIA_hydro_generation_2000_2014.csv``: Hydroelectricity net generation per country and year (`EIA <https://www.eia.gov/beta/international/data/browser/#/?pa=000000000000000000000000000000g&c=1028i008006gg6168g80a4k000e0ag00gg0004g800ho00g8&ct=0&ug=8&tl_id=2-A&vs=INTL.33-12-ALB-BKWH.A&cy=2014&vo=0&v=H&start=2000&end=2016>`_)

    .. image:: img/hydrogeneration.png
        :scale: 33 %

- ``resources/country_shapes.geojson``: confer :ref:`shapes`
- ``"cutouts/" + config["renewable"]['hydro']['cutout']``: confer :ref:`cutout`

Outputs
-------

- ``resources/profile_hydro.nc``:

    ===================  ================  =========================================================
    Field                Dimensions        Description
    ===================  ================  =========================================================
    inflow               countries, time   Inflow to the state of charge (in MW),
                                           e.g. due to river inflow in hydro reservoir.
    ===================  ================  =========================================================

    .. image:: img/inflow-ts.png
        :scale: 33 %

    .. image:: img/inflow-box.png
        :scale: 33 %

Description
-----------

.. seealso::
    :mod:`build_renewable_profiles`
"""

import logging
import xarray as xr
import atlite
import country_converter as coco
import geopandas as gpd
import pandas as pd
from _helpers import configure_logging

cc = coco.CountryConverter()


def get_eia_annual_hydro_generation(fn, countries):
    # in billion kWh/a = TWh/a
    df = pd.read_csv(fn, skiprows=2, index_col=1, na_values=[" ", "--"]).iloc[1:, 1:]
    df.index = df.index.str.strip()

    former_countries = {
        "Former Czechoslovakia": dict(
            countries=["Czech Republic", "Slovakia"], start=1980, end=1992
        ),
        "Former Serbia and Montenegro": dict(
            countries=["Serbia", "Montenegro"], start=1992, end=2005
        ),
        "Former Yugoslavia": dict(
            countries=[
                "Slovenia",
                "Croatia",
                "Bosnia and Herzegovina",
                "Serbia",
                "Montenegro",
                "North Macedonia",
            ],
            start=1980,
            end=1991,
        ),
    }

    for k, v in former_countries.items():
        period = [str(i) for i in range(v["start"], v["end"] + 1)]
        ratio = df.loc[v["countries"]].T.dropna().sum()
        ratio /= ratio.sum()
        for country in v["countries"]:
            df.loc[country, period] = df.loc[k, period] * ratio[country]

    baltic_states = ["Latvia", "Estonia", "Lithuania"]
    df.loc[baltic_states] = (
        df.loc[baltic_states].T.fillna(df.loc[baltic_states].mean(axis=1)).T
    )

    df.loc["Germany"] = df.filter(like="Germany", axis=0).sum()
    df.loc["Serbia"] += df.loc["Kosovo"].fillna(0.0)
    df = df.loc[~df.index.str.contains("Former")]
    df.drop(["Europe", "Germany, West", "Germany, East", "Kosovo"], inplace=True)

    df.index = cc.convert(df.index, to="iso2")
    df.index.name = "countries"

    df = df.T[countries] * 1e6  # in MWh/a

    return df


logger = logging.getLogger(__name__)

def replace_xarray_index_year(data_array):
   
    # New time index
    year = snakemake.config['snapshots'].get('year', '2013')
    boundary = snakemake.config['snapshots'].get('year_boundary', '01-01')
    new_index = pd.date_range(
            f"{year}-{boundary}",
            end=f"{int(year) + 1}-{boundary}",
            freq="h",
            inclusive="left",
        )

    # Handling the leap day
    new_has_leap = (new_index.month == 2) & (new_index.day == 29)
    data_has_leap = (data_array.time.dt.month == 2) & (data_array.time.dt.day == 29)

    # If the new year has a Feb 29, but data_array doesn't
    if new_has_leap.any() and not data_has_leap.any():
        # Fill missing values with data from the day before the gap
        leap_day = data_array.sel(time=(data_array.time.dt.month == 2) & (data_array.time.dt.day == 28))
        
        # Generate a new time index just for February 29th of the new year      
        leap_day = leap_day.assign_coords(time=pd.date_range(start=f"{year}-02-29", periods=24, freq='H'))        
        # Concatenate data arrays around the leap day
        data_array = xr.concat([data_array.sel(time=data_array.time.dt.month < 3 ), leap_day, data_array.sel(time=data_array.time.dt.month >= 3)], dim="time")
        
    # If data_array has a Feb 29, but the new year doesn't
    elif data_has_leap.any() and not new_has_leap.any():
        # Drop February 29 data from data_array
        data_array = data_array.where(~data_has_leap, drop=True)
    # Assign the new time index to the data
    new_data_array = data_array.assign_coords(time=new_index).chunk({"time": 100})
    return new_data_array

if __name__ == "__main__":
    if "snakemake" not in globals():
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("build_hydro_profile")
    configure_logging(snakemake)

    params_hydro = snakemake.params.hydro
    
    countries = snakemake.params.countries
    country_shapes = (
        gpd.read_file(snakemake.input.country_shapes)
        .set_index("name")["geometry"]
        .reindex(countries)
    )
    country_shapes.index.name = "countries"

    fn = snakemake.input.eia_hydro_generation
    eia_stats = get_eia_annual_hydro_generation(fn, countries)

    weather_year = snakemake.wildcards.weather_year
    norm_year = params_hydro.get("eia_norm_year")
    if norm_year:
        eia_stats.loc[weather_year] = eia_stats.loc[norm_year]
    elif weather_year and weather_year not in eia_stats.index:
        eia_stats.loc[weather_year] = eia_stats.median()
    cutout = atlite.Cutout(snakemake.input.cutout)
    
    inflow = cutout.runoff(
        shapes=country_shapes,
        smooth=True,
        lower_threshold_quantile=True,
        normalize_using_yearly=eia_stats,
    )

    if "clip_min_inflow" in params_hydro:
        inflow = inflow.where(inflow > params_hydro["clip_min_inflow"], 0)

    inflow = replace_xarray_index_year(inflow)

    inflow.to_netcdf(snakemake.output.profile)
