# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT

import itertools
def parse_year_wildcard(w):
    """
    Parse a {year} wildcard to a list of years. Parses a wildcard with years and ranges separated by '+' and '-'.
    
    Args:
    w (str): A string containing years and ranges like '1980+1990+2000-2002'.
    
    Returns:
    list: A list of years expanded from the ranges.
    """
    years = []
    parts = w.split('+')  # Split the string by '+'
    
    for part in parts:
        if '-' in part:
            # It's a range, split it into start and end, then generate the range of years.
            start, end = map(int, part.split('-'))
            # Ensure that the start year is less than the end year.
            if end < start:
                raise ValueError(f"Malformed range of years {part}.")
            years.extend(range(start, end + 1))  # Include the end year in the range.
        else:
            # It's a single year, just append it to the list.
            years.append(int(part))
    
    return sorted(years)

localrules:
    all,
    cluster_networks,
    extra_components_networks,
    prepare_elec_networks,
    prepare_sector_networks,
    solve_elec_networks,
    solve_sector_networks,
    plot_networks,


rule all:
    input:
        RESULTS + "graphs/costs.pdf",
    default_target: True


rule cluster_networks:
    input:
        expand(
            RESOURCES + "networks/elec{weather_year}_s{simpl}_{clusters}.nc",
            **config["scenario"]
        ),


rule extra_components_networks:
    input:
        expand(
            RESOURCES + "networks/elec{weather_year}_s{simpl}_{clusters}_ec.nc",
            **config["scenario"]
        ),


rule prepare_elec_networks:
    input:
        expand(
            RESOURCES
            + "networks/elec{weather_year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule prepare_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),

rule prepare_brownfield_sector_networks:
    input:
        expand(
            RESULTS
            + "prenetworks-brownfield/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            weather_year=parse_year_wildcard(config["scenario"]["weather_year"]),
            simpl=config["scenario"]["simpl"],
            clusters=config["scenario"]["clusters"],
            ll=config["scenario"]["ll"],
            opts=config["scenario"]["opts"],
            sector_opts=config["scenario"]["sector_opts"],
            planning_horizons=config["scenario"]["planning_horizons"]
        ),

rule solve_elec_networks:
    input:
        expand(
            RESULTS
            + "networks/elec{weather_year}_s{simpl}_{clusters}_ec_l{ll}_{opts}.nc",
            **config["scenario"]
        ),


rule solve_sector_networks:
    input:
        expand(
            RESULTS
            + "postnetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),


rule plot_networks:
    input:
        expand(
            RESULTS
            + "maps/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),
