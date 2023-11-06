# SPDX-FileCopyrightText: : 2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT


localrules:
    copy_config,
    copy_conda_env,


rule plot_network:
    params:
        foresight=config["foresight"],
        plotting=config["plotting"],
    input:
        network=RESULTS
        + "postnetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
        regions=RESOURCES
        + "regions_onshore_elec{weather_year}_s{simpl}_{clusters}.geojson",
    output:
        map=RESULTS
        + "maps/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
        today=RESULTS
        + "maps/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}-today.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    benchmark:
        (
            BENCHMARKS
            + "plot_network/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}"
        )
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_network.py"


rule copy_config:
    params:
        RDIR=RDIR,
    output:
        RESULTS + "config.yaml",
    threads: 1
    resources:
        mem_mb=1000,
    benchmark:
        BENCHMARKS + "copy_config"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/copy_config.py"


rule copy_conda_env:
    output:
        RESULTS + "config/environment.yaml",
    threads: 1
    resources:
        mem_mb=500,
    log:
        LOGS + "copy_conda_env.log",
    benchmark:
        BENCHMARKS + "copy_conda_env"
    conda:
        "../envs/environment.yaml"
    shell:
        "conda env export -f {output} --no-builds"


rule make_summary:
    params:
        foresight=config["foresight"],
        costs=config["costs"],
        snapshots=config["snapshots"],
        scenario=config["scenario"],
        RDIR=RDIR,
    input:
        networks=expand(
            RESULTS
            + "postnetworks/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}_{planning_horizons}.nc",
            **config["scenario"]
        ),
        costs="data/costs_{}.csv".format(config["costs"]["year"])
        if config["foresight"] == "overnight"
        else "data/costs_{}.csv".format(config["scenario"]["planning_horizons"][0]),
        plots=expand(
            RESULTS
            + "maps/elec{weather_year}_s{simpl}_{clusters}_l{ll}_{opts}_{sector_opts}-costs-all_{planning_horizons}.pdf",
            **config["scenario"]
        ),
    output:
        nodal_costs=RESULTS + "csvs/nodal_costs.csv",
        nodal_capacities=RESULTS + "csvs/nodal_capacities.csv",
        nodal_cfs=RESULTS + "csvs/nodal_cfs.csv",
        cfs=RESULTS + "csvs/cfs.csv",
        costs=RESULTS + "csvs/costs.csv",
        capacities=RESULTS + "csvs/capacities.csv",
        curtailment=RESULTS + "csvs/curtailment.csv",
        energy=RESULTS + "csvs/energy.csv",
        supply=RESULTS + "csvs/supply.csv",
        supply_energy=RESULTS + "csvs/supply_energy.csv",
        prices=RESULTS + "csvs/prices.csv",
        weighted_prices=RESULTS + "csvs/weighted_prices.csv",
        market_values=RESULTS + "csvs/market_values.csv",
        price_statistics=RESULTS + "csvs/price_statistics.csv",
        metrics=RESULTS + "csvs/metrics.csv",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        LOGS + "make_summary.log",
    benchmark:
        BENCHMARKS + "make_summary"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/make_summary.py"


rule plot_summary:
    params:
        countries=config["countries"],
        planning_horizons=config["scenario"]["planning_horizons"],
        sector_opts=config["scenario"]["sector_opts"],
        emissions_scope=config["energy"]["emissions"],
        eurostat_report_year=config["energy"]["eurostat_report_year"],
        plotting=config["plotting"],
        RDIR=RDIR,
    input:
        costs=RESULTS + "csvs/costs.csv",
        energy=RESULTS + "csvs/energy.csv",
        balances=RESULTS + "csvs/supply_energy.csv",
        eurostat=input_eurostat,
        co2="data/bundle-sector/eea/UNFCCC_v23.csv",
    output:
        costs=RESULTS + "graphs/costs.pdf",
        energy=RESULTS + "graphs/energy.pdf",
        balances=RESULTS + "graphs/balances-energy.pdf",
    threads: 2
    resources:
        mem_mb=10000,
    log:
        LOGS + "plot_summary.log",
    benchmark:
        BENCHMARKS + "plot_summary"
    conda:
        "../envs/environment.yaml"
    script:
        "../scripts/plot_summary.py"
