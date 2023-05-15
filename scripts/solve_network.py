"""Solve network."""

import logging
import re

import numpy as np
import pandas as pd
import pypsa
from helper import override_component_attrs
from pypsa.linopf import ilopf, network_lopf
from pypsa.linopt import define_constraints, get_var, join_exprs, linexpr
from vresutils.benchmark import memory_logger

logger = logging.getLogger(__name__)
pypsa.pf.logger.setLevel(logging.WARNING)


def add_land_use_constraint(n):
    if "m" in snakemake.wildcards.clusters:
        _add_land_use_constraint_m(n)
    else:
        _add_land_use_constraint(n)


def _add_land_use_constraint(n):
    # warning: this will miss existing offwind which is not classed
    # AC-DC and has carrier 'offwind'

    for carrier in ["solar", "onwind", "offwind-ac", "offwind-dc"]:
        existing = (
            n.generators.loc[n.generators.carrier == carrier, "p_nom"]
            .groupby(n.generators.bus.map(n.buses.location))
            .sum()
        )
        existing.index += " " + carrier + "-" + snakemake.wildcards.planning_horizons
        n.generators.loc[existing.index, "p_nom_max"] -= existing

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def _add_land_use_constraint_m(n):
    # if generators clustering is lower than network clustering,
    # land_use accounting is at generators clusters

    planning_horizons = snakemake.config["scenario"]["planning_horizons"]
    grouping_years = snakemake.config["existing_capacities"]["grouping_years"]
    current_horizon = snakemake.wildcards.planning_horizons

    for carrier in ["solar", "onwind", "offwind-ac", "offwind-dc"]:
        existing = n.generators.loc[n.generators.carrier == carrier, "p_nom"]
        ind = list(
            set(
                [
                    i.split(sep=" ")[0] + " " + i.split(sep=" ")[1]
                    for i in existing.index
                ]
            )
        )

        previous_years = [
            str(y)
            for y in planning_horizons + grouping_years
            if y < int(snakemake.wildcards.planning_horizons)
        ]

        for p_year in previous_years:
            ind2 = [
                i for i in ind if i + " " + carrier + "-" + p_year in existing.index
            ]
            sel_current = [i + " " + carrier + "-" + current_horizon for i in ind2]
            sel_p_year = [i + " " + carrier + "-" + p_year for i in ind2]
            n.generators.loc[sel_current, "p_nom_max"] -= existing.loc[
                sel_p_year
            ].rename(lambda x: x[:-4] + current_horizon)

    n.generators.p_nom_max.clip(lower=0, inplace=True)


def prepare_network(n, solve_opts=None, foresight="overnight", config=None):
    if "clip_p_max_pu" in solve_opts:
        for df in (
            n.generators_t.p_max_pu,
            n.generators_t.p_min_pu,
            n.storage_units_t.inflow,
        ):
            df.where(df > solve_opts["clip_p_max_pu"], other=0.0, inplace=True)

    if solve_opts.get("load_shedding"):
        n.add("Carrier", "load")
        n.madd(
            "Generator",
            n.buses.index,
            " load",
            bus=n.buses.index,
            carrier="load",
            sign=1e-3,  # Adjust sign to measure p and p_nom in kW instead of MW
            marginal_cost=1e2,  # Eur/kWh
            # intersect between macroeconomic and surveybased
            # willingness to pay
            # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
            p_nom=1e9,  # kW
        )

    if solve_opts.get("noisy_costs"):
        for t in n.iterate_components():
            # if 'capital_cost' in t.df:
            #    t.df['capital_cost'] += 1e1 + 2.*(np.random.random(len(t.df)) - 0.5)
            if "marginal_cost" in t.df:
                np.random.seed(174)
                t.df["marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df)) - 0.5
                )

        for t in n.iterate_components(["Line", "Link"]):
            np.random.seed(123)
            t.df["capital_cost"] += (
                1e-1 + 2e-2 * (np.random.random(len(t.df)) - 0.5)
            ) * t.df["length"]

    if solve_opts.get("nhours"):
        nhours = solve_opts["nhours"]
        n.set_snapshots(n.snapshots[:nhours])
        n.snapshot_weightings[:] = 8760.0 / nhours

    if foresight == "myopic":
        add_land_use_constraint(n)

    if n.stores.carrier.eq("co2 stored").any():
        limit = config["sector"].get("co2_sequestration_potential", 200)
        add_co2_sequestration_limit(n, limit=limit)

    return n


def add_battery_constraints(n):
    chargers_b = n.links.carrier.str.contains("battery charger")
    chargers = n.links.index[chargers_b & n.links.p_nom_extendable]
    dischargers = chargers.str.replace("charger", "discharger")

    if chargers.empty or ("Link", "p_nom") not in n.variables.index:
        return

    link_p_nom = get_var(n, "Link", "p_nom")

    lhs = linexpr(
        (1, link_p_nom[chargers]),
        (
            -n.links.loc[dischargers, "efficiency"].values,
            link_p_nom[dischargers].values,
        ),
    )

    define_constraints(n, lhs, "=", 0, "Link", "charger_ratio")


def add_chp_constraints(n):
    electric_bool = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("electric")
    )
    heat_bool = (
        n.links.index.str.contains("urban central")
        & n.links.index.str.contains("CHP")
        & n.links.index.str.contains("heat")
    )

    electric = n.links.index[electric_bool]
    heat = n.links.index[heat_bool]

    electric_ext = n.links.index[electric_bool & n.links.p_nom_extendable]
    heat_ext = n.links.index[heat_bool & n.links.p_nom_extendable]

    electric_fix = n.links.index[electric_bool & ~n.links.p_nom_extendable]
    heat_fix = n.links.index[heat_bool & ~n.links.p_nom_extendable]

    link_p = get_var(n, "Link", "p")

    if not electric_ext.empty:
        link_p_nom = get_var(n, "Link", "p_nom")

        # ratio of output heat to electricity set by p_nom_ratio
        lhs = linexpr(
            (
                n.links.loc[electric_ext, "efficiency"]
                * n.links.loc[electric_ext, "p_nom_ratio"],
                link_p_nom[electric_ext],
            ),
            (-n.links.loc[heat_ext, "efficiency"].values, link_p_nom[heat_ext].values),
        )

        define_constraints(n, lhs, "=", 0, "chplink", "fix_p_nom_ratio")

        # top_iso_fuel_line for extendable
        lhs = linexpr(
            (1, link_p[heat_ext]),
            (1, link_p[electric_ext].values),
            (-1, link_p_nom[electric_ext].values),
        )

        define_constraints(n, lhs, "<=", 0, "chplink", "top_iso_fuel_line_ext")

    if not electric_fix.empty:
        # top_iso_fuel_line for fixed
        lhs = linexpr((1, link_p[heat_fix]), (1, link_p[electric_fix].values))

        rhs = n.links.loc[electric_fix, "p_nom"].values

        define_constraints(n, lhs, "<=", rhs, "chplink", "top_iso_fuel_line_fix")

    if not electric.empty:
        # backpressure
        lhs = linexpr(
            (
                n.links.loc[electric, "c_b"].values * n.links.loc[heat, "efficiency"],
                link_p[heat],
            ),
            (-n.links.loc[electric, "efficiency"].values, link_p[electric].values),
        )

        define_constraints(n, lhs, "<=", 0, "chplink", "backpressure")


def basename(x):
    return x.split("-2")[0]


def add_pipe_retrofit_constraint(n):
    """Add constraint for retrofitting existing CH4 pipelines to H2 pipelines."""

    gas_pipes_i = n.links.query("carrier == 'gas pipeline' and p_nom_extendable").index
    h2_retrofitted_i = n.links.query(
        "carrier == 'H2 pipeline retrofitted' and p_nom_extendable"
    ).index

    if h2_retrofitted_i.empty or gas_pipes_i.empty:
        return

    link_p_nom = get_var(n, "Link", "p_nom")

    CH4_per_H2 = 1 / n.config["sector"]["H2_retrofit_capacity_per_CH4"]
    fr = "H2 pipeline retrofitted"
    to = "gas pipeline"

    pipe_capacity = n.links.loc[gas_pipes_i, "p_nom"].rename(basename)

    lhs = linexpr(
        (
            CH4_per_H2,
            link_p_nom.loc[h2_retrofitted_i].rename(index=lambda x: x.replace(fr, to)),
        ),
        (1, link_p_nom.loc[gas_pipes_i]),
    )

    lhs.rename(basename, inplace=True)
    define_constraints(n, lhs, "=", pipe_capacity, "Link", "pipe_retrofit")


def add_co2_sequestration_limit(n, limit=200):
    """
    Add a global constraint on the amount of Mt CO2 that can be sequestered.
    """
    n.carriers.loc["co2 stored", "co2_absorptions"] = -1
    n.carriers.co2_absorptions = n.carriers.co2_absorptions.fillna(0)

    limit = limit * 1e6
    for o in n.opts:
        if "seq" not in o:
            continue
        limit = float(o[o.find("seq") + 3 :]) * 1e6
        break

    n.add(
        "GlobalConstraint",
        "co2_sequestration_limit",
        sense="<=",
        constant=limit,
        type="primary_energy",
        carrier_attribute="co2_absorptions",
    )


def add_EQ_constraints(n, o, scaling=1e-1):
    # Sanity check for EQ option
    if not o[-1] == "c":
        logging.warning(
            "We only support country-level self-sufficiency "
            "constraints! Ignore EQ option."
        )
        return

    # Extract the self-sufficiency ratio
    float_regex = "[0-9]*\.?[0-9]+"
    level = float(re.findall(float_regex, o)[0])

    # For the self-sufficiency ratio, we limit the ratio of imported
    # energy to locally regenerated renewable energy. More
    # specifically, locally produced (renewable) energy includes:
    # - Variable renewables (wind, solar, hydro),
    # - Biogas
    # - Solid biomass
    # - Nuclear power (considered "renewable" here)
    # - Ambient heat (used in heat pumps)
    # On the other hand, imported energy includes:
    # - Electricity (grid, AC & DC)
    # - Hydrogen (pipeline)
    # - Gas
    # For all the above, it's the _net_ import that's consider; export
    # is counted as negative import. Gas has a special role
    # implementation-wise; there is a gas network with pipelines but
    # also "gas input" locations (implemented as PyPSA generators)
    # which can inject (import) as into the network. These import
    # locations represent both LNG terminals, pipeline imports and gas
    # production (on gas rigs). While gas production isn't technically
    # an import, we count it as such here. Rather, we don't count it
    # as local renewable production.

    # NB: We assume that biomass transport is enabled in the sector config! This
    # is needed to get a spatially resolved biomass supply.

    # The self-sufficiency constraint takes the form
    # `(1 - 1/level) * local_production + net_imports <= 0`

    local_factor = 1 - 1 / level

    # Start by building a collection of linear expressions
    # representing local renewable energy input / production for each
    # country.

    # Total variable renewable energy production (by country)
    var_renew_idx = n.generators.loc[
        n.generators.carrier.isin(
            [
                "onwind",
                "offwind-ac",
                "offwind-dc",
                "solar",
                "ror",
                "residential rural solar thermal",
                "services rural solar thermal",
                "residential urban decentral solar thermal",
                "services urban decentral solar thermal",
                "urban central solar thermal",
                "solar rooftop",
            ]
        )
    ].index
    var_renew = (
        linexpr(
            (
                local_factor * n.snapshot_weightings.generators,
                get_var(n, "Generator", "p").loc[:, var_renew_idx].T,
            )
        )
        .T.groupby(
            n.generators.bus.map(n.buses.location).map(n.buses.country), axis="columns"
        )
        .apply(join_exprs)
    )

    # Total hydro production per country (implemented as storage_unit, not generator)
    hydro_idx = n.storage_units.loc[n.storage_units.carrier == "hydro"].index
    hydro = (
        linexpr(
            (
                local_factor * n.snapshot_weightings.generators,
                get_var(n, "StorageUnit", "p_dispatch").loc[:, hydro_idx].T,
            )
        )
        .T.groupby(
            n.storage_units.bus.map(n.buses.location).map(n.buses.country),
            axis="columns",
        )
        .apply(join_exprs)
    )
    hydro = hydro.reindex(var_renew.index, fill_value="")

    # Total biogas and solid biomass production (by country).
    # NB: We assume that biomass transport is enabled! This is needed
    # to get a spatially resolved biomass supply.
    biomass_idx = n.stores.loc[
        n.stores.carrier.isin(
            [
                "biogas",
                "solid biomass",
            ]
        )
    ].index
    # Total biomass (solid and biogas) potential is inserted in the
    # network as initial store capacity. So we calculate the total
    # amount of biomass produced in the country as the difference
    # between the initial and final biomass "store" levels.
    biomass_first_e = get_var(n, "Store", "e").loc[n.snapshots[0], biomass_idx]
    biomass_last_e = get_var(n, "Store", "e").loc[n.snapshots[-1], biomass_idx]
    biomass = (
        pd.concat(
            [
                linexpr((local_factor, biomass_first_e)),
                linexpr((-local_factor, biomass_last_e)),
            ],
            axis="rows",
        )
        .groupby(n.stores.bus.map(n.buses.location).map(n.buses.country))
        .apply(join_exprs)
    )

    # Nuclear power
    nuclear_idx = n.links.loc[n.links.carrier == "nuclear"].index
    nuclear_coeffs = local_factor * pd.DataFrame(
        np.outer(
            n.snapshot_weightings.generators.values,
            n.links.loc[nuclear_idx, "efficiency"].values,
        ),
        index=n.snapshots,
        columns=nuclear_idx,
    )
    nuclear = (
        linexpr(
            (
                nuclear_coeffs,
                get_var(n, "Link", "p").loc[:, nuclear_idx],
            )
        )
        .groupby(n.links.bus1.map(n.buses.country), axis="columns")
        .apply(join_exprs)
    )

    # Ambient heat for heat pumps
    heat_pump_idx = n.links.filter(like="heat pump", axis="rows").index
    # To get the ambient heat extracted, we subtract 1 from the
    # efficiency of the heat pump (where "efficiency" is really COP
    # for heat pumps).
    from_ambient = n.links_t["efficiency"].loc[:, heat_pump_idx] - 1
    ambient_heat_coeffs = local_factor * from_ambient.mul(
        n.snapshot_weightings.generators, axis="rows"
    )
    ambient_heat = (
        linexpr(
            (
                ambient_heat_coeffs,
                get_var(n, "Link", "p").loc[:, heat_pump_idx],
            )
        )
        .groupby(
            n.links.bus1.map(n.buses.location).map(n.buses.country), axis="columns"
        )
        .apply(join_exprs)
    )

    # Total local renewable energy production
    local_energy_prod = var_renew + hydro + biomass + nuclear + ambient_heat

    # Now, we build a collection of linear expressions representing
    # imported energy for each country.

    # Electricity imports
    # Lines crossing borders:
    lines_in = {
        c: n.lines.loc[
            (n.lines.bus1.map(n.buses.country) == c)
            & (n.lines.bus0.map(n.buses.country) != c)
        ].index
        for c in n.buses.country
    }
    lines_out = {
        c: n.lines.loc[
            (n.lines.bus0.map(n.buses.country) == c)
            & (n.lines.bus1.map(n.buses.country) != c)
        ].index
        for c in n.buses.country
    }

    elec_imports = {
        c: linexpr((n.snapshot_weightings.generators, get_var(n, "Line", "s").T))
        .T.loc[:, lines_in[c]]
        .sum()
        .sum()
        if len(lines_in[c]) > 0
        else ""
        for c in local_energy_prod.index
    }
    elec_exports = {
        c: linexpr((-n.snapshot_weightings.generators, get_var(n, "Line", "s").T))
        .T.loc[:, lines_out[c]]
        .sum()
        .sum()
        if len(lines_out[c]) > 0
        else ""
        for c in local_energy_prod.index
    }
    elec_net_imports = {
        c: elec_imports[c] + elec_exports[c] for c in local_energy_prod.index
    }
    elec_net_imports = pd.Series(elec_net_imports)

    # Pipeline imports
    pipeline_carriers = [
        "H2 pipeline",
        "H2 pipeline retrofitted",
        "gas pipeline",
        "gas pipeline new",
        "solid biomass transport",  # Not a pipeline but functionally equivalent
        "DC",  # Ditto
        "Fischer-Tropsch",
        "biomass to liquid",
        "residential rural oil boiler",
        "services rural oil boiler",
        "residential urban decentral oil boiler",
        "services urban decentral oil boiler",
    ]
    into_country = lambda c: (
        (n.links.bus1.map(n.buses.location).map(n.buses.country) == c)
        & (n.links.bus0.map(n.buses.location).map(n.buses.country) != c)
    )
    out_of_country = lambda c: (
        (n.links.bus0.map(n.buses.location).map(n.buses.country) == c)
        & (n.links.bus1.map(n.buses.location).map(n.buses.country) != c)
    )
    links_in = {
        c: n.links.loc[
            into_country(c) & (n.links.carrier.isin(pipeline_carriers))
        ].index
        for c in n.buses.country
    }
    links_out = {
        c: n.links.loc[
            out_of_country(c) & (n.links.carrier.isin(pipeline_carriers))
        ].index
        for c in n.buses.country
    }

    pipeline_imports = {
        c: linexpr((n.snapshot_weightings.generators, get_var(n, "Link", "p").T))
        .T.loc[:, links_in[c]]
        .sum()
        .sum()
        if len(links_in[c]) > 0
        else ""
        for c in local_energy_prod.index
    }
    pipeline_exports = {
        c: linexpr((-n.snapshot_weightings.generators, get_var(n, "Link", "p").T))
        .T.loc[:, links_out[c]]
        .sum()
        .sum()
        if len(links_out[c]) > 0
        else ""
        for c in local_energy_prod.index
    }
    pipeline_net_imports = {
        c: pipeline_imports[c] + pipeline_exports[c] for c in local_energy_prod.index
    }
    pipeline_net_imports = pd.Series(pipeline_net_imports)

    # Gas imports (LNG, pipeline, production)
    gas_import = (
        linexpr(
            (
                n.snapshot_weightings.generators,
                get_var(n, "Generator", "p")
                .loc[:, n.generators.loc[n.generators.carrier == "gas"].index]
                .T,
            )
        )
        .T.groupby(
            n.generators.bus.map(n.buses.location).map(n.buses.country), axis="columns"
        )
        .apply(join_exprs)
    )
    gas_import = gas_import.reindex(local_energy_prod.index).fillna("")

    # Total net imports
    net_imports = elec_net_imports + pipeline_net_imports + gas_import

    # Construct and add the final constraint
    lhs = local_energy_prod + net_imports
    name = "self_sufficiency"
    sense = "<="

    n.add(
        "GlobalConstraint",
        name,
        sense=sense,
        constant=0,
        type=np.nan,
        carrier_attribute=np.nan,
    )

    define_constraints(n, lhs, sense, 0, "GlobalConstraint", "mu", spec=name)


def extra_functionality(n, snapshots):
    add_battery_constraints(n)
    add_pipe_retrofit_constraint(n)

    for o in n.opts:
        if o.startswith("EQ"):
            add_EQ_constraints(n, o)


def solve_network(n, config, opts="", snapshots=None, **kwargs):
    solver_options = config["solving"]["solver"].copy()
    solver_name = solver_options.pop("name")
    cf_solving = config["solving"]["options"]
    track_iterations = cf_solving.get("track_iterations", False)
    min_iterations = cf_solving.get("min_iterations", 4)
    max_iterations = cf_solving.get("max_iterations", 6)
    keep_shadowprices = cf_solving.get("keep_shadowprices", True)

    # add to network for extra_functionality
    n.config = config
    n.opts = opts

    if snapshots is None:
        snapshots = n.snapshots

    if cf_solving.get("skip_iterations", False):
        network_lopf(
            n,
            snapshots,
            solver_name=solver_name,
            solver_options=solver_options,
            extra_functionality=extra_functionality,
            keep_shadowprices=keep_shadowprices,
            **kwargs
        )
    else:
        ilopf(
            n,
            snapshots,
            solver_name=solver_name,
            solver_options=solver_options,
            track_iterations=track_iterations,
            min_iterations=min_iterations,
            max_iterations=max_iterations,
            extra_functionality=extra_functionality,
            keep_shadowprices=keep_shadowprices,
            **kwargs
        )
    return n


if __name__ == "__main__":
    if "snakemake" not in globals():
        from helper import mock_snakemake

        snakemake = mock_snakemake(
            "solve_network",
            weather_year="",
            simpl="",
            opts="",
            clusters="37",
            lv=1.0,
            sector_opts="168H-T-H-B-I-A-solar+p3-dist1",
            planning_horizons="2030",
        )

    logging.basicConfig(
        filename=snakemake.log.python, level=snakemake.config["logging_level"]
    )

    tmpdir = snakemake.config["solving"].get("tmpdir")
    if tmpdir is not None:
        from pathlib import Path

        Path(tmpdir).mkdir(parents=True, exist_ok=True)
    opts = snakemake.wildcards.sector_opts.split("-")
    solve_opts = snakemake.config["solving"]["options"]

    fn = getattr(snakemake.log, "memory", None)
    with memory_logger(filename=fn, interval=30.0) as mem:
        overrides = override_component_attrs(snakemake.input.overrides)
        n = pypsa.Network(snakemake.input.network, override_component_attrs=overrides)

        n = prepare_network(n, solve_opts, snakemake.config["foresight"], snakemake.config)

        n = solve_network(
            n,
            config=snakemake.config,
            opts=opts,
            solver_dir=tmpdir,
            solver_logfile=snakemake.log.solver,
        )

        if "lv_limit" in n.global_constraints.index:
            n.line_volume_limit = n.global_constraints.at["lv_limit", "constant"]
            n.line_volume_limit_dual = n.global_constraints.at["lv_limit", "mu"]

        n.export_to_netcdf(snakemake.output[0])

    logger.info("Maximum memory usage: {}".format(mem.mem_usage))
