"""Build population-weighted energy and heat totals."""

import pandas as pd

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_population_weighted_energy_totals',
            weather_year='',
            simpl='',
            clusters=37,
        )

    config = snakemake.config["energy"]
    data_year = int(config["energy_totals_year"])
    if snakemake.wildcards.weather_year and snakemake.wildcards.kind == 'heat':
        data_year = int(snakemake.wildcards.weather_year)

    pop_layout = pd.read_csv(snakemake.input.clustered_pop_layout, index_col=0)

    totals = pd.read_csv(snakemake.input.totals, index_col=[0,1])
    totals = totals.xs(data_year, level='year')

    # Since our model nodes may consist of multiple countries (i.e.
    # pop_layout.country can look like "AA_BB_CC"), we need to
    # aggregate the energy totals to this level before proceeding. The
    # aggregation is by summing up the energy totals, except the
    # district heating share column, which should be averaged.
    node_countries = totals.index.map(lambda c: next(x for x in pop_layout.country if c in x))
    agg_funcs = {c: "sum" for c in totals.columns}
    if "district heat share" in agg_funcs:
        agg_funcs["district heat share"] = "mean"
    totals = totals.groupby(node_countries).agg(agg_funcs)

    nodal_totals = totals.loc[pop_layout.country].fillna(0.)
    nodal_totals.index = pop_layout.index
    nodal_totals = nodal_totals.multiply(pop_layout.fraction, axis=0)

    nodal_totals.to_csv(snakemake.output[0])
