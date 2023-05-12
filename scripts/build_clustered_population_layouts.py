"""Build clustered population layouts."""

import geopandas as gpd
import xarray as xr
import pandas as pd
import atlite

if __name__ == '__main__':
    if 'snakemake' not in globals():
        from helper import mock_snakemake
        snakemake = mock_snakemake(
            'build_clustered_population_layouts',
            weather_year='',
            simpl='',
            clusters=48,
        )

    cutout_name = snakemake.input.cutout
    year = snakemake.wildcards.weather_year
    if year: cutout_name = cutout_name.format(weather_year=year)
    cutout = atlite.Cutout(cutout_name)

    # Load the onshore model regions
    clustered_regions = gpd.read_file(snakemake.input.regions_onshore).set_index('name')
    # Simplify the geometry
    clustered_regions.geometry = clustered_regions.geometry.buffer(0).squeeze()

    I = cutout.indicatormatrix(clustered_regions)

    pop = {}
    for item in ["total", "urban", "rural"]:
        pop_layout = xr.open_dataarray(snakemake.input[f'pop_layout_{item}'])
        pop[item] = I.dot(pop_layout.stack(spatial=('y', 'x')))

    pop = pd.DataFrame(pop, index=clustered_regions.index)

    pop["country"] = clustered_regions["country"]
    country_population = pop.total.groupby(pop.country).sum()
    pop["fraction"] = pop.total / pop.country.map(country_population)

    pop.to_csv(snakemake.output.clustered_pop_layout)
