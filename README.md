# Fork of PyPSA-Eur-Sec

This is a fork of PyPSA-Eur-Sec, with an updated verion of the `multiyear` branch.

## Implementation notes

The multiyear branch is at the current time not 100% functional out of the box.
1. As noted at https://github.com/PyPSA/pypsa-eur-sec/pull/75, this branch needs the 2021 version of eurostat data, which isn't included in the Zenodo dataset and must be downloaded manually. To do this, download [this](https://ec.europa.eu/eurostat/documents/38154/4956218/Energy-balance-sheets-February-2021-edition.zip/4b1d6665-f303-be7d-a7e5-1e0da16ec0d9?t=1612709565471) zip file and extract the contents into `data/eurostat-energy_balances-june_2021_edition`. This is needed for the `build_energy_totals` rule. However, note that the `build_industrial_production_per_country` also needs eurostat data, but _has not been updated_ to use the 2021 version. Therefore, both the 2021 version (downloaded manually) _and_ the 2018 version (in the Zenodo dataset) are needed for now. (Update: 2021 version of eurostat data is included in git repo for now.)
