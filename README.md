# Inceoglu 2015: grand solar minima and maxima from 14C and 10Be

This repository reproduces the downstream analysis in Inceoglu et al. (2015). The study compared solar-modulation-potential reconstructions derived from IntCal13 radiocarbon and GRIP beryllium-10 over their common interval, identified simultaneous low- and high-activity episodes, and examined their durations, waiting times, and 20–40 year variability.

The event catalogue contains the 160 paired peaks and dips identified visually for the study. The Python pipeline fits the two-component Gaussian amplitude models, applies the documented amplitude and duration criteria, reconstructs Tables 1–5 and Figures 1–8, and runs numerical tests against the publisher paper. It does not replace the visual catalogue with a new automated detector.

## Citation

F. Inceoglu, R. Simoniello, M. F. Knudsen, C. Karoff, J. Olsen, S. Turck-Chièze, and B. H. Jacobsen, “Grand solar minima and maxima deduced from 10Be and 14C: magnetic dynamo configuration and polarity reversal,” *Astronomy & Astrophysics* **577**, A20 (2015). [https://doi.org/10.1051/0004-6361/201424212](https://doi.org/10.1051/0004-6361/201424212)

## Reproduce

Create an environment from the exact pins in `requirements.lock`, then run from a fresh clone:

```bash
python3 scripts/reproduce.py
```

That one command verifies the two committed aligned solar-modulation-potential inputs, creates every published table and figure under `outputs/`, writes a numerical validation report, and runs the test suite.

To retrieve the authoritative public upstream archives and record response metadata and checksums separately:

```bash
python3 scripts/download_public_data.py
```

The downloader is an optional public-source snapshot tool. It writes only to ignored `data/raw/` and `data/processed/` paths and never replaces the committed SMP inputs. It stops before the physical inversion and chronology-alignment steps, which are not implemented in this repository. See [reproduction notes](docs/reproduction_notes.md) and [data provenance](docs/data_provenance.md).

## Repository contents

- `catalog/manual_event_catalog.csv`: the study's 160 final visual event pairs.
- `data/derived/`: the two aligned annual SMP series actually consumed by the paper reproduction.
- `src/grand_solar_reproduction/`: Python implementation.
- `scripts/reproduce.py`: complete local reproduction and validation command.
- `scripts/download_public_data.py`: fail-closed public-source downloader.
- `outputs/`: eight PDF figures, five CSV tables, and validation reports.
- `docs/`: data provenance, reproduction method, and known scientific differences.

Raw observational downloads are available from the cited upstream sources and
are not committed. Original software and repository documentation are licensed
under the [MIT License](LICENSE). That grant does not cover scientific data,
catalogues, or generated research outputs except where
[`DATA_NOTICE.md`](DATA_NOTICE.md) expressly says otherwise. The MIT entries in
`pyproject.toml` and `CITATION.cff` describe the software package, not every
file in the repository.
