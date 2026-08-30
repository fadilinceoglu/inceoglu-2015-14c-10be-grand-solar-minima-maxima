# Data provenance

## Released analysis inputs

The reproduction consumes the three committed inputs below. The two
solar-modulation-potential (SMP) series share an annual grid from 305 to 8552
years BP relative to 1950 CE, with 8,248 rows each.

| File | SHA-256 | Contents |
| --- | --- | --- |
| `data/derived/smp_14c_aligned.dat` | `404219e6723480b7b9166e2185874c33954e63c97c9c9dbc7a9598cae2fe9538` | Author-derived annual 14C SMP on the common analysis grid |
| `data/derived/smp_10be_aligned.dat` | `563da1692bfdf6279370cfcf1c57038fa01f364d6a36b785707fd05ab2e74215` | Author-derived, chronology-aligned annual 10Be SMP on the common analysis grid |
| `catalog/manual_event_catalog.csv` | `a61021bd9588a58f5da65544b7c0bdf0273e50fed12a47586afe3546e960694b` | The 160 visually identified paired events used by the study |

The catalogue has 75 low-activity events and 85 high-activity events. Each row
records the onset, end, center, duration, standardized amplitude, and sign for
both nuclides. The rows are paired in the catalogue's declared order.

These files are analysis inputs, not hardcoded table or figure values. The
pipeline verifies their checksums before calculating the published outputs.

## Scientific lineage and scope

The 14C series is derived from the IntCal13 radiocarbon record, and the 10Be
series is derived from the GRIP beryllium-10 record. Knudsen et al. (2009), DOI
[`10.1029/2009GL039439`](https://doi.org/10.1029/2009GL039439), describes the
chronology-alignment method: the radiocarbon chronology is used as the
reference and the GRIP timescale is adjusted using maxima of running
cross-correlations.

This repository reproduces the downstream analysis from the released aligned
SMP inputs. It does not implement the complete physical inversion from raw
observations to SMP or the original running-cross-correlation alignment. The
parameters and software needed to reproduce those upstream transformations are
not included, so the public observations listed below must not be substituted
for the committed analysis inputs.

## Public upstream sources

`scripts/download_public_data.py` downloads the following public files,
verifies their expected checksums, and records response metadata:

| File | Authoritative source | Expected SHA-256 | Role |
| --- | --- | --- | --- |
| `intcal13.14c` | <https://www.intcal.org/curves/intcal13.14c> | `0af281d2b559143ea9230c45e3c709137d868b8dbc09abf1ad16739fdb424cbd` | IntCal13 atmospheric calibration curve named in the paper |
| `intcal04.14c` | <https://www.intcal.org/curves/intcal04.14c> | `e119d6235e11a08414aa75aae66db5d2ccfa6dfd75c77ce874f6af1d1c7838f8` | IntCal04 comparison used by the Knudsen processing chain |
| `grip_10be.txt` | [NOAA/NCEI GRIP archive](https://www.ncei.noaa.gov/pub/data/paleo/icecore/greenland/summit/grip/cosmoiso/grip_10be.txt) | `ad70f6be5c11adee1db1484c51969dc7497c538580f69d56fc3db744ae521ddb` | Upstream GRIP measurements |
| `grip_10be.xls` | [NOAA/NCEI GRIP archive](https://www.ncei.noaa.gov/pub/data/paleo/icecore/greenland/summit/grip/cosmoiso/grip_10be.xls) | `cb0d30b0b70932a91197fba59c27cb481475cfc83ef342afaa7dcac56f08ad2d` | Spreadsheet form of the GRIP archive |

Relevant primary citations are Reimer et al. (2013), DOI
[`10.2458/azu_js_rc.55.16947`](https://doi.org/10.2458/azu_js_rc.55.16947), and
the NCEI dataset DOI
[`10.25921/hj8n-xc03`](https://doi.org/10.25921/hj8n-xc03). IntCal04 is
associated with DOI
[`10.1017/S0033822200032999`](https://doi.org/10.1017/S0033822200032999).

The downloader selects 300–8550 BP from IntCal13 and linearly interpolates
Delta14C and its uncertainty to annual resolution. It then stops before the
carbon-cycle production calculation, solar-modulation inversion, and
chronology alignment. Optional downloads are written only below ignored
`data/raw/` and `data/processed/` paths and never overwrite `data/derived/`.

## Reproduction boundary

Regenerating the committed SMP inputs from the public observations would
require the processed GRIP age, accumulation, and flux series; the
nuclide-production and inversion implementation; the geomagnetic-dipole input;
the carbon-cycle production model; and the cross-correlation alignment
parameters. These are outside the scope of this repository. The causal boundary
is explicit: public observations are upstream references, while the committed
SMP files are the inputs to the reproduced analysis.

## Rights

The IntCal download page does not state an express dataset license. The NCEI
GRIP record requires the original dataset citation, the downloaded subset and
access date, and the archive DOI. No separate public license is stated for the
author-derived aligned arrays or event catalogue. The repository's MIT license
covers the original Python software and documentation only; path-level terms
and credits are listed in `DATA_NOTICE.md`.
