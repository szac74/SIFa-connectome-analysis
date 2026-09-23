# SIFa connectome analysis

Reproducibility code for the connectome analyses supporting the SIFamide (SIFa) study.

## Scope and provenance

The original connectome analysis was conducted interactively using custom Cypher queries in neuPrint and stepwise Python-based mesh-proximity calculations. A single monolithic analysis script was not archived at the time of the original analysis.

This repository consolidates the final query logic, neuron identifiers, filtering rules, representative-site metadata, and mesh-distance calculations into reproducible scripts corresponding to the analyses reported in the manuscript.

The repository reproduces:

1. High-precision SIFa connectivity summaries for the female Hemibrain (`hemibrain:v1.2.1`) and male MaleCNS (`male-cns:v1.0`) datasets.
2. SIFa-to-target synapse counts used for Supplementary Tables 3 and 4.
3. MaleCNS mesh-proximity calculations used for representative SIFa-Kenyon-cell sites in Fig. 7.

Final EM interpretation of representative sites was performed manually in Neuroglancer using orthogonal views, segmentation overlays, and consecutive raw EM sections. Mesh-derived distances were used only as proximity-screening metrics and should not be interpreted as exact membrane-to-membrane distances or, by themselves, as evidence of synaptic connectivity.

## Datasets

- Female connectome: `hemibrain:v1.2.1` on `neuprint.janelia.org`
- Male connectome: `male-cns:v1.0` on `neuprint.janelia.org`
- MaleCNS segmentation mesh: `precomputed://gs://flyem-male-cns/v1.0/segmentation`
- MaleCNS native voxel size: 8 x 8 x 8 nm

## Repository layout

```text
SIFa-connectome-analysis/
|-- README.md
|-- CODE_AVAILABILITY_TEXT.md
|-- requirements.txt
|-- requirements-lock.txt
|-- .gitignore
|-- config/
|   |-- neuron_ids.csv
|   `-- expected_summary_counts.csv
|-- metadata/
|   `-- fig7_representative_sites.csv
|-- scripts/
|   |-- 01_synapse_counts.py
|   |-- 02_mesh_proximity.py
|   `-- 03_validate_expected_counts.py
`-- outputs/
    |-- connectivity_summary_all.csv
    |-- fig7M_mesh_proximity.json
    |-- fig7O_site_metadata_reproduced.csv
    |-- fig7O_soma_associated_proximity.csv
    |-- grouped_summary_hemibrain_v1.2.1.csv
    |-- grouped_summary_male-cns_v1.0.csv
    |-- neuprint_hp_thresholds.csv
    |-- raw_connections_hemibrain_v1.2.1.csv
    |-- raw_connections_male-cns_v1.0.csv
    |-- supplementary_table3_reproduced.csv
    |-- supplementary_table4_reproduced.csv
    `-- validation_report.csv
```

## Installation

The validated environment used for the reproducibility rerun was Python 3.13.11.

Create a virtual environment and install the required packages:

```bash
python -m venv .venv
pip install -r requirements.txt
```

`requirements-lock.txt` records the complete package environment used for the validated rerun.

## neuPrint authentication

Do not place a neuPrint token in the scripts or commit it to GitHub.

Set the token as an environment variable.

Linux/macOS:

```bash
export NEUPRINT_APPLICATION_CREDENTIALS='YOUR_NEUPRINT_TOKEN'
```

Windows PowerShell:

```powershell
$env:NEUPRINT_APPLICATION_CREDENTIALS='YOUR_NEUPRINT_TOKEN'
```

## 1. Reproduce connectivity summaries

Run:

```bash
python scripts/01_synapse_counts.py
```

The script:

- queries the four SIFa neurons in each dataset;
- retrieves both SIFa-to-target and target-to-SIFa `ConnectsTo` relationships;
- retains connections with `weightHP > 0`;
- records neuPrint high-precision threshold information;
- groups PDF-positive lateral neurons and mushroom-body-associated neuronal classes;
- reports male left and right sides separately and also provides bilateral sums;
- writes raw connection tables and grouped summaries to `outputs/`.

Then run:

```bash
python scripts/03_validate_expected_counts.py
```

A successful validation ends with:

```text
VALIDATION PASSED: all configured key totals matched.
```

## 2. Reproduce MaleCNS mesh-proximity calculations

First, optionally verify that the required mesh manifests are accessible:

```bash
python scripts/02_mesh_proximity.py --probe
```

Reproduce the full-resolution Fig. 7M analysis:

```bash
python scripts/02_mesh_proximity.py --panel M --lod 0
```

Reproduce the full-resolution Fig. 7O analysis:

```bash
python scripts/02_mesh_proximity.py --panel O --lod 0
```

LOD values greater than 0 are intended only for diagnostic testing and should not be used for manuscript distance measurements.

Full-resolution LOD0 calculations can require substantial RAM and computation time.

## Reproducibility validation

The consolidated workflows were rerun and validated against the values reported in the manuscript.

### neuPrint connectivity analysis

All configured key totals were reproduced exactly.

These included:

- female Hemibrain SIFa output to PDF-positive lateral neurons: 16
- male MaleCNS right-side SIFa output to PDF-positive lateral neurons: 15
- male MaleCNS left-side SIFa output to PDF-positive lateral neurons: 3
- female Hemibrain SIFa output to alpha/beta KCs: 169
- female Hemibrain SIFa output to gamma KCs: 37
- female Hemibrain SIFa output to alpha-prime/beta-prime KCs: 9
- male MaleCNS bilateral SIFa output to alpha/beta KCs: 99
- male MaleCNS bilateral SIFa output to gamma KCs: 13
- male MaleCNS bilateral SIFa output to alpha-prime/beta-prime KCs: 5

### Fig. 7M mesh proximity

Using full-resolution LOD0 meshes:

- SIFa body ID: `512245`
- KC body ID: `535798`
- reproduced minimum vertex-to-vertex distance: `19.595918 nm`

This value exactly reproduced the value obtained in the original analysis.

### Fig. 7O soma-associated proximity

For SIFa body ID `10342` and KC body ID `535798`, KC mesh vertices within 3 um of the annotated soma coordinate were defined as the soma-associated region.

The reproduced minimum distances between local SIFa mesh regions and the KC soma-associated region were:

- 0.5-um local radius: 2.566929 um
- 1.0-um local radius: 2.356626 um
- 1.5-um local radius: 2.201874 um

The distance between the selected SIFa-associated center and the annotated KC soma coordinate was:

- 5.936550 um

These values are proximity-screening measurements and do not establish direct soma contact or synaptic connectivity.

## Representative Fig. 7 sites

`metadata/fig7_representative_sites.csv` records body IDs and coordinates used to relocate representative sites.

- Fig. 7K: neuPrint-annotated SIFa-to-l-LNv synaptic site.
- Fig. 7M: mesh-screened close apposition; not classified as a synapse.
- Fig. 7O: SIFa process in the vicinity of a KC soma-associated region; not a demonstrated direct soma contact.

## Reproducibility note

The scripts in this repository consolidate the final analysis logic after the original interactive neuPrint/Cypher workflow. The consolidated workflows were subsequently rerun and verified against the reported connectivity totals and Fig. 7 mesh-proximity measurements.

## Archival recommendation

For publication, a tagged GitHub release can additionally be archived in a DOI-minting repository such as Zenodo.