# kgt_mlst (mlst tool with kmergenetyper)

# Installation

kgt_mlst typer requires kma (>=1.4.9) and kmergenetyper (>=1.0.16) to be installed and usable directly from path (i.e. $kma).
This can be done by installing kma from source or using conda:

`conda create -n kgt_mlst -f kgt_mlst.yml`

`conda activate kgt_mlst`

`pip install kgt_mlst`

# Usage

`kgt_mlst -h`

# Database download

The following command will download the kgt_mlst database to your current working directory:

`kgt_mlst --download_db`