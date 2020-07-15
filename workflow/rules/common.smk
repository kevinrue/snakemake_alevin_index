from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

##### print date and time #####

# prints date and time, e.g. '2020-07-14T10:03:08'
DATETIME = "date +'%Y-%m-%dT%H:%M:%S'"
