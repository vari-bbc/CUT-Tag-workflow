# CUT&Tag workflow

## Usage

**NOTE** This workflow is optimized for the HPC @ Van Andel Institute.

### Step 1: Configure the workflow

* Move your fastq sequences to `raw_data/`
* Look over the config file under bin/, make sure the references are correct. Do not change modules. 
* Set up the sample sheet (samples.tsv). You may find it helpful to run `./make_samples_template.sh` from inside the `bin/` subdirectory to get a template file based on the files in `raw_data/`. The columns are:
    * **sample**        - Name of the sample. Fastqs will be renamed to this.
    * **control**         - You can  use 'NA'. Leaving it blank should also be fine. 
    * **fq1**           - Name of read1 fastq
    * **fq2**           - Name of read2 fastq
    * **sample_group**            - You can use the same information as sample column.
* If you used make_samples_template.sh, you need to rename the generated template file to "samples.tsv". 

### Step 2: Run the workflow

* From the root directory of the project (where Snakefile is located), run `sbatch bin/run_snakemake.sh`.
