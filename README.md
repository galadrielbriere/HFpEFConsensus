# HFpEFConsensus
Consensus clustering for HFpEF samples

## Data preprocessing and input clusterings generation
Run R scripts `data_preprocessing.R` and `make_input_clusterings.R` to preprocess the data and compute input clusterings

## Installing and running ClustOmics

Visit [https://github.com/galadrielbriere/ClustOmics](https://github.com/galadrielbriere/ClustOmics) to learn more about ClustOmics and install the tool. 
Read the associated publication here:
> Brière, G., Darbo, É., Thébault, P. et al. Consensus clustering applied to multi-omics disease subtyping. BMC Bioinformatics 22, 361 (2021). https://doi.org/10.1186/s12859-021-04279-1

In ClustOmics folder, update the configuration file `config.yaml` using the configuration file provided in this repository. Create a `data/HFpEF` folder and copy/paste all files from the `input_clusterings` folder provided in this repository inside the `data/HFpEF` folder.

Run the command:

`snakemake out/HFpEF.HFpEF.FuseClusterings.10_supports.log`

The file `consensusClustering.clst` from this repository store the consensus clustering obtained from running ClustOmics.
