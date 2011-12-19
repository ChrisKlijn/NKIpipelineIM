# Pipeline for insertional mutagenesis mapping

This pipeline is designed for the mapping of insertional mutagenesis sequencing data from the Roche 454 platform. It currently supports MMTV, Maloney, Sleeping Beauty and PiggyBAC systems.

## Example run

The basic structure of a analysis is a single .fna file per directory. for example, a Sleeping Beauty experiment in the following directory:

    ~/analysis_dir/library.fna

Then run the pipeline in that subdirectory as follows (assuming the pipeline code is located in `pipelineDir`, and the barcode files are located in `barcodeDir`):

    cd ~/analysis_dir/
    pipelineDir/pipeline_SB.sh library.fna analysis1 1 barcodeDir/MMTVbarcodes.fa

The first argument is the library fna file with the 454 reads, the second argument is a basename that is prefixed to all the output files, the 1 indicated that an sqLite database will be constructed with the output data and the last argument is a FASTA file containing the barcodes of the experiment.

These commands will create several output files. The `insert.db` file is the sqLite database with the mapped information, `analysis1.export.txt`  will be the text-based output file that can be further processed using the R code.

The R code contains two files, `functions_preprocessing.R' contains all the code for preprocessing an export file and `example_workflow.R` contains an example how to link the different functions. The result is a text file that can be uploaded to the iMDB.

## To Do and Issues

* Currently there are 4 implementation for each system (PS, SB, MMTV and MuLV). This should be consolidated into 2 files (one .sh and one .pl) with the system being a commandline argument
* There is no other name than `insert.db` for the database. This limits the pipeline to one fna file per directory.

    