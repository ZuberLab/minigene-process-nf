# minigene-process-nf
Process minigene screening data

## Installation

### Nextflow
Install `nextflow` (DSL2) following the instructions at https://www.nextflow.io/docs/latest/getstarted.html

### Apptainer
Install `apptainer` following the instructions at
https://apptainer.org/docs/user/latest/quick_start.html#installation

### minigene-process-nf pipeline
The most convenient way is to install `minigene-process-nf` is to use `nextflow`'s built-in `pull` command
```bash
nextflow pull zuberlab/minigene-process-nf
```

## Test
Before you start, make sure `nextflow` and `apptainer` are properly installed on your system.

Clone git repository from Github and run the pipeline using the provided test data.
```bash
git clone https://github.com/ZuberLab/minigene-process-nf.git
cd minigene-process-nf
./test
```

## Documentation
```bash
nextflow run zuberlab/minigene-process-nf --help
```