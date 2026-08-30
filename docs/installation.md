# Installation



## Recommended (conda or mamba)

**1. Install conda (skip this step if you already have conda installed)**

If you do not already have `conda` installed, you can check out the documentation [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html). We recommend you install [`miniforge`](https://github.com/conda-forge/miniforge)


**2. Create and install bohra**

```
mamba (or conda) create -n bohra -c bioconda bohra
```

**3. Install dependencies**

```
conda activate bohra
bohra deps install
```

**4. Set environment variables (optional)**

_kraken2_
You will need to ensure that you have a kraken2 database available. If you do not set this as an environment variable you will need to provide the path to the kraken2 database you wish to use at runtime. **IMPORTANT - if you do not have a kraken2 database set species specific functions of `bohra` will not work.**

You may already have a kraken2 database environment variable set on your system. To check this 
```
echo $KRAKEN2_DEFAULT_DB
```
If this outputs a path - this is the default kraken2 database.

If this is empty or if you wish to use a different database you can set it at runtime OR set the path in your `bohra` environment.

```
conda env config vars set KRAKEN2_DEFAULT_DB=/path/to/kraken2_db
```
it will tell you to deactivate and re-activate your environment

```
conda deactivate bohra
conda activate bohra
```

Now you can test that the environment variable has been set
```
echo $KRAKEN2_DEFAULT_DB
```
should print 
```
/path/to/kraken2_db
```
_mlst_

If you have your own `mlst` database you can provide this to `bohra`. Similar to setting the `KRAKEN2_DEFAULT_DB`.


You may already have a kraken2 database environment variable set on your system. To check this 
```
echo $MLST_DBDIR
```
If this outputs a path - this is the default mlst database.

If this is empty or if you wish to use a different database you can set it at runtime OR set the path in your `bohra` environment.

```
conda env config vars set MLST_DBDIR=/path/to/mlst_db
```
it will tell you to deactivate and re-activate your environment

```
conda deactivate bohra
conda activate bohra
```

Now you can test that the environment variable has been set
```
echo $MLST_DBDIR
```
should print 
```
/path/to/mlst_db
```


**5. Testing your installation**

If you would like to make sure that all the dependencies are installed properly and that the pipeline will run on your system you can run


```
bohra test --cpus X
```

(where X is the number of cpus you want to use)

This will download some fastq files and run the full pipeline.


## Install from source


**Setup the `bohra` environment and install the pipeline co-ordinator**

1. Clone the `bohra` repository

```
git clone https://github.com/MDU-PHL/bohra.git
```


. Create the environment and install `bohra`


```
cd bohra
conda env create -f environment.yml -n bohra
conda activate bohra
pip3 install .
```

**Install dependencies**


This step will setup conda environments under the `path/to/conda/envs/bohra`(depending on how you have configured your `conda` installation). 


1. Activate your environment from step 2 above
```
conda activate bohra
```

2. Run the bohra dependency installation

```
bohra  deps install
```

## Updating `bohra`

Updates to the `bohra` pipeline  will often include 2 steps/

1. **Update `bohra` env**

```
conda update -n bohra --all
```

This will update the bohra pipeline and running scripts. But not any dependencies

2. **Update any changed dependencies**

```
conda activate bohra
bohra deps update
```

This will update any dependencies that have pinned versions or any new dependencies and should be run each time you update `bohra`.

### Other update options

You may also update a specific tool or dependency environment. This should be done with caution as updates that have not been adeqautely incorporated into the pipeline may cause unexpected behaviour or breakages. 

```
conda activate bohra
bohra deps update --tool <env_name> --force
```

If you would like another tool or different versions please leave an issue and we will endeavour to facillitate it.

If you need to re-install your dependencies for some reason you can run

```
conda activate bohra
bohra deps update --force
```
This will overwrite all existing dependencies and may take some time.