# Pipelines overview


`bohra` has 6 pipelines which will use different combinations of modules to generate specific results. 

**if full table not visible please scroll**

| Module | [basic](../pipelines/basic.md) | [preview](../pipelines/preview.md)| [assemble](../pipelines/assemble.md) | [amr_typing](../pipelines/amr_typing.md) | [comparative](../pipelines/comparative.md) | [full](../pipelines/full.md) |[tb](../pipelines/tb.md)|
| :---: |:---: |:---: |:---: |:---: |:---: |:---: |:---: |
|Sequence assessment|Yes|Yes|Yes|Optional|Yes|Yes|Yes|
|Speciation|Optional|Optional|Optional|Optional|Optional|Optional|Optional|
|Assemble| No | No | Yes | Yes (if required) | No | Yes | No|
|MLST| No | No | No | Yes | No | Yes | No |
|Serotype | No | No | No | Yes (if speciation done) | No | Yes | No |
|AMR | No | No | No | Yes | No | Yes | Yes |
| Plasmid detection | No | No | No | Yes | No | Yes | No |
| Tree | No | Yes (default quicktree)| No | No | Yes (default iqtree) |Yes (default iqtree)| Yes (default iqtree)|
|Distance/comparative analysis| No | Yes (mash) |No | No | Yes (default snippy) | Yes (default snippy) | Yes (default snippy) | 
| Pangenome | No | No | No | No | No | Yes  | No  