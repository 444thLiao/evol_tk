
### NCBI converter
All converter functions actually has been embedded into a single object namely `NCBI_convertor`.

For retrieving metadata from bioproject and biosample of NCBI, you could use below command

`python3 ncbi_convertor/pid2bio.py -i assembly_ids.list -o metadata.csv -f -redo`

Except for the normally and regular parameters like `-i` `-o`. There are some special parameteters.

like
> `-s`, You could retrieve metadata directly from protein id by passing it.
> `-f` force to overwrite output
> `-redo` This script will use some cache for quick re-analysis. You could **pass it to avoid some bugs**.

The above one is the major one.
Besides that, It has some extra toolkits for extract some other info from NCBI.
Such as,
`pid2GI.py` with proteins ids, it will retrieve GI of its. GI is the unique ID defined by NCBI.
`pid2tax.py` with proteins ids, it will retrieve taxonomy ID of its relative organism.
`pid2genome.py` with proteins ids, it will retrieve genome assembly id and it could countinue pass it to `pid2bio.py`
`pid2basic.py` with proteins ids, it will do all three above scripts staff and combined them into one tab.
`pid2all.py`, beside the basic information, it will do the `pid2bio.py` staffs.


TODO: It doesn't try using nucleotide directly...need to test





