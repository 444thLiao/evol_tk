# batch run scripts

Normally, it just accept some basically parameters and run some usual commands.

These scripts are mainly for reducing the effort to adjust and re-write some commands.

For example, using `batch_any.py` for batch `prokka` command as follow:

`python batch_any.py -i ./fna_files -o ./prokka_o -s fna -ns '' -np 10 -get_name -cmd "prokka {infile} --outdir {ofile} --force --prefix {name} --locustag {name} --cpus 10 "  `

Be careful that the `cmd` parameter accept a plain text command and doesn't recognize it correctness. The command should be tested with single IO in terminal first.
More detailed descriptions about the parameters of `batch_any.py` could see the output of `python batch_any.py --help`