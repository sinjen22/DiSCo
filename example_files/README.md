DiSCo - examples
=============


Example 1:

```bash
  perl DiSCo.pl -i A.vinosum.faa -o  
```

This command will run DiSCo against the file [A.vinosum.faa](input/A.vinosum.faa "Allochomatium vinosum") and automatically extract the sequences of all [DiSCo hits](output/A.vinosum.DiSCo.filtered.txt "Allochomatium vinosum - hits") in FASTA format ([A.vinosum.DiSCo.faa](output/A.vinosum.DiSCo.faa "Allochomatium vinosum - seq"))

 
- - - -


Example 2:

```bash
  perl DiSCo.pl -i D.vulgaris.faa -s 2 -p D.vulgaris_sep_2 -o  
```

This command will run DiSCo against the file [D.vulgaris.faa](input/D.vulgaris.faa "Desulfovibrio vulgaris") and automatically extract the sequences of all [DiSCo hits](output/D.vulgaris_sep_2.DiSCo.filtered.txt "Desulfovibrio vulgaris - hits") in FASTA format. This time  the separator is changed to ; and a file name prefix is used, which is also used for the [extracted sequences file](output/D.vulgaris_sep_2.DiSCo.faa "Desulfovibrio vulgaris - seq").


- - - -


Example 3:

```bash
  perl DiSCo.pl -i E.coli.faa 
```




This command will run DiSCo against the file [E.coli.faa](input/E.coli.faa "Escherichia coli"). *E. coli* does not have any Dsr-dependent dissimilatory sulfur metabolism related proteins and the DiSCo tool will stop with a warning: 

```bash
  E.coli.DiSCo.txt No significant hits
```
