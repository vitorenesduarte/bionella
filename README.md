# bionella
Legionella genome annotation

### Prerequesites

- [pandas](http://pandas.pydata.org/)
```bash
$ conda install pandas
```

To start the Jupyter Notebook:
```bash
$ make start-nb
```

To generate the site:

```bash
$ make site
```

To run each of the steps __first__, __second__ and __third__:

```bash
$ make first
$ make second
$ make third
```


### Build docker image
```bash
$ cd Dockerfiles/
$ docker build -t vitorenesduarte/swissprot_blast -f swissprot_blast.
```

This will build a docker image with blast.
This docker image only has one database called __swissprot__.
Check [here](http://www.uniprot.org/downloads) more details
about UniProtKB Swiss-Prot.

### Push the docker image
```
$ docker push vitorenesduarte/swissprot_blast
```

### Run blast in docker

Place in a folder, for example __.query_dir__, a set of
FASTA files and run:

```bash
$ docker run -e QUERY_DIR=.query_dir \
             -e DB=swissprot -v $PWD/.query_dir:/.query_dir \
             -ti vitorenesduarte/swissprot_blast
```

In the end, for each of the FASTA files, you should have
a correspondent XML file.

### Useful links:

- [http://fenyolab.org/ibb2014/tutorials/week3/BLAST%20with%20BioPython.pdf](http://fenyolab.org/ibb2014/tutorials/week3/BLAST%20with%20BioPython.pdf)
- [https://www.biostars.org/p/82059/](https://www.biostars.org/p/82059/)
- [http://www.uniprot.org/help/programmatic_access](http://www.uniprot.org/help/programmatic_access)
- [https://github.com/jdrudolph/uniprot](https://github.com/jdrudolph/uniprot)
- [http://web.expasy.org/blast/blast+_help.html#prog_access](http://web.expasy.org/blast/blast+_help.html#prog_access)
