# bionella
Legionella genome annotation

### Build docker image
```bash
$ cd Dockerfiles/
$ docker build -t vitorenesduarte/local_blast -f local_blast .
```

### Run docker image
```bash
$ docker run -it vitorenesduarte/local_blast
```

### Local Database alias
```bash
$ blastdb_aliastool -dblist "nr.00 nr.01 nr.02" -dbtype prot -title nr -out nr
```

### Useful links:

- [http://fenyolab.org/ibb2014/tutorials/week3/BLAST%20with%20BioPython.pdf](http://fenyolab.org/ibb2014/tutorials/week3/BLAST%20with%20BioPython.pdf)
- [https://www.biostars.org/p/82059/](https://www.biostars.org/p/82059/)
- [http://www.uniprot.org/help/programmatic_access9](http://www.uniprot.org/help/programmatic_access)
