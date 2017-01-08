# bionella
Legionella genome annotation

In the following you can replace `$IMAGE` by:
- __swissprot_blast__: docker image with with Swiss-Prot database
- __uniprot_blast__: docker image with Swiss-Prot and TrEMBL database

See the difference between these databases [here](http://www.uniprot.org/downloads).

### Build docker image
```bash
$ cd Dockerfiles/
$ docker build -t vitorenesduarte/$IMAGE -f $IMAGE.
```

### Run docker image
```bash
$ docker run -e QUERY_DIR=/root/query_dir \
             -v $PWD/query_dir:/root/query_dir \
             -ti vitorenesduarte/swissprot_blast
```

### Useful links:

- [http://fenyolab.org/ibb2014/tutorials/week3/BLAST%20with%20BioPython.pdf](http://fenyolab.org/ibb2014/tutorials/week3/BLAST%20with%20BioPython.pdf)
- [https://www.biostars.org/p/82059/](https://www.biostars.org/p/82059/)
- [http://www.uniprot.org/help/programmatic_access9](http://www.uniprot.org/help/programmatic_access)
