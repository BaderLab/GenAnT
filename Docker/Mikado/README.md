Instance installs all the required tools and programs to easily run the [mikado pipeline](https://github.com/EI-CoreBioinformatics/mikado).

* Python - 3.10
* [Portucullis](https://github.com/EI-CoreBioinformatics/portcullis)
* [Prodigal](https://github.com/hyattpd/Prodigal)
* [Samtools](https://github.com/samtools/samtools/releases/tag/1.11)
* [SQLite3](https://docs.python.org/3/library/sqlite3.html)

To run the full pipeline you also need to use the [blast docker]()

Uses Ubuntu 22 with Mikado v2.3.2

Primary function of instance is to run mikado.  To use this instance to run mikado run the following:
```
docker run -it -v  "$(pwd)":/global baderlab/mikado:ubuntu22_mikado2.3.2 mikado configure [...]
```

For detailed instructions running mikado see - [mikado documentation](https://mikado.readthedocs.io/en/stable/Tutorial/)}