STORI
=====

Selectable Taxon Ortholog Retrieval Iteratively

User Guide: [https://tinyurl.com/mr395m6z](https://tinyurl.com/mr395m6z)  
Thesis: [http://linkd.in/1fZO63l](http://linkd.in/1fZO63l)

### Introducing single-node STORI

Now STORI should run on any box with bash, Perl, and Python.
The initial 2013 release of STORI (this repo) requires a cluster using the job
scheduler Moab, RHEL 6.5, and assumes that sequences are identified using only their GI number. But the latest release (linked below) runs on a single node, and is compatible with sequence IDs using the accession.version format.
I tested it on CentOS 7.

https://github.com/jgstern/STORI_singlenode
