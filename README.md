# radBarCoder
scripts to extract, align, and type mtDNA data from restriction site associated DNA

---

## How it works

1. map rad data to mtDNA genome using [dDocentHPC mkBAM](https://github.com/cbirdlab/dDocentHPC)
  * obtain reference genome from [NCBI GenBank](https://www.ncbi.nlm.nih.gov/genbank/)
    * reference genome should be a `fasta` formatted file and can be composed of 1, several, or all loci in the mtGenome
  * prepare your `fastq` files for mapping as described in the [dDocentHPC README.md](https://github.com/cbirdlab/dDocentHPC)

2.
