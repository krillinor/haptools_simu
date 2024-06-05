library(docopt)
library(data.table)
library(stringr)

"Liftover HM3 SNPs to hg38.

Usage:
  hm3_to_hg38.R --dir=<dir> [--resource_dir=<resource_dir>]
  hm3_to_hg38.R (-h | --help)

Options:
  --dir=<dir>                    Working directory.
  --resource_dir=<resource_dir>  Directory with resource. Assumes working dir. if not specified [default: NULL].
  -h --help                      Show this screen.

" -> doc

a <- docopt(doc)

if (a$resource_dir == "NULL") {
    a$resource_dir <- a$dir
}

d <- fread(str_glue("{a$resource_dir}/snpinfo_mult_1kg_hm3"))
setnames(d, old = c("CHR", "BP"), new = c("chr", "pos"))

dd <- bigsnpr::snp_modifyBuild(d, str_glue("{a$dir}/liftOver"), from = "hg19", to = "hg38")

fwrite(dd[!is.na(pos), ], str_glue("{a$dir}/snpinfo_mult_1kg_hm3_hg38"), sep = "\t")
