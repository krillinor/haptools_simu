library(docopt)
library(data.table)
library(stringr)

"Filter 1KG by HM3 SNPs.

Usage:
  get_1kg_snplist.R --dir=<dir> [--resource_dir=<resource_dir>]
  get_1kg_snplist.R (-h | --help)

Options:
  --dir=<dir>                    Working directory.
  --resource_dir=<resource_dir>  Directory with resources, f.x. 1KG pgen files. Assumes working dir. if not specified [default: NULL].
  -h --help                      Show this screen.

" -> doc

a <- docopt(doc)

if (a$resource_dir == "NULL") {
    a$resource_dir <- a$dir
}

hm3 <- fread(str_glue("{a$resource_dir}/snpinfo_mult_1kg_hm3_hg38"))
hm3$chr <- as.character(hm3$chr)

pvar_1kg <- fread(str_glue("grep -v \"##\" {a$resource_dir}/all_hg38_qc.pvar"))

hm3_a <- hm3[pvar_1kg, on = c(chr = "#CHROM", pos = "POS", A1 = "REF", A2 = "ALT"), nomatch = 0]
hm3_b <- hm3[pvar_1kg, on = c(chr = "#CHROM", pos = "POS", A1 = "ALT", A2 = "REF"), nomatch = 0]

use <- rbindlist(list(hm3_a, hm3_b))

fwrite(list(use$ID), str_glue("{a$dir}/hm3_use"))
