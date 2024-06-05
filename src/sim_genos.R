library(docopt)
library(data.table)
library(stringr)
library(fs)

"Simulate genotypes.

Usage:
  sim_genos.R --dir=<dir> --ids=<ids> --chrom=<chrom> [--haptools_seed=<haptools_seed> --haptools_path=<haptools_path> --resource_dir=<resource_dir>]
  sim_genos.R (-h | --help)

Options:
  --dir=<dir>                      Working directory.
  --ids=<ids>                      IDs from model files, e.g. ceu,yri.
  --chrom=<chrom>                  Chromosome (1-22).
  --haptools_seed=<seed>           haptools seed [default: 5].
  --haptools_path=<haptools_path>  haptools path [default: haptools].
  --resource_dir=<resource_dir>    Directory with resources, f.x. 1KG pgen files. Assumes working dir. if not specified [default: NULL].
  -h --help                        Show this screen.

" -> doc

a <- docopt(doc)
a$ids <- str_split(a$ids, ",")[[1]]
a$chrom <- as.numeric(a$chrom)
a$haptools_seed <- as.numeric(a$haptools_seed)

if (a$resource_dir == "NULL") {
    a$resource_dir <- a$dir
}

regions <- fread(str_glue("{a$dir}/region_map"), col.names = "region_id", header = FALSE)
regions <- regions[str_detect(region_id, str_glue("^{a$chrom}:")), ]

map_path <- str_glue("{a$resource_dir}/map")
model_path <- str_glue("{a$dir}/model_{a$id}.dat")
pgen_path <- str_glue("{a$resource_dir}/all_hg38_qc.pgen")
sample_info_path <- str_glue("{a$resource_dir}/1000genomes_sampleinfo.tsv")

sim_genos <- function(region, out) {
    cmd <- str_glue("{a$haptools_path} simgenotype --model {model_path} --region {region} --mapdir {map_path} --ref_vcf {pgen_path} --sample_info {sample_info_path} --pop_field --out {out} --seed {a$haptools_seed}")
    system(cmd)
}

for (id in a$ids) {
    out_dir <- str_glue("{a$dir}/{id}_genos")
    dir_create(out_dir)
    for (i in 1:nrow(regions)) {
        regions_i <- regions[i]
        cat(str_glue("{id}_{regions_i}"), "\n", sep = "")
        out <- str_glue("{out_dir}/{id}_{regions_i}_hseed{a$haptools_seed}.pgen")
        sim_genos(region = regions_i, out = out)
    }
}
