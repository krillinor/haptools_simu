library(docopt)
library(data.table)
library(stringr)
library(fs)

"Simulate phenotypes.

Usage:
  sim_phenos.R --dir=<dir> --ids=<ids> [--p=<p> --h2=<h2> --seed=<seed> --haptools_seed=<haptools_seed> --haptools_path=<haptools_path>]
  sim_phenos.R (-h | --help)

Options:
  --dir=<dir>                      Working directory.
  --ids=<ids>                      IDs from model files, e.g. ceu,yri.
  --p=<p>                          Polygenicity [default: 0.001,0.01,0.1,0.5].
  --h2=<h2>                        Heritability [default: 0.05,0.1,0.25,0.5,0.8].
  --seed=<seed>                    R seed [default: 159]
  --haptools_seed=<seed>           haptools seed [default: 5]
  --haptools_path=<haptools_path>  haptools path [default: haptools]
  -h --help                        Show this screen.

" -> doc

a <- docopt(doc)
a$ids <- str_split(a$ids, ",")[[1]]
a$p <- str_split(a$p, ",")[[1]]
a$h2 <- str_split(a$h2, ",")[[1]]
a$seed <- as.numeric(a$seed)
a$haptools_seed <- as.numeric(a$haptools_seed)

g <- expand.grid(p = a$p, h2 = a$h2)

# could also use 2*f_j*(1-f_j) in the denominator
# where f_j is the freq. of variant j
# since rarer variants should have larger effect sizes
# maybe complicated, since freq. not same across genetic ancestry groups...
sim_phenos <- function(dir = dir, p = p, h2 = h2, id = id, seed = 159, haptools_seed = 5) {
    out_causal <- str_glue("{dir}/betas/causal_p{p}_h2{h2}_seed{seed}.snplist")
    betas <- fread(out_causal, col.names = c("snp", "beta"), header = FALSE)
    out_causal2 <- str_glue("{dir}/betas/{id}_causal_p{p}_h2{h2}_seed{seed}.snplist")
    pgen_prefix <- str_glue("{dir}/{id}_genos/{id}_hseed{haptools_seed}")
    pvar_ids <- fread(str_glue("grep -v \"##\" {pgen_prefix}.pvar"))
    betas2 <- betas[snp %in% pvar_ids$ID, ]
    fwrite(betas2, out_causal2, sep = "\t", col.names = FALSE)

    pgen <- str_glue("{pgen_prefix}.pgen")
    out_dir <- str_glue("{dir}/{id}_phenos")
    dir_create(out_dir)
    out_pheno <- str_glue("{out_dir}/{id}_p{p}_h2{h2}_seed{seed}_hseed{haptools_seed}.pheno")
    print(out_pheno)
    cmd <- str_glue("{a$haptools_path} simphenotype {pgen} {out_causal2} -o {out_pheno} --seed {haptools_seed} --heritability {h2} --chunk-size 100000")
    system(cmd)
    cmd <- str_glue("sed '1c #IID\tpheno' {out_pheno} > /tmp/sim_phenos; mv /tmp/sim_phenos {out_pheno}")
    system(cmd)
}

for (id in a$ids) {
    for (i in 1:nrow(g)) {
        p_i <- g$p[i]
        h2_i <- g$h2[i]
        cat(str_glue("id={id}, p={p_i}, h2={h2_i}"), "\n", sep = "")
        sim_phenos(dir = a$dir, p = p_i, h2 = h2_i, id = id, seed = a$seed, haptools_seed = a$haptools_seed)
    }
}
