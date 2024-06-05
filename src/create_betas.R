library(docopt)
library(data.table)
library(stringr)
library(fs)

"Create betas (causal effects).

Usage:
  create_betas.R --dir=<dir> [--p=<p> --h2=<h2> --seed=<seed> --resource_dir=<resource_dir>]
  create_betas.R (-h | --help)

Options:
  --dir=<dir>                    Working directory.
  --p=<p>                        Polygenicity [default: 0.001,0.01,0.1,0.5].
  --h2=<h2>                      Heritability [default: 0.05,0.1,0.25,0.5,0.8].
  --seed=<seed>                  R seed [default: 159].
  --resource_dir=<resource_dir>  Directory with resources, f.x. 1KG pgen files. Assumes working dir. if not specified [default: NULL].
  -h --help                      Show this screen.

" -> doc

a <- docopt(doc)
a$p <- str_split(a$p, ",")[[1]]
a$h2 <- str_split(a$h2, ",")[[1]]
a$seed <- as.numeric(a$seed)

g <- expand.grid(p = a$p, h2 = a$h2)

if (a$resource_dir == "NULL") {
    a$resource_dir <- a$dir
}
pvar <- fread(cmd = str_glue("grep -v \"##\" {a$resource_dir}/all_hg38_qc.pvar"))

create_betas <- function(dir = dir, vs = vs, p = p, h2 = h2, seed = 159) {
    set.seed(seed)

    m <- length(vs)
    p <- p |>
        as.character() |>
        as.numeric()
    h2 <- h2 |>
        as.character() |>
        as.numeric()

    is_causal <- rbinom(m, 1, p)
    not_causal <- which(is_causal == 0)
    betas_sd <- h2 / (m * p)
    betas <- rnorm(m, mean = 0, sd = betas_sd)
    betas[not_causal] <- 0
    betas <- data.table(snp = vs, beta = betas)

    out_dir <- str_glue("{dir}/betas")
    dir_create(out_dir)
    out_causal <- str_glue("{out_dir}/causal_p{p}_h2{h2}_seed{seed}.snplist")
    betas_filtered <- betas[beta != 0, ]
    fwrite(betas_filtered, out_causal, sep = "\t", col.names = FALSE)
}

for (i in 1:nrow(g)) {
    p_i <- g$p[i]
    h2_i <- g$h2[i]
    cat(str_glue("p={p_i}, h2={h2_i}"), "\n", sep = "")
    create_betas(dir = a$dir, vs = pvar$ID, p = p_i, h2 = h2_i, seed = a$seed)
}
