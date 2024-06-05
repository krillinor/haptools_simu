library(docopt)
library(data.table)
library(stringr)

"Create a region map.

Usage:
  create_region_map.R --dir=<dir> [--size=<size> --resource_dir=<resource_dir>]
  create_region_map.R (-h | --help)

Options:
  --dir=<dir>                    Working directory.
  --size=<size>                  Region size [default: 25e6].
  --resource_dir=<resource_dir>  Directory with resources, f.x. 1KG pgen files. Assumes working dir. if not specified [default: NULL].
  -h --help                      Show this screen.

" -> doc

a <- docopt(doc)
a$size <- as.numeric(a$size)

if (a$resource_dir == "NULL") {
    a$resource_dir <- a$dir
}

make_map <- function(dir, chr, size) {
    d <- fread(str_glue("{a$resource_dir}/map/plink.chr{chr}.GRCh38.map"))
    s <- seq(1, max(d$V4), by = size)
    s <- c(s, max(d$V4))
    ss <- c()
    for (i in 1:(length(s) - 1)) {
        s_i <- as.integer(s[i])
        s_ip1 <- as.integer(s[i + 1] - 1)
        tmp <- str_glue("{chr}:{s_i}-{s_ip1}")
        ss <- c(ss, tmp)
    }
    ss
}

m <- c()
for (i in 1:22) {
    tmp <- make_map(a$dir, i, a$size)
    m <- c(m, tmp)
}

out <- str_glue("{a$dir}/region_map")
cat(str_glue("Writing to {out}"), "", "\n")
fwrite(list(m), out)
