# haptools_simu

```bash
git clone https://github.com/krillinor/haptools_simu.git
cd haptools_simu
```

Simulate genotypes from 1KG (possibly admixed genomes) and phenotypes.

NB: `haptools` is really high memory usage when >=50K samples... f.x., preallocates up to 200GB...

## Requirements

1. [haptools](https://haptools.readthedocs.io)
2. `R` $\geq$ 4.1 and `Rscript` in path
3.  `R` libraries
	- `bigsnpr` (`remotes::install_github("privefl/bigsnpr"`)
	- `data.table`
	- `docopt`
	- `fs`
	- `stringr`
4. `plink2`
	```bash
	wget https://s3.amazonaws.com/plink2-assets/alpha5/plink2_linux_x86_64_20240526.zip
	unzip plink2_linux_x86_64_20240526.zip
	rm plink2_linux_x86_64_20240526.zip
	```

## Example simulation

Use `--resource_dir` in the `R` scripts below to point to the resource directory. Otherwise, follow ([these steps](#setup-download-resources)) to download them.

### 1. Create a model file

For example, 10 samples, 100% CEU (EUR), 1 generation.

```bash
ID="ceu"
n_samples=10
echo -e "${n_samples} Admixed CEU YRI\n1 0 1 0" > model_${ID}.dat
```

### 2. Create a region map

Assign variables for the next steps (must modify based on current environment).

```bash
dir=${PWD}
ids="ceu"  # corresponding to the model_${ID}.dat, possibly a list, f.x., "ceu,yri"
resource_dir="/home/kristjan/work/haptools_simu"
haptools_path="/home/kristjan/.local/bin/haptools"
```

`haptools` uses a lot of memory when simulating from whole chromosomes, so we partition each chromosome by f.x. 25MB.

```bash
Rscript ${dir}/src/create_region_map.R --dir=${dir} --size 25e6 --resource_dir=${resource_dir}
```

### 3. Simulate genotypes

Uses `haptools simgenotype` to simulate the genotypes ([link](https://haptools.readthedocs.io/en/stable/commands/simgenotype.html)).

Default (not specified):
- `--haptools_seed=159`

```bash
for chrom in {1..22}; do
	Rscript ${dir}/src/sim_genos.R --dir=${dir} --ids=${ids} --chrom=${chrom} --haptools_path=${haptools_path} --resource_dir=${resource_dir}
done
```

### 4. Merge genotypes into a single file

Uses `haptools simphenotype` ([link](https://haptools.readthedocs.io/en/stable/commands/simphenotype.html)) to simulate phenotypes from the causal effects created with `./src/create_betas.R`. `haptools simphenotypes` adds the environmental effect and scales according to the specified heritability.

```bash
region_map="${dir}/region_map"
plink2="${dir}/plink2"

for id in $(echo ${ids} | tr "," " "); do
	out_pmerge="${dir}/pmerge_${id}"
	out_pgen="${dir}/${id}_genos/${id}_hseed${haptools_seed}"
	awk -v OFS='' -v id=${id} -v haptools_seed=${haptools_seed} '{print id,"_genos/",id,"_",$0,"_hseed",haptools_seed}' ${region_map} > ${out_pmerge}
	${plink2} --pmerge-list ${out_pmerge} --make-pfile --out ${out_pgen}
	rm ${out_pgen}-merge.*
	${plink2} --pfile ${id} --out ${id} --freq
done
```

Remove the intermediate files (but keep the `*.bp` files)

```bash
for id in $(echo ${ids} | tr "," " "); do
	rm ${dir}/${id}_genos/${id}_*:*.pgen
	rm ${dir}/${id}_genos/${id}_*:*.pvar
	rm ${dir}/${id}_genos/${id}_*:*.psam
done
```

### 5. Create causal effects (betas)

Currently simulates effects like this:

$$
\begin{equation}
\beta_{j}=\left \{\begin{array}{ll}
\textrm{N} \Bigl( 0, \frac{h_g^2}{Mp_{\textrm{causal}}} \Bigr), \textrm{with probability } p_{\textrm{causal}},\\
0, \textrm{with probability } 1 - p_{\textrm{causal}}
\end{array}
\right.
\end{equation}
$$

(Might be interesting to add local haplotype effects / cross-ancestry genetic correlations / larger effects for rarer variants...)

Default (not specified):
- `--p="0.001,0.01,0.1,0.5"`
- `--h2="0.05,0.1,0.25,0.5,0.8"`
- `--seed=5` (`R` seed)

```bash
Rscript ${dir}/src/create_betas.R --dir=${dir} --resource_dir=${resource_dir}
```

### 6. Simulate phenotypes

Default (not specified):
- `--p="0.001,0.01,0.1,0.5"`
- `--h2="0.05,0.1,0.25,0.5,0.8"`
- `--haptools_seed=159`
- `--seed=5` (`R` seed)

```bash
Rscript ${dir}/src/sim_phenos.R --dir=${dir} --ids=${ids} --haptools_path=${haptools_path}
```

## Setup (download resources)

It's not necessary to run the following setup steps if they exist somewhere.

Use `--resource_dir` in the `R` scripts in the example to point to the directory with the resources.

### Get map files

```bash
wget https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip plink.GRCh38.map.zip
rm plink.GRCh38.map.zip
mkdir map
mv plink.* map
```

### Get 1KG pgen

```bash
wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst
wget https://www.dropbox.com/s/ngbo2xm5ojw9koy/all_hg38_noannot.pvar.zst
wget https://www.dropbox.com/s/2e87z6nc4qexjjm/hg38_corrected.psam
./plink2 --zst-decompress all_hg38.pgen.zst > all_hg38.pgen
./plink2 --zst-decompress all_hg38_noannot.pvar.zst > all_hg38.pvar
mv hg38_corrected.psam all_hg38.psam
rm *.zst
```

### Get 1KG sampleinfo

```bash
wget https://raw.githubusercontent.com/CAST-genomics/haptools/main/example-files/1000genomes_sampleinfo.tsv
```

### Get HM3 SNPs and liftover to hg38

```bash
wget https://www.dropbox.com/s/rhi806sstvppzzz/snpinfo_mult_1kg_hm3
    
# linux
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver

Rscript ./src/hm3_to_hg38.R --dir=${PWD}
```
### Filter 1KG

```bash
Rscript ./src/get_1kg_snplist.R --dir=${PWD}

./plink2 --pfile all_hg38 --out all_hg38_qc --make-pfile --allow-extra-chr --extract hm3_use
```

## Possible TODOs

- Use local haplotypes (and admixed people) (in `./src/sim_genos.R`)
- Simulate causal effects from MVN with cross-ancestry genetic correlations and/or weight by variant freq. (in `./src/create_betas.R`)
