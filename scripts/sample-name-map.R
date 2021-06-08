#each element is saved as list element
samples <- snakemake@input[["sample"]]
name <- snakemake@params[["name"]]
outdir <- snakemake@output[[1]]
map_df <- unique(do.call(rbind, Map(data.frame, A=name, B=samples)))

write.table(map_df, outdir, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
