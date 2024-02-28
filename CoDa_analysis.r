library("tidyverse")
library("DESeq2")

# Remove low count genes, as it could just be noise:
count_matrix <- read_delim("abun.csv") |> 
  `rownames<-`(c("1001","1011","8020","8024")) |> 
  t()

count_matrix <- count_matrix[rowSums(count_matrix) > 10, ]

meta_data <- data.frame(country = as.factor(c("DK","DK","NL","NL"))) |> 
  `rownames<-`(c("1001","1011","8020","8024"))

gene_length <- read_delim(file = "ResFinder_gene_lengths.csv", col_names = F) |> 
  `colnames<-`(c("AMR_gene", "count")) |> 
  mutate(AMR_gene = as.factor(AMR_gene))

gene_length <- gene_length |> 
  filter(AMR_gene %in% rownames(count_matrix))

gene_length_longer <- gene_length |> pivot_wider(names_from = AMR_gene, values_from = count) |> 
  select(mcols(dds)@rownames) |> rownames_to_column(var = "rownames") |> 
  pivot_longer(cols = -c(rownames))


all(rownames(meta_data) == colnames(count_matrix))

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = meta_data,
                              design = ~ country)

mcols(dds)$basepairs = gene_length_longer$value

normalized_count_matrix <- fpkm(dds)
ALR <- log2(normalized_count_matrix+1) |> as.data.frame()

gmean <- function(row) {
  exp((1 / length(row)) * sum(log(row)))
}

data_CLR <- t(normalized_count_matrix)
CLR <- apply(data_CLR+1, 1, function(x) log2(x / gmean(x))) |> as.data.frame()

ALR |> rownames_to_column(var = "AMR_GENE") |>pivot_longer(cols = -c("AMR_GENE")) |> 
  ggplot(aes(x = reorder(AMR_GENE, -value), y = value)) + 
  geom_bar(aes(fill = name),stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  labs(title = "Total sum of ALR for AMR genes ",
       x = "AMR gene")

data_CLR |> t() |> 
  data.frame() |> 
  rownames_to_column(var = "AMR_GENE") |> 
  pivot_longer(cols = -c("AMR_GENE")) |> 
  ggplot(aes(x = reorder(AMR_GENE, -value), y = value)) + 
  geom_bar(aes(fill = name),stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  labs(title = "Total sum of CLR for AMR genes ", 
       x = "AMR gene")

normalized_count_matrix |> 
  as.data.frame() |> 
  rownames_to_column(var = "AMR_GENE") |>
  pivot_longer(cols = -c("AMR_GENE")) |> 
  ggplot(aes(x = reorder(AMR_GENE, -value), y = value)) + 
  geom_bar(aes(fill = name),stat = "identity") + 
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) + 
  labs(title = "Total sum of ALR for AMR genes ",
       x = "AMR gene")

