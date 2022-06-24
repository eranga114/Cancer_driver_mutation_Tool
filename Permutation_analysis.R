#install.packages("tidyverse")
library(tidyverse)

# Read in the files, exclude unnecessary columns, rename columns.

ptm_data = read_csv("/Users/srp33/Downloads/human_PTMs_phosphositeplus_03_07.csv")

ptm_data = select(ptm_data, GENE, residue, pro_len) %>%
  rename(Protein = GENE, Position = residue, Length = pro_len)

mutation_data = read_csv("/Users/srp33/Downloads/missense_driver_mutations (2).csv") %>%
  select(gene_name, mutproteinPosStart) %>%
  rename(Protein = gene_name, Position = mutproteinPosStart)

# Exclude PTMs at the start or end of a gene

ptm_data = filter(ptm_data, Position >= 4 & Position <= (Length - 3))

# Remove proteins that have multiple (ambiguous) lengths. I believe these have different ACC_IDs.

duplicate_proteins = select(ptm_data, Protein, Length) %>%
  distinct() %>%
  group_by(Protein) %>%
  summarize(Count = n()) %>%
  filter(Count > 1) %>%
  pull(Protein)

ptm_data = filter(ptm_data, !(Protein %in% duplicate_proteins))

# Limit the data to proteins that have at least one PTM and one mutation.

proteins = intersect(pull(ptm_data, Protein), pull(mutation_data, Protein)) %>%
  sort()

ptm_data = filter(ptm_data, Protein %in% proteins) %>%
  arrange(Protein, Position)

mutation_data = filter(mutation_data, Protein %in% proteins) %>%
  arrange(Protein, Position)

# Identify unique proteins

protein_data = select(ptm_data, Protein, Length) %>%
  distinct()

# Build a tibble with 7-amino-acid windows around PTMs.

positions = pull(ptm_data, Position)
lengths = pull(ptm_data, Length)
ptm_data_expanded = tibble(Protein = rep(pull(ptm_data, Protein), 6),
                           Position = c(positions - 3, positions - 2, positions - 1,
                                        positions + 1, positions + 2, positions + 3),
                           Length = rep(lengths, 6)) %>%
  bind_rows(ptm_data) %>%
  arrange(Protein, Position)

# Calculate proportion of mutations that occurred within 7 amino-acid windows.

ptm_mutation_data = inner_join(ptm_data_expanded, mutation_data) %>%
  distinct()

proportion_mutations_in_ptm_window = nrow(ptm_mutation_data) / nrow(mutation_data)

# Perform permutation analysis

permutation_out_file_path = "Permutation_Proportions.tsv"

if (!file.exists(permutation_out_file_path)) {
  set.seed(0)
  num_permutations = 10000
  
  protein_summary = group_by(mutation_data, Protein) %>%
    summarize(Mutation_Count = n()) %>%
    inner_join(protein_data)
  
  all_proportion_mutations_in_ptm_window = tibble(iteration = 0, proportion_mutations_in_ptm_window)
  
  generate_random_mutations = function(mutation_count, length) {
    return(Position = sample(1:length)[1:mutation_count])
  }
  
  for (iteration in 1:num_permutations) {
    print(iteration)
    
    permuted_ptm_mutation_data = group_by(protein_summary, Protein) %>%
      summarize(Position = generate_random_mutations(Mutation_Count, Length)) %>%
    inner_join(ptm_data_expanded) %>%
    distinct() %>%
    suppressMessages()
    
    permuted_proportion_mutations_in_ptm_window = nrow(permuted_ptm_mutation_data) / nrow(mutation_data)
    
    all_proportion_mutations_in_ptm_window = rbind(all_proportion_mutations_in_ptm_window, c(iteration, permuted_proportion_mutations_in_ptm_window))
  }
  
  write_tsv(all_proportion_mutations_in_ptm_window, "Permutation_Proportions.tsv")
}