library(ggplot2)
library(tidyr)
library(dplyr)

vcf_1 <- read.table("C:/Users/mikcha1/Downloads/vcf_C_pyg_26_reverse_indelRm_intersect_annotated_SNP_is_alt.bed")[,c(1,3,6,20)]
colnames(vcf_1) <- c("chr", "pos", "chCADD", "FORMAT")

# policz warianty
n_variants <- nrow(vcf_1)
vcf_1 <- vcf_1 %>%
  separate(FORMAT, into = "genotype", sep = ":", extra = "drop") %>%
  mutate(genotype = ifelse(genotype == "0/1", "HET", 
                           ifelse(genotype == "1/1", "HOM_ALT", 
                                  ifelse(genotype == "0/0", "HOM_REF", "error"))))

# policz heterozygoty i homozygoty:
table(vcf_1$genotype)
# policz ponownie:
vcf_1 %>%
  group_by(genotype) %>%
  tally()
n_HET <- vcf_1 %>%
  group_by(genotype) %>%
  tally() %>%
  filter(genotype == "HET") %>%
  pull()
n_HOM <- vcf_1 %>%
  group_by(genotype) %>%
  tally() %>%
  filter(genotype == "HOM_ALT" | genotype == "HOM_REF") %>%
  pull()
paste0("W danych zidentyfikowano ", n_HET, " heterozygot i ", n_HOM, " homozygot.")

ggplot(vcf_1, aes(pos, chCADD, color = genotype)) +
  geom_point() +
  scale_color_manual(values = c("HOM_ALT" = "blue", "HOM_REF" = "blue", "HET" = "red")) +
  ggtitle("chCADD, scaffold 1, individual C_pyg_26", 
          subtitle = paste0(n_variants, " variants"))

cor.test(vcf_1$pos, vcf_1$chCADD, method = "spearman")





vcf_2 <- read.table("C:/Users/mikcha1/Downloads/vcf_C_ruf_09_reverse_indelRm_intersect_annotated_SNP_is_alt.bed")[,c(1,3,6,20)]
colnames(vcf_2) <- c("chr", "pos", "chCADD", "FORMAT")
# policz warianty
n_variants <- nrow(vcf_2)
# zmiana 0/0 na hom_ref etc. 
vcf_2 <- vcf_2 %>%
  separate(FORMAT, into = "genotype", sep = ":", extra = "drop") %>%
  mutate(genotype = ifelse(genotype == "0/1", "HET", 
                           ifelse(genotype == "1/1", "HOM_ALT", 
                                  ifelse(genotype == "0/0", "HOM_REF", "error"))))
# policz heterozygoty i homozygoty:
table(vcf_2$genotype)

# policz ponownie:
vcf_2 %>%
  group_by(genotype) %>%
  tally()

n_HET <- vcf_2 %>%
  group_by(genotype) %>%
  tally() %>%
  filter(genotype == "HET") %>%
  pull()

n_HOM <- vcf_2 %>%
  group_by(genotype) %>%
  tally() %>%
  filter(genotype == "HOM_ALT" | genotype == "HOM_REF") %>%
  pull()

paste0("W danych zidentyfikowano ", n_HET, " heterozygot i ", n_HOM, " homozygot.")

ggplot(vcf_2, aes(pos, chCADD, color = genotype)) +
  geom_point() +
  scale_color_manual(values = c("HOM_ALT" = "blue", "HOM_REF" = "blue", "HET" = "red")) +
  ggtitle("chCADD, scaffold 1, individual C_ruf_09", 
          subtitle = paste0(n_variants, " variants"))

cor.test(vcf_2$pos, vcf_2$chCADD, method = "spearman")


vcf_2 %>%
  ggplot(aes(chCADD)) +
  geom_histogram() + 
  facet_grid(row = "genotype")
vcf_2 %>%
  ggplot(aes(log(chCADD) + 1)) +
  geom_histogram() + 
  facet_grid(row = "genotype")
vcf_2 %>%
  ggplot(aes(sqrt(chCADD))) +
  geom_histogram() + 
  facet_grid(row = "genotype")
x <- vcf_2 %>%
  filter(genotype == "HET") %>%
  mutate(chCADD_sqrt = sqrt(chCADD))
shapiro.test(x$chCADD_sqrt) 

y <- vcf_2 %>%
  filter(genotype == "HOM_ALT" | genotype == "HOM_REF") %>%
  mutate(chCADD_sqrt = sqrt(chCADD))
shapiro.test(y$chCADD_sqrt) 

# zaróWno heterozygoty jak i homozygoty nie mają rozkładu normalnego
# Do zbadania róznicy w wartości chCADD wykorzystamy test nieparametryczny: test Manna-Whitneya (Wilcoxona dla dwóch niezależnych grup)

df <- data.frame(wartość = vcf_2$chCADD,
                      grupa = vcf_2$genotype)
wilcox.test(wartość~grupa, data=df)

obs_diff <- mean(df$wartość[df$grupa == "HET"]) - mean(df$wartość[df$grupa == "HOM_ALT"])
paste0("Różnica średniej wartości wynosi ", round(obs_diff, 3))
perm_diffs <- vector("numeric", 1000)
for(i in 1:1000) {
  # Losowo zmieniaj przyporządkowanie wartości do grupy:
  perm_grupa <- sample(df$grupa)
  perm_diffs[i] <- mean(df$wartość[perm_grupa == "HET"]) - mean(df$wartość[perm_grupa == "HOM_ALT"])
}

p_value <- sum(perm_diffs <= obs_diff) / length(perm_diffs)
print(p_value)
# wizualizacja:
perm_diffs_df <- data.frame(perm_diffs = perm_diffs)
ggplot(perm_diffs_df, aes(perm_diffs)) +
  geom_histogram() +
  geom_vline(xintercept = obs_diff, linetype = "dashed")




####### porównane wartości między gatunkami

vcf_1 %>%
  ggplot(aes(chCADD)) +
  geom_histogram()
vcf_1 %>%
  ggplot(aes(log(chCADD) + 1)) +
  geom_histogram()

vcf_1 %>%
  ggplot(aes(sqrt(chCADD))) +
  geom_histogram()

vcf_2 %>%
  ggplot(aes(chCADD)) +
  geom_histogram() 

vcf_2 %>%
  ggplot(aes(log(chCADD) + 1)) +
  geom_histogram()

vcf_2 %>%
  ggplot(aes(sqrt(chCADD))) +
  geom_histogram()

x <- vcf_1 %>%
  mutate(chCADD_sqrt = sqrt(chCADD))
shapiro.test(x$chCADD_sqrt) 

y <- vcf_2 %>%
  mutate(chCADD_sqrt = sqrt(chCADD))
shapiro.test(y$chCADD_sqrt) 

# wartość CADD dla obu gatunków nie mają rozkładu normalnego
# Do zbadania różnicy w wartości chCADD wykorzystamy test nieparametryczny: test Manna-Whitneya (Wilcoxona dla dwóch niezależnych grup)

df <- data.frame(wartość = c(vcf_1$chCADD, vcf_2$chCADD),
                 grupa = c(rep("C_pyg_26", nrow(vcf_1)), rep("C_ruf_09", nrow(vcf_2))))
wilcox.test(wartość~grupa, data=df)




obs_diff <- mean(df$wartość[df$grupa == "C_pyg_26"]) - mean(df$wartość[df$grupa == "C_ruf_09"])
paste0("Różnica średniej wartości wynosi ", round(obs_diff, 3))
perm_diffs <- vector("numeric", 1000)
for(i in 1:1000) {
  # Losowo zmieniaj przyporządkowanie wartości do grupy:
  perm_grupa <- sample(df$grupa)
  perm_diffs[i] <- mean(df$wartość[perm_grupa == "C_pyg_26"]) - mean(df$wartość[perm_grupa == "C_ruf_09"])
}

p_value <- sum(perm_diffs <= obs_diff) / length(perm_diffs)
print(p_value)
# wizualizacja:
perm_diffs_df <- data.frame(perm_diffs = perm_diffs)
ggplot(perm_diffs_df, aes(perm_diffs)) +
  geom_histogram() +
  geom_vline(xintercept = obs_diff, linetype = "dashed")

