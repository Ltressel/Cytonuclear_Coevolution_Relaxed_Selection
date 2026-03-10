# generate_boxplots.R
# This script generates boxplots for dN and dS values

# Load necessary libraries
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)  # To combine plots

# Load the substitution rate data from the Excel file
data <- read_excel("Mitochondrial_Substitution_Rates.xlsx")

# 1. dN Boxplot
# Filter data and reshape for dN
data_dN <- data %>%
  filter(Gene_group %in% c("mt", "N-mt OXPHOS", "Glycolysis", "Cell Cycle", "Cytosolic ribosomal")) %>%
  pivot_longer(cols = starts_with("dN_"), names_to = "Species_Group", values_to = "dN_Value") %>%
  filter(Species_Group %in% c("dN_other_Legumes", "dN_50KB")) %>%
  mutate(Species_Group = recode(Species_Group, `dN_other_Legumes` = "Other Legumes", `dN_50KB` = "50-kb Legumes"))

# Perform Wilcoxon test for dN
p_values_dN <- data_dN %>%
  group_by(Gene_group) %>%
  summarise(p_value = wilcox.test(dN_Value ~ Species_Group, data = .)$p.value) %>%
  mutate(p_label = ifelse(p_value < 0.001, "< 0.001", round(p_value, 3)))

# Plot dN boxplot
p_dN <- ggplot(data_dN, aes(x = Gene_group, y = dN_Value, fill = Species_Group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "Nonsynonymous Substitutions (dN)", x = "Gene Group", y = "dN Values") +
  scale_fill_manual(values = c("Other Legumes" = "blue", "50-kb Legumes" = "red")) +
  theme_minimal() +
  geom_text(data = p_values_dN, aes(x = Gene_group, y = -0.01, label = p_label), size = 3, vjust = 1)

# 2. dS Boxplot
# Filter data and reshape for dS
data_dS <- data %>%
  pivot_longer(cols = starts_with("dS_"), names_to = "Species_Group", values_to = "dS_Value") %>%
  filter(Species_Group %in% c("dS_other_Legumes", "dS_50KB")) %>%
  mutate(Species_Group = recode(Species_Group, `dS_other_Legumes` = "Other Legumes", `dS_50KB` = "50-kb Legumes"))

# Perform Wilcoxon test for dS
p_values_dS <- data_dS %>%
  group_by(Gene_group) %>%
  summarise(p_value = wilcox.test(dS_Value ~ Species_Group, data = .)$p.value) %>%
  mutate(p_label = ifelse(p_value < 0.001, "< 0.001", round(p_value, 3)))

# Plot dS boxplot
p_dS <- ggplot(data_dS, aes(x = Gene_group, y = dS_Value, fill = Species_Group)) +
  geom_boxplot(position = position_dodge(0.8)) +
  labs(title = "Synonymous Substitutions (dS)", x = "Gene Group", y = "dS Values") +
  scale_fill_manual(values = c("Other Legumes" = "blue", "50-kb Legumes" = "orange")) +
  theme_minimal() +
  geom_text(data = p_values_dS, aes(x = Gene_group, y = -0.01, label = p_label), size = 3, vjust = 1)

# Combine the dN and dS plots
combined_plot <- (p_dN / p_dS) + plot_layout(ncol = 1)

# Save the combined plot as a PNG
ggsave("Combined_dN_dS_Boxplots.png", combined_plot, width = 10, height = 10, dpi = 300)

# Print the combined plot
print(combined_plot)
