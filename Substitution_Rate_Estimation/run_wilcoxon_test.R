# run_wilcoxon_test.R
# Perform Wilcoxon rank sum test on dN and dS values

# Load necessary packages
library(dplyr)

# Load the data (update this with the actual file paths)
data <- read.csv("dn_ds_data.csv")

# Perform Wilcoxon rank sum test between gene types
wilcox_test <- wilcox.test(data$dN ~ data$GeneType)
p_value <- wilcox_test$p.value

# Adjust P values using false discovery rate (FDR) correction
adjusted_p <- p.adjust(p_value, method = "fdr")

# Output results
cat("P value:", p_value, "\n")
cat("Adjusted P value (FDR):", adjusted_p, "\n")

# Save the results
write.csv(adjusted_p, file = "adjusted_p_values.csv")

# NOTE: Ensure R and necessary packages (e.g., dplyr) are installed.
