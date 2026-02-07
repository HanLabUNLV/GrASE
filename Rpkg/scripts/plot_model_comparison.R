library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork) # For combining plots (optional, install if needed)

# --- CONFIGURATION: Update these paths to your result files ---
files <- list(
  Prior    = "/tmp/test_glmmTMB_prior.sort",
  FixedEB  = "/tmp/test_glmmTMB_fixed_EB.sort",
  VGAM     = "/tmp/test_VGAM.sort"
)

# Function to read and clean data
read_results <- function(path, name) {
  if (!file.exists(path)) return(NULL)
  df <- read.table(path, header=FALSE, stringsAsFactors=FALSE)
  
  # Heuristic column naming based on your provided output
  # Adjust column indices if your file format differs slightly
  # Expected: gene, event, LRT, p.value, model, phi, effect_size
  cols <- c("gene", "event", "LRT", "p.value", "model_name", "phi", "effect_size")
  colnames(df)[1:length(cols)] <- cols
  
  df$method <- name
  return(df)
}

# 1. Load Data
all_data <- bind_rows(lapply(names(files), function(n) read_results(files[[n]], n)))

# Handle Binomial Fallback (Phi = NA or Inf)
# We calculate an "effective phi" for plotting.
# If Phi is NA (Binomial), set it to a ceiling value (e.g., 20,000) for visualization
MAX_PHI <- 20000 
all_data <- all_data %>%
  mutate(
    phi_plot = ifelse(is.na(phi) | phi > MAX_PHI, MAX_PHI, phi),
    log_phi = log10(phi_plot),
    log_p   = -log10(p.value + 1e-300) # Avoid Inf
  )

# Reshape for pairwise comparison (Pivot Wider)
wide_df <- all_data %>%
  select(gene, event, method, log_phi, log_p, effect_size) %>%
  pivot_wider(
    names_from = method, 
    values_from = c(log_phi, log_p, effect_size),
    names_sep = "_"
  ) %>%
  filter(!is.na(log_p_Prior)) # Keep only events processed by Prior model (your reference)

# --- PLOT 1: Phi Comparison (Prior vs Fixed) ---
p1 <- ggplot(wide_df, aes(x=log_phi_Prior, y=log_phi_FixedEB)) +
  geom_point(alpha=0.3, size=0.8) +
  geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
  theme_bw() +
  labs(
    title = "Comparison of Dispersion (Phi)",
    subtitle = "Points on Right/Top edge are Binomial (Phi -> Inf)",
    x = "Log10 Phi (Prior Method)",
    y = "Log10 Phi (Fixed EB Method)"
  ) +
  coord_fixed(xlim=c(0, log10(MAX_PHI)), ylim=c(0, log10(MAX_PHI)))

# --- PLOT 2: P-value Comparison (Prior vs Fixed) ---
p2 <- ggplot(wide_df, aes(x=log_p_Prior, y=log_p_FixedEB)) +
  geom_point(alpha=0.3, size=0.8) +
  geom_abline(intercept=0, slope=1, color="red", linetype="dashed") +
  theme_bw() +
  labs(
    title = "Comparison of Significance (-Log10 P)",
    subtitle = "Below line: Fixed EB is Conservative.\nAbove line: Fixed EB is Liberal/Anti-Conservative.",
    x = "-Log10 P-value (Prior Method)",
    y = "-Log10 P-value (Fixed EB Method)"
  )

# --- PLOT 3: Impact of Phi on Significance ---
# Did forcing the Fixed Phi kill significance?
# Compare Phi difference vs P-value difference.
wide_df <- wide_df %>%
  mutate(
    p_diff = log_p_Prior - log_p_FixedEB, # Positive = Prior is more significant
    phi_diff = log_phi_Prior - log_phi_FixedEB # Positive = Prior found less dispersion (cleaner)
  )

p3 <- ggplot(wide_df, aes(x=phi_diff, y=p_diff)) +
  geom_point(alpha=0.3) +
  geom_vline(xintercept=0, color="gray") +
  geom_hline(yintercept=0, color="gray") +
  theme_bw() +
  labs(
    title = "Impact of Dispersion Estimation on Power",
    x = "Log(Phi Prior) - Log(Phi Fixed)\n(Right = Prior found gene was cleaner)",
    y = "Log(P Prior) - Log(P Fixed)\n(Top = Prior found more significance)"
  ) +
  annotate("text", x=2, y=50, label="Prior Correctly\nDetected Clean Gene\n(boosted power)", hjust=0, color="blue")


# Save Plots
pdf("model_comparisons.pdf", width=12, height=4)
# If patchwork is installed:
# p1 + p2 + p3
# Else separate pages:
print(p1)
print(p2)
print(p3)
dev.off()

cat("Plots saved to model_comparisons.pdf\n")
