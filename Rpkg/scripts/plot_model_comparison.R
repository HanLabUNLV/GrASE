library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)

# Define paths
f_prior <- "~/DICE/split_exoncnts/test_glmmTMB_MAP_prior.txt"
f_fixed <- "~/DICE/split_exoncnts/test_glmmTMB_fixed_EB.txt"
f_vgam  <- "~/DICE/split_exoncnts/test_VGAM_MLE_EB.txt"

cat("Reading data...\n")
t_fixed <- read.table(f_fixed, header=TRUE, stringsAsFactors=FALSE)
t_prior <- read.table(f_prior, header=TRUE, stringsAsFactors=FALSE)
t_vgam  <- read.table(f_vgam,  header=TRUE, stringsAsFactors=FALSE)

cat("Processing...\n")
# Select and label
df_fixed <- t_fixed %>% 
  select(gene, event, p.value, phi) %>% 
  mutate(method = "FixedEB")

df_prior <- t_prior %>% 
  select(gene, event, p.value, phi) %>% 
  mutate(method = "Prior")

df_vgam <- t_vgam %>% 
  select(gene, event, p.value, phi) %>% 
  mutate(method = "VGAM")

# Combine
all_data <- bind_rows(df_fixed, df_prior, df_vgam)

# Clean and transform
# We limit phi to visualization range since VGAM goes to 10^69
MAX_PHI_PLOT <- 20000 

all_data <- all_data %>%
  mutate(
    p.value = suppressWarnings(as.numeric(p.value)),
    # Handle extremely small p-values (avoid log(0))
    p_clamped = pmax(p.value, 1e-300),
    log_p = -log10(p_clamped),
    
    phi = suppressWarnings(as.numeric(phi)),
    # Clamp phi for plotting visualization
    phi_plot = ifelse(is.na(phi) | phi > MAX_PHI_PLOT, MAX_PHI_PLOT, phi),
    log_phi = log10(phi_plot + 0.1)
  )

# Widen for direct comparison
wide_df <- all_data %>%
  select(gene, event, method, log_phi, log_p) %>%
  pivot_wider(
    names_from = method,
    values_from = c(log_phi, log_p),
    names_sep = "_"
  )

cat("Plotting...\n")

# P1: Phi Comparison (Prior vs FixedEB)
p1 <- ggplot(wide_df, aes(x=log_phi_Prior, y=log_phi_FixedEB)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Dispersion (Log10 Phi): Prior vs Fixed EB", 
       subtitle="Red line = 1:1 identity")

# P2: Phi Comparison (Prior vs VGAM)
p2 <- ggplot(wide_df, aes(x=log_phi_Prior, y=log_phi_VGAM)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Dispersion (Log10 Phi): Prior vs VGAM",
       subtitle="VGAM often estimates infinite (Binomial) phi")

# P3: P-value Comparison (Prior vs VGAM)
p3 <- ggplot(wide_df, aes(x=log_p_Prior, y=log_p_VGAM)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): Prior vs VGAM")

# P4: P-value Comparison (Prior vs Fixed EB)
p4 <- ggplot(wide_df, aes(x=log_p_Prior, y=log_p_FixedEB)) +
  geom_point(alpha=0.2, size=0.5) +
  geom_abline(color="red", linetype="dashed") +
  theme_bw() + 
  labs(title="Significance (-Log10 P): Prior vs Fixed EB")

pdf("model_comparisons.pdf", width=12, height=10)
grid.arrange(p1, p2, p3, p4, nrow=2)
dev.off()

cat("Success. Plots saved to model_comparisons.pdf\n")
