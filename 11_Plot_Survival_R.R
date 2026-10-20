# --- PIPELINE SCRIPT 11: SURVIVAL PLOTTING ---
library(survival)
library(survminer)
library(readr)

# DYNAMIC PATH CONFIGURATION 
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0) { script_dir <- getwd() }
data_dir <- file.path(script_dir, "..", "data")
fig_dir <- file.path(script_dir, "..", "figures")

dir.create(file.path(fig_dir, "Fig 7"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(fig_dir, "Fig 9"), showWarnings = FALSE, recursive = TRUE)

plot_km <- function(data_file, output_file, title_text, xlab_text) {
  df <- read_csv(file.path(data_dir, data_file), show_col_types = FALSE)
  fit <- survfit(Surv(Time, Event) ~ Group, data = df)
  
  splots <- ggsurvplot(
    fit, data = df, size = 1.5, palette = c("#E63946", "#1D3557"), 
    pval = TRUE, pval.coord = c(5, 0.1), risk.table = TRUE, risk.table.col = "strata", 
    risk.table.height = 0.25, legend.title = "Biomarker Panel Expression", 
    legend.labs = c("High", "Low"), xlab = xlab_text, ylab = "Survival Probability", 
    title = title_text, ggtheme = theme_classic(base_family = "sans", base_size = 12)
  )
  
  cairo_pdf(file.path(fig_dir, output_file), width = 7, height = 6)
  print(splots)
  dev.off()
}

plot_km("Figure_7B_Survival_Data.csv", "Fig 7/Figure_7B_TCGA_Survival.pdf", "Overall Survival in EGFR-Mutant LUAD", "Overall Survival (Months)")
plot_km("Figure_9_Survival_Data.csv", "Fig 9/Figure_9_TCGA_EGFR_PFS_Quartiles.pdf", "Progression-Free Survival in EGFR-Mutant LUAD", "Time to Progression (Months)")
print("--- SCRIPT 11 COMPLETE ---")
