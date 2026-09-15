# --- PIPELINE SCRIPT 10: SURVIVAL PLOTTING (FIGURE 8) ---
library(survival)
library(survminer)
library(readr)

# 1. PATH CONFIGURATION (Dynamic relative paths)
script_dir <- dirname(sys.frame(1)$ofile)
if (length(script_dir) == 0) { script_dir <- getwd() }
data_path <- file.path(script_dir, "..", "data", "Figure_8_Survival_Data.csv")
fig_dir <- file.path(script_dir, "..", "figures", "Fig 8")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# 2. LOAD DATA & MATH
print("Loading Python exported survival data...")
df <- read_csv(data_path)
fit <- survfit(Surv(Time, Event) ~ Group, data = df)

# 3. GENERATE PLOT
print("Rendering Kaplan-Meier curve and Risk Tables...")
splots <- ggsurvplot(
  fit, data = df, size = 1.5,
  palette = c("#E63946", "#1D3557"), 
  pval = TRUE, pval.coord = c(5, 0.1),
  risk.table = TRUE, risk.table.col = "strata", risk.table.height = 0.25,
  legend.title = "Biomarker Panel Expression", legend.labs = c("High", "Low"),
  xlab = "Time to Progression (Months)", ylab = "Progression-Free Probability",
  title = "Progression-Free Survival in EGFR-Mutant LUAD",
  ggtheme = theme_classic(base_family = "sans", base_size = 12)
)

# 4. EXPORT TO PDF
output_path <- file.path(fig_dir, "Figure_8_R_Version.pdf")
cairo_pdf(output_path, width = 7, height = 6)
print(splots)
dev.off()

print(paste("--- SCRIPT 10 COMPLETE. Figure saved to:", output_path, "---"))
