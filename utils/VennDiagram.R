# =============================================================================
# Venn Diagram for Cross Cohort Analysis
# =============================================================================
# This script creates Venn diagrams to visualize the overlap between different
# cohorts in terms of brain regions and baseline variables.
# 
# Required packages:
# - ggvenn: For creating Venn diagrams
# =============================================================================

# Load required library
library(ggvenn)

# =============================================================================
# PART 1: Brain Regions Analysis
# =============================================================================

# Brain regions for UCT study
uct <- c("Ventricles", "Supratentorial Tissue", "Thalamus", "Supratentorial CSF", "Corpus Callosum", "Putamen")

# Brain regions for PRISMA study
prisma <- c("Corpus Callosum", "Ventricles", "Supratentorial CSF", "Thalamus", "Putamen", "Supratentorial Tissue", "Cerebellum")

# Brain regions for MINE study
mine <- c("Putamen", "Ventricles", "Corpus Callosum", "Cerebellum", "Thalamus", "Supratentorial Tissue", "Supratentorial CSF")

# Brain regions for PRIMES study
primes <- c("Corpus Callosum", "Putamen", "Ventricles", "Caudate", "Supratentorial CSF", "Thalamus")

# Combine all brain region datasets into a list for Venn diagram
brain_sets <- list(
  UCT = uct,
  MINE = mine,
  PRISMA = prisma, 
  PRIMES = primes
)

# Create Venn diagram for brain regions
ggvenn(
  brain_sets,
  show_elements = TRUE,        # Show individual elements (brain region names) in each section
  label_sep = "\n",            # Break each region name onto a new line for better readability
  fill_color = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#CAB2D6"),  # Color scheme for each set
  stroke_size = 0.6,           # Thickness of the circle borders
  set_name_size = 4,           # Size of the set labels (UCT, MINE, etc.)
  text_size = 3                # Size of the element text (brain region names)
)

# Save the brain regions Venn diagram
ggsave("img/venn_brain_regions.png", width = 7, height = 7)

# =============================================================================
# PART 2: Baseline Variables Analysis
# =============================================================================
# Define baseline variables included in each study/dataset:
# - WAZ: Weight-for-age Z-score
# - Birth Length: Length at birth
# - Age at Scan: Age when MRI scan was performed
# - Gestation Weeks: Gestational age at birth
# - Maternal Age: Age of the mother

# Baseline variables for UCT study
uct <- c("WAZ")

# Baseline variables for PRISMA study
prisma <- c("Birth Length", "Age at Scan", "Gestation Weeks")

# Baseline variables for MINE study
mine <- c("Birth Length", "Gestation Weeks", "WAZ", "Maternal Age")

# Baseline variables for PRIMES study
primes <- c("WAZ", "Age at Scan")

# Combine all baseline variable datasets into a list for Venn diagram
base_sets <- list(
  UCT = uct,
  MINE = mine,
  PRISMA = prisma,
  PRIMES = primes
)

# Create Venn diagram for baseline variables
ggvenn(
  base_sets,
  show_elements = TRUE,        # Show individual elements (variable names) in each section
  label_sep = "\n",            # Break each variable name onto a new line for better readability
  fill_color = c("#A6CEE3", "#B2DF8A", "#FB9A99", "#CAB2D6"),  # Same color scheme as brain regions
  stroke_size = 0.6,           # Thickness of the circle borders
  set_name_size = 4,           # Size of the set labels (UCT, MINE, etc.)
  text_size = 3                # Size of the element text (variable names)
)

# Save the baseline variables Venn diagram
ggsave("img/venn_baseline.png", width = 7, height = 7)