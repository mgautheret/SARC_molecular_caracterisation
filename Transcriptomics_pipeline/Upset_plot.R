# This script characterizes the distribution and co-occurrence of MIBC molecular classes and histological subtypes across patients, using pie charts to summarize their distribution and UpSet plots to visualize subtype/class combinationsobserved within individual patients.


# Load packages
library(dplyr)
library(readr)
library(tibble)
library(ggplot2)
library(tidyr)
library(UpSetR)

# Load metadata
metadata <- read.csv2("~/Stage_M2/proj_SARC/results/01_RNA_seq/0_preprocessing/metadata_clean.csv",sep = ",")



# Select relevant variables and remove duplicates
# Keep only the variables required for the molecular and  histological classification analyses.
metadata_small <- metadata %>%
  select(
    patient,
    MIBC_classif,
    Histological_subtype
  )

# Remove duplicated patient/classification combinations.
# This ensures that each patient is counted once for a given
# molecular class and histological subtype.
metadata_small_unique <- metadata_small %>%
  distinct()


# ============================================================
# 1. MIBC class distribution
# ============================================================

# Display the number of unique patient/classification combinations for each MIBC class
table(metadata_small_unique$MIBC_classif)

# Count the number of patients associated with each MIBC class and calculate their percentage among all observations
molClass_counts <- metadata_small_unique %>%
  count(MIBC_classif) %>%
  mutate(
    percentage = n / sum(n) * 100
  )

# Create a pie chart showing the proportion of MIBC molecular classes.
pie <- ggplot(
  molClass_counts,
  aes(x = "", y = n, fill = MIBC_classif)
) +
  geom_bar(
    stat = "identity",
    width = 1
  ) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(
    aes(
      label = paste0(
        round(percentage, 1),
        "%"
      )
    ),
    position = position_stack(vjust = 0.5)
  ) +
  
  labs(
    fill = "Molecular Class",
    title = "Distribution of Molecular Classes"
  )

print(pie)

# Save the pie chart 
ggsave(
  "~/Stage_M2/proj_SARC/results/01_RNA_seq/1_mol_class/MIBC_pie_chart.pdf",
  pie
)


# ============================================================
# MIBC class co-occurrence
# ============================================================


# For each patient, retrieve all distinct MIBC classes observed across their samples.
# Classes are concatenated into a comma-separated character string to represent the combination of classes observed for each patient
patient_classes <- metadata %>%
  group_by(patient) %>%
  summarize(
    classes = paste(
      sort(unique(MIBC_classif)),
      collapse = ","
    )
  ) %>%
  ungroup()

# Count how frequently each combination of molecular classes
# occurs across patients.
combination_counts <- patient_classes %>%
  count(classes) %>%
  arrange(desc(n))



# create one row per patient/class combination.
link_data <- patient_classes %>%
  separate_rows(
    classes,
    sep = ","
  ) %>%
  group_by(patient) %>%
  mutate(
    order = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    patient_index = as.numeric(factor(patient))
  )



# data are converted to wide format:
#   rows    = patients
#   columns = MIBC classes
#   values  = 1 if the patient has the class, otherwise 0
patient_classes <- metadata %>%
  select(
    patient,
    MIBC_classif
  ) %>%
  distinct() %>%
  mutate(
    value = 1
  ) %>%
  pivot_wider(
    names_from = MIBC_classif,
    values_from = value,
    values_fill = 0
  )

# Convert the result to a data.frame
patient_classes_df <- as.data.frame(patient_classes)

# Assign color to each MIBC class
class_colors <- c(
  "Ba/Sq"       = "#E41A1C",
  "LumU"        = "#377EB8",
  "Stroma-rich" = "#FFFF33",
  "NE-like"     = "#984EA3",
  "LumP"        = "#4DAF4A"
)


# Generate MIBC molecular class UpSet plot
upset(
  patient_classes_df,
  sets = colnames(patient_classes_df)[-1],
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = "black",
  sets.bar.color = class_colors[
    colnames(patient_classes_df)[-1]
  ],
  point.size = 3.5,
  line.size = 1,
  text.scale = c(
    1.3, 1.3, 1, 1, 1.3, 1.2
  ),
  mainbar.y.label = "Number of class co-occurrences",
  sets.x.label = "Consensus class"
)


# ============================================================
# Histological subtype distribution
# ============================================================


# Count the number of unique observations for each histological subtype and calculate their percentage.
Hist_sub_counts <- metadata_small_unique %>%
  count(Histological_subtype) %>%
  mutate(
    percentage = n / sum(n) * 100
  )

# Create a pie chart representing the proportion of histological subtypes.
ggplot(
  Hist_sub_counts,
  aes(
    x = "",
    y = n,
    fill = Histological_subtype
  )
) +
  geom_bar(
    stat = "identity",
    width = 1
  ) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(
    aes(
      label = paste0(
        round(percentage, 1),
        "%"
      )
    ),
    position = position_stack(vjust = 0.5)
  ) +
  
  labs(
    fill = "Histological_subtype",
    title = "Distribution of Histological subtype"
  )


# ============================================================
# Histological subtype co-occurrence
# ============================================================


# For each patient, retrieve all distinct histological subtypes observed across their samples.
patient_classes <- metadata %>%
  group_by(patient) %>%
  summarize(
    classes = paste(
      sort(unique(Histological_subtype)),
      collapse = ","
    )
  ) %>%
  ungroup()


# Count the frequency of each combination of histological subtypes across patients.
combination_counts <- patient_classes %>%
  count(classes) %>%
  arrange(desc(n))



# Split the comma-separated subtype lists into individual rows to obtain one row per patient/subtype combination.
link_data <- patient_classes %>%
  separate_rows(
    classes,
    sep = ","
  ) %>%
  group_by(patient) %>%
  mutate(
    order = row_number()
  ) %>%
  ungroup() %>%
  mutate(
    patient_index = as.numeric(factor(patient))
  )



# Create a binary patient-by-histological-subtype matrix.
#   1 = subtype observed in the patient
#   0 = subtype not observed in the patient
patient_classes <- metadata %>%
  select(
    patient,
    Histological_subtype
  ) %>%
  distinct() %>%
  mutate(
    value = 1
  ) %>%
  pivot_wider(
    names_from = Histological_subtype,
    values_from = value,
    values_fill = 0
  )

# Convert to data frame
patient_classes_df <- as.data.frame(patient_classes)


# Assign color to each histological subtype.
Histological_subtype <- c(
  "sarcomatoid"                 = "#98FB98",
  "conventional"                = "#9B5DE5",
  "giant cells"                 = "#E9C46A",
  "glandular differentiation"   = "tomato",
  "micropapillary"              = "skyblue",
  "neuroendocrine"              = "#264653",
  "poorly differentiated"       = "violetred",
  "squamous differentiation"    = "grey78",
  "clear cell"                  = "slategrey",
  "lipid-rich"                  = "blue",
  "trophoblastic differentiation" = "olivedrab",
  "tubular and microcystic"     = "hotpink",
  "plasmacytoid"                = "red",
  "nested"                      = "#2A9D8F",
  "admixed subtypes"            = "gold1",
  "large nested"                = "chartreuse"
)


# Generate histological subtype UpSet plot

upset(
  patient_classes_df,
  sets = colnames(patient_classes_df)[-1],
  order.by = "freq",
  keep.order = TRUE,
  main.bar.color = "black",
  sets.bar.color = Histological_subtype[
    colnames(patient_classes_df)[-1]
  ],
  
  point.size = 2,
  line.size = 1,
  text.scale = c(
    1.3, 1.3, 1, 1, 1.3, 1.2
  ),
  
  mainbar.y.label = "Number of subtype co-occurrences",
  sets.x.label = "Histological subtype"
)


# ============================================================
# End of script
# ============================================================
