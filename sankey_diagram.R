# This script creates a sankey diagram from metadata table of CaveSARC + VESPER RNA seq data
# It compares the MIBC class of SARC zones and non SARC zones from the same patient 

# Load packages
library(networkD3)
library(dplyr)
library(tidyr)
library(htmlwidgets)
library(webshot)


# Load metadata
metadata <- read.csv2(
  "~/Stage_M2/proj_SARC/results/01_RNA_seq/0_preprocessing/metadata_clean.csv",
  sep = ","
)

metadata <- as.data.frame(metadata)


links <- metadata %>%
  dplyr::select(
    Sample.ID,
    patient,
    MIBC_classif,
    Sarcomatoid
  )


# Reshape the data so that each patient has separate columns for the Non Sarcomatoid and Sarcomatoid samples
df_pairs <- links %>%
  dplyr::select(
    patient,
    MIBC_classif,
    Sarcomatoid
  ) %>%
  pivot_wider(
    names_from = Sarcomatoid,
    values_from = MIBC_classif,
    values_fn = list
  ) %>%
  unnest(`Non Sarcomatoid`) %>%
  unnest(Sarcomatoid)


# Create a unique label representing the changes between the MIBC class of the Non Sarcomatoid component and the Sarcomatoid component
df_pairs$pairs <- paste0(
  df_pairs$`Non Sarcomatoid`,
  "_",
  df_pairs$Sarcomatoid
)


# Count how many patient/sample pairs show each MIBC class change 
df_count<- df_pairs %>%
  group_by(pairs) %>%
  count()



# Define Sankey nodes
# "_NS" = Non Sarcomatoid
# "_S"  = Sarcomatoid
# The "group" variable is used to assign the same color to the same MIBC subtype on both sides of the Sankey 
nodes <- data.frame(
  label = c(
    "Ba/Sq_NS",
    "LumU_NS",
    "LumP_NS",
    "NE-like_NS",
    "Stroma-rich_NS",
    "Ba/Sq_S",
    "NE-like_S",
    "Stroma-rich_S"
  ),
  
  group = c(
    "BaSq",
    "LumU",
    "LumP",
    "NElike",
    "Stroma",
    "BaSq",
    "NElike",
    "Stroma"
  )
)



# Build Sankey links
# Split each subtype transition into a source and target.
# Source = MIBC subtype of the Non Sarcomatoid component
# Target = MIBC subtype of the Sarcomatoid component
# The transition count is stored in the "value" column
links <- df_count %>%
  separate(
    pairs,
    into = c("source", "target"),
    sep = "_"
  ) %>%
  mutate(value = n) %>%
  dplyr::select(-n)


# networkD3 identifies nodes using zero-based numeric indices. Therefore, source subtype names are mapped to the corresponding node positions in nodes df 
links$source[links$source == "Ba/Sq"] <- 0
links$source[links$source == "LumU"] <- 1
links$source[links$source == "LumP"] <- 2
links$source[links$source == "NE-like"] <- 3
links$source[links$source == "Stroma-rich"] <- 4


# Convert target subtype names to Sankey node indices

# Sarcomatoid nodes occupy positions 5-7 in nodes2.
links$target[links$target == "Ba/Sq"] <- 5
links$target[links$target == "NE-like"] <- 6
links$target[links$target == "Stroma-rich"] <- 7


# Convert source and target indices to numeric values
links <- as.data.frame(links)
links$source <- as.numeric(links$source)
links$target <- as.numeric(links$target)


# Combine node and link information into a single object used to generate the Sankey
plot_list<- list(
  nodes = nodes,
  links = links
)



# Define color scale
colors <- 'd3.scaleOrdinal()
  .domain(["BaSq", "LumU", "LumP", "NElike", "Stroma"])
  .range(["#E41A1C", "#377EB8", "green", "#984EA3", "#FFFF33"])'



# Generate Sankey diagram

sankey <- sankeyNetwork(
  Links = plot_list$links,
  Nodes = plot_list$nodes,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "label",
  NodeGroup = "group",
  fontSize = 0,
  nodeWidth = 30,
  colourScale = colors
)

sankey

# Export Sankey
setwd("~/Stage_M2/proj_SARC/results/01_RNA_seq/1_mol_class")

# Save interactive version
saveWidget(
  p,
  "sankey.html",
  selfcontained = TRUE
)

# Convert html to png 
webshot(
  "sankey.html",
  file = "sankey.png",
  vwidth = 400,
  vheight = 600,
  zoom = 3
)

# ============================================================
# End of script
# ============================================================

