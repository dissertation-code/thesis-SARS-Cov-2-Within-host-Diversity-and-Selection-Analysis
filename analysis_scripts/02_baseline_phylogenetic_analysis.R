################################################################################
### Constructing a baseline phylogenetic tree for LUNAR and Usual Care
# Contents:
# 1) Settings
# 2) Loading metadata and IQ tree file
# 3) Function to add LINEAGE + Persistent data to base tree
# 4) Circular (fan) tree
# 5) Vertical (rectangular) tree [Figure 11]

## Includes figures:
# 1) Figure 11. Baseline phylogenetic tree for a subset of both cohorts (n = 69)



# Load packages 
library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(tibble)
library(ggnewscale)
library(phangorn)

# ----------------------------------------------------------------------------
# 0. Settings 
# ----------------------------------------------------------------------------

tree_path         <- "path/to/your/data/directory/iqt.treefile"   
metadata_path     <- "path/to/your/data/directory/metadata_edited.tsv"
output_circular   <- "Figure_circular.pdf"
output_vertical   <- "Figure_vertical.pdf"

lineage_levels <- c("BA.1.1" ,"BA.1.1.15" ,"BA.2", "BA.2.1", "BA.2.10", "BA.5", "BA.5.1", "BA.5.2", "BA.5.3", 
                    "BA.5.3.1", "XE", "BF.1", "XQ")

lineage_cols <- c(
  "BA.1.1"    = "#08519C",  # BA.1 family - blues
  "BA.1.1.15" = "#6BAED6",
  "BA.2"      = "#00441B",  # BA.2 family - greens
  "BA.2.1"    = "#41AB5D",
  "BA.2.10"   = "#A1D99B",
  "BA.5"      = "#7F2704",  # BA.5 family - oranges
  "BA.5.1"    = "#A63603",
  "BA.5.2"    = "#D94801",
  "BA.5.3"    = "#F16913",
  "BA.5.3.1"  = "#FDAE6B",
  "BF.1"      = "#FDD0A2",  # BA.5.2.1 descendant
  "XE"        = "#6A51A3",  # recombinants - purples
  "XQ"        = "#DD3497"
)

persistent_cols <- c("NO" = "#8FA8BF", "YES" = "firebrick")

# ----------------------------------------------------------------------------
# 1. Load metadata
# ----------------------------------------------------------------------------
metadata <- read.delim(metadata_path, stringsAsFactors = FALSE)

metadata$Persistent <- factor(metadata$Persistent, levels = c("NO", "YES"))

metadata$LINEAGE <- factor(metadata$LINEAGE, levels = lineage_levels)

rownames(metadata) <- metadata$SampleID

lineage_df    <- metadata["LINEAGE"]
persistent_df <- metadata["Persistent"]


# ----------------------------------------------------------------------------
# 2. Load the tree (straight from IQ-TREE, rooted in R)
# ----------------------------------------------------------------------------

tree <- read.tree(tree_path)      

cat("Node labels found:", length(tree$node.label), "\n")
print(head(tree$node.label))


tree <- phangorn::midpoint(tree)

# Bootstrap support -> numeric
if (length(tree$node.label) == tree$Nnode) {
  if (any(grepl("/", tree$node.label))) tree$node.label <- sub(".*/", "", tree$node.label)
  tree$node.label <- suppressWarnings(as.numeric(tree$node.label))
} else {
  tree$node.label <- NULL
}
have_support <- !is.null(tree$node.label) && any(!is.na(tree$node.label))

cat("Tips:", Ntip(tree), " Nodes:", tree$Nnode, " Support:", sum(!is.na(tree$node.label)), "\n")


support_layers <- if (have_support) {
  list(
    geom_nodepoint(aes(colour = as.numeric(label),
                       size   = as.numeric(label)), na.rm = TRUE),
    scale_colour_gradient(low = "#8B8000", high = "darkgreen",
                          name = "Bootstrap\nsupport (%)", limits = c(0, 100)),
    scale_size_continuous(range = c(0.65, 4), limits = c(0, 100),
                          name = "Bootstrap\nsupport (%)")
  )
} else list()

tree_depth <- max(node.depth.edgelength(tree))

# ----------------------------------------------------------------------------
# 3. Reusable function: adds the LINEAGE + Persistent rings to a base tree
# ----------------------------------------------------------------------------
add_annotation_rings <- function(base_tree_plot, ring_width, ring_offset_step) {
  
  p <- base_tree_plot %<+% metadata
  
  p1 <- gheatmap(p, lineage_df,
                 width = ring_width, offset = 0,
                 colnames = FALSE, color = NA) +
    scale_fill_manual(values = lineage_cols, name = "Lineage",
                      na.value = "grey90", drop = FALSE) +
    new_scale_fill()
  
  p2 <- gheatmap(p1, persistent_df,
                 width = ring_width, offset = ring_offset_step,
                 colnames = FALSE, color = NA) +
    scale_fill_manual(values = persistent_cols, name = "Persistent",
                      na.value = "grey95", drop = FALSE) +
    theme(legend.position   = "right",
          legend.title      = element_text(size = 10, face = "bold"),
          legend.text       = element_text(size = 8),
          legend.key.size   = unit(0.4, "cm"))
  
  p2
}

# ----------------------------------------------------------------------------
# 4. Circular (fan) version
# ----------------------------------------------------------------------------
p_circular_base <- ggtree(td, size = 0.6, layout = "fan", open.angle = 10) +
  geom_treescale(x = 0, y = -1) +
  support_layers

ring_frac <- 0.05   # ring width as fraction of plot width

p_circular <- add_annotation_rings(
  p_circular_base,
  ring_width       = ring_frac,
  ring_offset_step = ring_frac * tree_depth + tree_depth * 0.02
)

ggsave(plot = p_circular, filename = output_circular, height = 10, width = 12, dpi = 300)
cat("Saved", output_circular, "\n")

# ----------------------------------------------------------------------------
# 5. Vertical (rectangular) version 
# Figure 11. Baseline phylogenetic tree for a subset of both cohorts 
# ----------------------------------------------------------------------------
p_vertical_base <- ggtree(tree, size = 0.6, layout = "rectangular") +
  geom_tiplab(size = 3, align = FALSE, offset = tree_depth * 0.004) +  
  geom_treescale(x = 0, y = -2) +
  support_layers

bar_frac   <- 0.02                   # bar width as fraction of tree width
bar_units  <- bar_frac * tree_depth  # same bar in x-axis units
lab_space  <- tree_depth * 0.25      # room for tip labels before bar 1
gap        <- tree_depth * 0.01      # gap between the two bars

add_rings_rect <- function(base_plot) {
  p <- base_plot %<+% metadata
  p1 <- gheatmap(p, lineage_df, width = bar_frac, offset = lab_space,
                 colnames = FALSE, color = NA) +
    scale_fill_manual(values = lineage_cols, name = "Lineage",
                      na.value = "grey90", drop = FALSE) +
    new_scale_fill()
  gheatmap(p1, persistent_df, width = bar_frac,
           offset = lab_space + bar_units + gap,
           colnames = FALSE, color = NA) +
    scale_fill_manual(values = persistent_cols, name = "Persistent",
                      na.value = "grey95", drop = FALSE) +
    theme(legend.position    = "right",
          legend.title       = element_text(size = 10, face = "bold"),
          legend.text        = element_text(size = 8),
          legend.key.size    = unit(0.4, "cm"),
          legend.box.spacing = unit(2, "pt"),
          legend.margin      = margin(0, 0, 0, 0),
          legend.box.margin  = margin(0, 0, 0, -5),
          plot.margin        = margin(5, 1, 5, 5))
}

p_vertical <- add_rings_rect(p_vertical_base) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.02)))

# Save output as a pdf
ggsave(plot = p_vertical, filename = output_vertical,
       height = 12, width = 9, dpi = 300)

# Save output as a png 
ggsave(plot = p_vertical, "Vertical Tree_final.png", width = 9, height = 12, dpi = 300)










