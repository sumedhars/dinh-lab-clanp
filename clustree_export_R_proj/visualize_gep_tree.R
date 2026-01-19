
# ============================================================================
# R code to visualize GEP hierarchy
# Run this after loading the exported CSV files
# ============================================================================

# install.packages(c("clustree", "ggraph", "igraph", "dplyr", "ggplot2"))

library(clustree)
library(ggraph)
library(igraph)
library(dplyr)
library(ggplot2)

# --- OPTION 1: Use clustree directly ---
df <- read.csv("clustree_input.csv")
p1 <- clustree(df, prefix = "K")
ggsave("gep_clustree.pdf", p1, width = 12, height = 10)

# --- OPTION 2: Custom ggraph with similarity shown ---
edges <- read.csv("gep_edges_for_r.csv")
nodes <- read.csv("gep_nodes_for_r.csv")

g <- graph_from_data_frame(edges, vertices = nodes, directed = TRUE)

# Custom layout with K levels
layout_df <- nodes %>%
  group_by(K) %>%
  mutate(x = row_number(),
         x = (x - 0.5) / max(1, n()),
         y = -K) %>%
  ungroup()

layout_matrix <- as.matrix(layout_df[match(V(g)$name, layout_df$name), c("x", "y")])

p2 <- ggraph(g, layout = layout_matrix) +
  geom_edge_link(aes(edge_alpha = similarity, edge_color = color),
                 arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                 end_cap = circle(4, 'mm')) +
  geom_node_point(aes(size = size), alpha = 0.9, color = "steelblue") +
  geom_node_text(aes(label = cluster), size = 3.5, color = "white", fontface = "bold") +
  scale_edge_color_identity() +
  scale_edge_alpha(range = c(0.3, 1), name = "Similarity") +
  scale_size(range = c(8, 18), name = "# Pathways") +
  labs(title = "GEP Hierarchical Tree",
       y = "K (number of programs)") +
  theme_graph() +
  theme(axis.text.y = element_text())

ggsave("gep_hierarchy_ggraph.pdf", p2, width = 14, height = 10)

# --- OPTION 3: With edge labels ---
p3 <- ggraph(g, layout = layout_matrix) +
  geom_edge_link(aes(edge_alpha = similarity, edge_color = color,
                     label = sprintf("%.2f", similarity)),
                 arrow = arrow(length = unit(2, 'mm'), type = "closed"),
                 end_cap = circle(4, 'mm'),
                 angle_calc = 'along', label_dodge = unit(2.5, 'mm'), label_size = 2) +
  geom_node_point(aes(size = size), alpha = 0.9, color = "steelblue") +
  geom_node_text(aes(label = cluster), size = 3.5, color = "white", fontface = "bold") +
  scale_edge_color_identity() +
  scale_edge_alpha(range = c(0.3, 1)) +
  scale_size(range = c(8, 18)) +
  theme_graph()

ggsave("gep_hierarchy_labeled.pdf", p3, width = 16, height = 12)
