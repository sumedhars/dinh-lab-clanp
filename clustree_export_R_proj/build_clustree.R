library(clustree)

df <- read.csv("clustree_input_pathways_by_bestGEP.csv",
               stringsAsFactors = FALSE)

p <- clustree(df, prefix = "K")
print(p)

ggplot2::ggsave("clustree_pathways_by_bestGEP.png",
                p, width = 14, height = 8, dpi = 300)

