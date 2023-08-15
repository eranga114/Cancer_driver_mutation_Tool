library(tidyverse)
library(ggrepel)

# Read in data from CSV file
tnk2_exp <- read_csv("stad_coad.csv")

# Convert data to tidy format and calculate log2 expression values
tnk2_exp_tidy <- tnk2_exp %>%
  pivot_longer(cols = c("STAD_No_trunc", "STAD_trunc", "COAD_No_trunc", "COAD_trunc"), 
               names_to = "tnk2_exp", 
               values_to = "expression") %>%
  mutate(log2_expression = log2(expression),
         group = ifelse(tnk2_exp %in% c("STAD_No_trunc", "COAD_No_trunc"), "Group1", "Group2"))

ggplot(tnk2_exp_tidy, aes(x = tnk2_exp, y = log2_expression, color = group)) +
  geom_boxplot() +
  geom_point() +
  geom_text_repel(aes(label = label), 
                  color = "black",
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = "black",
                  size = 3) +
  stat_summary(fun.data = function(x) data.frame(y = max(x), label = paste0("n=", length(x))), 
               geom = "text", vjust = -2.5, color = "black", size = 4) +
  labs(x = "", y = "Log2 Expression") +
  theme_classic() +
  scale_x_discrete(labels = c("STAD_No_trunc" = "STAD\nno truncations",
                              "STAD_trunc" = "STAD\ntruncated",
                              "COAD_No_trunc" = "COAD\nno truncations",
                              "COAD_trunc" = "COAD\ntruncated")) +
  scale_color_manual(values = c("Group1" = "#69b3a2", "Group2" = "grey")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 14, color = "black", face = "bold", margin = margin(0, 0, 20, 0)),
        axis.title = element_text(size = 16, face = "bold", margin = margin(20, 0, 0, 0)),
        plot.title = element_text(size = 20, face = "bold"))

ggsave("tnk2_exp_plot.pdf", width = 8, height = 9.5, dpi = 300)
