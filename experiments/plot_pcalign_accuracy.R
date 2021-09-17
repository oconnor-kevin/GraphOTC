library(tidyverse)
library(stringr)
library(latex2exp)

homedir <- # redacted for anonymity
options(stringsAsFactors=FALSE)

## Read Data.
accuracy_df <- read.csv(file.path(homedir, "alignment_accuracies.txt"))
edge_accuracy_df <- read.csv(file.path(homedir, "edge_alignment_accuracies.txt"))

## Wrangle data
acc_df <- accuracy_df %>%
  group_by(Algorithm, Lambda, Sigma) %>%
  summarise("mean_acc"=mean(Accuracy),
            "max_acc"=max(Accuracy),
            "min_acc"=min(Accuracy))
acc_df$Algorithm[acc_df$Algorithm == "OTC"] <- "GraphOTC"
acc_df$Sigma <- factor(acc_df$Sigma)
levels(acc_df$Sigma) <- c("High Overlap", "Moderate Overlap", "Low Overlap")

edge_acc_df <- edge_accuracy_df %>%
  group_by(Algorithm, Lambda, Sigma) %>%
  summarise("mean_acc"=mean(Accuracy),
            "max_acc"=max(Accuracy),
            "min_acc"=min(Accuracy))
edge_acc_df$Algorithm[edge_acc_df$Algorithm == "OTC"] <- "GraphOTC"
edge_acc_df$Sigma <- factor(edge_acc_df$Sigma)
levels(edge_acc_df$Sigma) <- c("High Overlap", "Moderate Overlap", "Low Overlap")

## Make plots
plot_id <- format(Sys.time(), "%m-%d-%y_%H-%M-%S")
### Vertex accuracy
ggplot(acc_df, 
       aes(x=log10(Lambda), y=mean_acc, ymin=min_acc, ymax=max_acc, group=interaction(factor(Algorithm), factor(Sigma)), color=factor(Algorithm))) + 
  facet_wrap(~Sigma) + 
  geom_line(size=1) + 
  #geom_point(size=2) + 
  geom_ribbon(alpha = 0.1, size=0.05) +
  theme_minimal() +
  theme(legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.title.align=0.5,
        legend.box.background = element_rect(colour="black"),
        plot.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=8),
        axis.title = element_text(size=12, face="bold"))+
  ggtitle("Vertex Alignment Accuracy") + 
  ylim(0, 1) +
  labs(group = TeX('$\\log_{10}\\lambda$'),
       color = "Algorithm",
       x = TeX('$\\log_{10}\\lambda$'),
       y = "Accuracy")
ggsave(file.path(homedir, paste0("vertex_accuracy_", plot_id, ".png")), width=8, height=3)

### Edge accuracy
ggplot(edge_acc_df, 
       aes(x=log10(Lambda), y=mean_acc, ymin=min_acc, ymax=max_acc, group=interaction(factor(Algorithm), factor(Sigma)), color=factor(Algorithm))) + 
  facet_wrap(~Sigma) + 
  geom_line(size=1) + 
  geom_ribbon(alpha = 0.1, size=0.05) +
  theme_minimal() +
  theme(legend.title = element_text(size=12), 
        legend.text = element_text(size=10),
        legend.title.align=0.5,
        legend.box.background = element_rect(colour="black"),
        plot.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=8),
        axis.title = element_text(size=12, face="bold"))+
  ggtitle("Edge Alignment Accuracy") + 
  ylim(0, ceiling(max(edge_acc_df$max_acc*10))/10) +
  labs(group = TeX('$\\log_{10}\\lambda$'),
       color = "Algorithm",
       x = TeX('$\\log_{10}\\lambda$'),
       y = "Accuracy")
  ggsave(file.path(homedir, paste0("edge_accuracy_", plot_id, ".png")), width=8, height=3)

