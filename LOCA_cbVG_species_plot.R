#!/usr/bin/env Rscript

library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(tidyr)
library(RColorBrewer)
library(tibble)
library(reshape2)

args <- commandArgs(trailingOnly = TRUE)
data <- read_delim(args[1], delim=" ", escape_double=FALSE, trim_ws=TRUE, show_col_types = FALSE)
GENOME_LENGTH = args[2]
SAMPLE = args[3]
TYPE = args[4]

# 1. generate new report with species ID based on the mode B_R (new column 'mode')
data_mode_list <- lapply(split(data, f = data$SPECIES),
                         `[.function` <- `[.closure` <- function(df) {
                           uniqv <- unique(df$BREAK_POSITION)
                           mode <- uniqv[which.max(tabulate(match(df$BREAK_POSITION, uniqv)))]
                           rejoin <- df[df$BREAK_POSITION == mode, 4][1,]
                           mode <- paste(mode, rejoin, sep="_")
			   rule <- df$LENGTH%%as.integer(RULE)
			   cbind(df, mode, rule)
                           })

loca_data <- do.call(rbind, data_mode_list)

mode_file <- str_replace(args[1], ".txt$", paste0("_mode_rule",RULE,".txt"))
write_tsv(loca_data, file = mode_file)

# 2. plot species (based on mode)
all_reads = data.frame(table(loca_data$mode), stringsAsFactors = F)
loca_data$N = sapply(loca_data$mode, function(w){all_reads$Freq[all_reads$Var1 == w]})
loca_data = loca_data[!duplicated(loca_data$mode),]
data <- select(loca_data,c("BREAK_POSITION","REJOIN_POSITION","N"))
data$color <- ifelse(TYPE == "CB", '#b6429f', '#4293b6')
data$color[data$N < 3 ] <- 'lightgray'

tiff(str_replace(args[1], ".txt$", "_mode_plot.tiff"), units = 'in', width = 5, height = 4, res = 300)
ggplot(data, aes(x = BREAK_POSITION, y = REJOIN_POSITION))+
  geom_point(aes(size = N), color = data$color, stroke=0, alpha=0.5)+
  scale_size_continuous(limits=c(1,max(data$N)),
                        breaks=unique(round(seq(1,max(data$N),length.out=5))))+
  theme_linedraw()+
  theme(panel.border = element_rect(colour = "black", linewidth=1),
        axis.ticks = element_line(linewidth = 0.5),
        axis.text = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"),
        legend.text = element_text(face="bold"),
        legend.title = element_text(face="bold"),
        plot.title = element_text(face="bold"),
        text = element_text(family = "Arial"))+
  ylab('Rejoin position')+
  xlab('Break position')+
  xlim(c(0,as.numeric(GENOME_LENGTH)))+
  ylim(c(0,as.numeric(GENOME_LENGTH)))+
  labs(size = 'junction reads')+
  ggtitle(paste(SAMPLE,TYPE))+
  geom_abline(linetype="dashed", alpha=0.2)
dev.off()



