#set up for ggplot
require(ggplot2)
require(devtools)
#devtools::install_github('bbc/bbplot')
require(bbplot)
require(grid)
require(gridExtra)
require(DescTools)
require(ggpubr)
require(ggstance)
require(ggprism)
require(patchwork)
require(scales)
require(data.table)
require(pals)
require(RColorBrewer)

# DIAMETER ACROSS TRYPTOPHAN CONDITIONS ##### 
# import data & process
setwd("~/Desktop/EleftheriaGitHub")
data <- read.csv("FigureS1_data.csv", stringsAsFactors = F)
data <- subset(data, is.na(Diameter) == F)

data$mean <- NA
data$sd <- NA
for(i in seq(1, length(data$Strain))){
  data$mean[i] <- mean(data$Diameter[data$Tryptophan.concentration == data$Tryptophan.concentration[i] & data$Strain == data$Strain[i]])
  data$sd[i] <- sqrt(var(data$Diameter[data$Tryptophan.concentration == data$Tryptophan.concentration[i] & data$Strain == data$Strain[i]]))
}

# stats ##### 
model <- aov(Diameter ~ as.factor(Tryptophan.concentration) * Strain, data=data)
summary(model)
stats <- TukeyHSD(model)
stats <- as.data.frame(stats$`as.factor(Tryptophan.concentration):Strain`)

stats$strain1 <- gsub( "(.*)-(.*)", "\\1", rownames(stats))
stats$tryp1 <- gsub( "(.*):(.*)", "\\1", stats$strain1)
stats$strain1 <- gsub( "(.*):(.*)", "\\2", stats$strain1)
stats$strain2 <- gsub( "(.*)-(.*)", "\\2", rownames(stats))
stats$tryp2 <- gsub( "(.*):(.*)", "\\1", stats$strain2)
stats$strain2 <- gsub( "(.*):(.*)", "\\2", stats$strain2)

stats <- subset(stats, tryp1 == tryp2)

stats$stars <- "ns"
for(i in seq(1, length(stats$diff))){
  if(stats$`p adj`[i] < 0.05){ stats$stars[i] <- '*' }
  if(stats$`p adj`[i] < 0.01){ stats$stars[i] <- '**' }
  if(stats$`p adj`[i] < 0.001){ stats$stars[i] <- '***' }
  if(stats$`p adj`[i] < 0.0001){ stats$stars[i] <- '****' }
}

# final plot ##### 
ggplot(data=data, aes(x=as.factor(Tryptophan.concentration), y=Diameter)) +
  geom_point(aes(x=as.factor(Tryptophan.concentration), y=Diameter, colour=factor(Strain, levels=c('PK',"ΔtrpC","pTrpC"))),
             position = position_jitterdodge(jitter.width=.1, dodge.width=.6), size=2) +
  geom_boxplot(aes(x=as.factor(Tryptophan.concentration), y=Diameter, fill=factor(Strain, levels=c('PK',"ΔtrpC","pTrpC"))), 
               position = position_dodge(width=.6), width=.5, alpha=.25) +
  geom_errorbar(aes(x=as.factor(Tryptophan.concentration), ymin=mean-sd, ymax=mean+sd, group=factor(Strain, levels=c('PK',"ΔtrpC","pTrpC"))),
                position = position_dodge(width=.6), width=.25, alpha=.5) +
  scale_y_continuous(name='Diameter (cm)', limits=c(0,3.5)) +
  scale_x_discrete(name='Tryptophan concentration (mM)') +
  scale_color_manual(values=c("#F8766D","#00BA38","#619CFF"), breaks=c('PK',"ΔtrpC","pTrpC"),
                     labels=c('PK307-4',expression(Delta*italic("trpC")),expression('' %up% italic("trpC")))) +
  scale_fill_manual(values=c("#F8766D","#00BA38","#619CFF"), breaks=c('PK',"ΔtrpC","pTrpC"),
                    labels=c('PK307-4',expression(Delta*italic("trpC")),expression('' %up% italic("trpC")))) +
  #stats tryptophan 0
  geom_segment(aes(x = .8, y = 3.3, xend = 1.2, yend = 3.3)) +
  annotate('text',x=1, y=3.5, label=stats$stars[stats$tryp1 == 0 & stats$strain1 == 'pTrpC' & stats$strain2 == 'PK'], size=5) +
  geom_segment(aes(x = .8, y = 3, xend = 1, yend = 3)) +
  annotate('text',x=.9, y=3.1, label=stats$stars[stats$tryp1 == 0 & stats$strain1 == 'ΔtrpC' & stats$strain2 == 'PK'], size=5) +
  #stats tryptophan 0.1
  geom_segment(aes(x = 1.8, y = 2.7, xend = 2.2, yend = 2.7)) +
  annotate('text',x=2, y=2.9, label=stats$stars[stats$tryp1 == 0.1 & stats$strain1 == 'pTrpC' & stats$strain2 == 'PK'], size=5) +
  geom_segment(aes(x = 1.8, y = 2.4, xend = 2, yend = 2.4)) +
  annotate('text',x=1.9, y=2.5, label=stats$stars[stats$tryp1 == 0.1 & stats$strain1 == 'ΔtrpC' & stats$strain2 == 'PK'], size=5) +
  #stats tryptophan 24.5
  geom_segment(aes(x = 2.8, y = 1.3, xend = 3.2, yend = 1.3)) +
  annotate('text',x=3, y=1.5, label=stats$stars[stats$tryp1 == 0.1 & stats$strain1 == 'pTrpC' & stats$strain2 == 'PK'], size=5) +
  geom_segment(aes(x = 2.8, y = 1, xend = 3, yend = 1)) +
  annotate('text',x=2.9, y=1.1, label=stats$stars[stats$tryp1 == 0.1 & stats$strain1 == 'ΔtrpC' & stats$strain2 == 'PK'], size=5) +
  theme_classic() +
  theme(plot.margin = unit(c(.1,.1,.1,.1),'cm'),
        plot.title = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14, margin = margin(10,0,0,0)),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=14),
        axis.text.x=element_text(size=14),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.minor.y = element_blank(),
        legend.position='top',
        legend.title = element_blank(),
        legend.text = element_text(size=14,margin=margin(r=1,unit='cm')))
