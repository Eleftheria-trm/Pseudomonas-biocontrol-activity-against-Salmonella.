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
data <- read.csv("Figure4B_data.csv", stringsAsFactors = F)

data$mean <- NA
data$sd <- NA
for(i in seq(1, length(data$Strain))){
  data$mean[i] <- mean(data$Diameter[data$Tryptophan_concentration == data$Tryptophan_concentration[i] & data$Strain == data$Strain[i]])
  data$sd[i] <- sqrt(var(data$Diameter[data$Tryptophan_concentration == data$Tryptophan_concentration[i] & data$Strain == data$Strain[i]]))
}

# stats ##### 
model <- aov(Diameter ~ as.factor(Tryptophan_concentration) * Strain, data=data)
summary(model)
stats <- TukeyHSD(model)
stats <- as.data.frame(stats$`as.factor(Tryptophan_concentration):Strain`)

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
ggplot(data=data, aes(x=as.factor(Tryptophan_concentration), y=Diameter)) +
  geom_point(aes(x=as.factor(Tryptophan_concentration), y=Diameter, colour=Strain),
             position = position_jitterdodge(jitter.width=.1, dodge.width=.6), size=2) +
  geom_boxplot(aes(x=as.factor(Tryptophan_concentration), y=Diameter, fill=Strain), position = position_dodge(width=.6), width=.5, alpha=.25) +
  geom_errorbar(aes(x=as.factor(Tryptophan_concentration), ymin=mean-sd, ymax=mean+sd, group=Strain),
                position = position_dodge(width=.6), width=.25, alpha=.5) +
  scale_y_continuous(name='Diameter (cm)', limits=c(0,4), expand=c(0,0)) +
  scale_x_discrete(name='Tryptophan concentration (mM)') +
  scale_color_manual(values=c("#00BA38","#619CFF"), breaks=c('PK',"ΔtrpC"), labels=c('PK307-4',expression(Delta*italic("trpC")))) +
  scale_fill_manual(values=c("#00BA38","#619CFF"), breaks=c('PK',"ΔtrpC"), labels=c('PK307-4',expression(Delta*italic("trpC")))) +
  geom_segment(aes(x = .85, y = 3.8, xend = 1.15, yend = 3.8)) +
  annotate('text',x=1, y=3.9, label=stats$stars[stats$tryp1 == 0], size=5) +
  geom_segment(aes(x = 1.85, y = 2.1, xend = 2.15, yend = 2.1)) +
  annotate('text',x=2, y=2.2, label=stats$stars[stats$tryp1 == 0.1], size=5) +
  geom_segment(aes(x = 2.85, y = 1.2, xend = 3.15, yend = 1.2)) +
  annotate('text',x=3, y=1.3, label=stats$stars[stats$tryp1 == 24.5], size=3) +
  theme_classic() +
  theme(plot.margin = unit(c(.5,0,0,.1),'cm'),
        plot.title = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14, margin = margin(10,0,0,0)),
        axis.text.y=element_text(margin = margin(2,2,2,2), size=14),
        axis.text.x=element_text(size=14),
        axis.line = element_line(color = "#cbcbcb"),
        axis.ticks = element_line(color = "#cbcbcb"),
        panel.grid.minor.y = element_blank(),
        legend.position='top',
        legend.title = element_text(size=14),
        legend.text = element_text(size=14))
