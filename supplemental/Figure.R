library(data.table)
library(optparse)
library(dplyr)
library(stringr)
library(nnls)
library(ggplot2)
library(cowplot)
library("ggbeeswarm")
library("grid")

# Figure 2
fontsize = 6
plot.simu <- function(tag, tagtext, sizelist, eff="Constant") {

dat0 <- fread(paste0(tag, ".fig2.tsv"))
colnames(dat0)[3] <- "test"
mult <- dat0[dat0$method == "Multi.pop" & 
             (dat0$valid==dat0$test | (dat0$valid=="AMR" & dat0$test=="tbd")| (dat0$valid=="tbd" & dat0$test=="AMR")), ]
disco.mult <- dat0[dat0$method == "Disco.multi", ]
dat <- merge(mult, disco.mult, by=c("i", "test", "Nc"))
dat$increase <- (dat$R.sqr.y - dat$R.sqr.x)/ dat$R.sqr.x
dat <- dat[dat$test == "tbd", ]

print(tag)
dat$sig <- "NS"
dat$p <- 1
for (TEST in c("tbd")){
for (NC in c(100, 300, 1000, 10000)){
old <- dat[dat$test == TEST & dat$Nc == NC, ]$R.sqr.x
new <- dat[dat$test == TEST & dat$Nc == NC, ]$R.sqr.y
p = t.test(old, new, paired = TRUE)$p.value
print(paste(TEST, NC, p, mean(dat[dat$test == TEST & dat$Nc == NC, ]$increase)))
dat[dat$test == TEST & dat$Nc == NC, ]$sig <- ifelse(p<0.05 & mean(old) < mean(new), "+*", dat[dat$test == TEST & dat$Nc == NC, ]$sig)
dat[dat$test == TEST & dat$Nc == NC, ]$sig <- ifelse(p<0.05 & mean(old) > mean(new), "-*", dat[dat$test == TEST & dat$Nc == NC, ]$sig)
dat[dat$test == TEST & dat$Nc == NC, ]$sig <- ifelse(p<0.01 & mean(old) < mean(new), "+**", dat[dat$test == TEST & dat$Nc == NC, ]$sig)
dat[dat$test == TEST & dat$Nc == NC, ]$sig <- ifelse(p<0.01 & mean(old) > mean(new), "-**", dat[dat$test == TEST & dat$Nc == NC, ]$sig)
dat[dat$test == TEST & dat$Nc == NC, ]$sig <- ifelse(p<0.001 & mean(old) < mean(new), "+***", dat[dat$test == TEST & dat$Nc == NC, ]$sig)
dat[dat$test == TEST & dat$Nc == NC, ]$sig <- ifelse(p<0.001 & mean(old) > mean(new), "-***", dat[dat$test == TEST & dat$Nc == NC, ]$sig)
}}
dat1 <- dat %>%
group_by(Nc, test) %>%
summarise(
    y_pos = max(increase, na.rm = TRUE) + 0.001,
    sig = unique(sig),
    .groups = "drop"
  ) %>%
ungroup()
dat1$Nc <- as.factor(dat1$Nc)
print(min(dat$increase))
print(max(dat$increase))
cat("\n\n\n") 

dat$Nc <- as.factor(dat$Nc)
plot.a <- ggplot(dat, aes(x= Nc, y=increase)) + 
  geom_violin(linewidth=0.15, position=position_dodge(), width=0.5, scale="width") +
  geom_point(size=1, alpha=0.5) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.95, size=0.4) +
  scale_color_manual(values = c("yes.pos" = "blue", "yes.neg" = "darkred", "yes2.pos" = "cyan", "yes2.neg" = "red", "no"="darkgrey"), guide="none") +
  xlab("Number of causal SNPs") + 
  ylab(bquote('Relative '~R^2~'increase')) +
  scale_y_continuous(labels = scales::percent, limits = c(-0.020, 0.08), expand =c(0, 0), breaks = seq(-0.020, 0.08, by = 0.01)) + 
  geom_abline(linewidth=0.5, alpha=0.3, slope =0, intercept=0) +
  geom_text(
    data=dat1, 
    aes(x = Nc, y = y_pos, label = sig),
    inherit.aes = FALSE, 
    vjust = 0, 
    size = 3 ) +
  theme_classic() +
  theme(text = element_text(size=fontsize), 
        title = element_text(size=fontsize),
        axis.title = element_text(size=fontsize + 1),
        axis.text = element_text(size=fontsize + 1),
        # legend.text = element_blank(),
        # legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust=1, hjust = 1, size=fontsize),
        strip.background = element_rect( fill="lightgrey", linetype = "blank"),
        panel.spacing = unit(1, "lines")
  ) +
  labs(tag = tagtext) +
  annotate("text", x = 0.55, y=0.077, label= paste("AFR:", sizelist[1] ), hjust = 0, size =fontsize/.pt) + 
  annotate("text", x = 0.55, y=0.074, label= paste("EAS:", sizelist[2] ), hjust = 0, size =fontsize/.pt) + 
  annotate("text", x = 0.55, y=0.071, label= paste("EUR:", sizelist[3] ), hjust = 0, size =fontsize/.pt) + 
  annotate("text", x = 0.55, y=0.068, label= paste("SAS:", sizelist[4] ), hjust = 0, size =fontsize/.pt) + 
  annotate("text", x = 0.55, y=0.065, label= paste("Causal SNP effect size:", eff ), hjust = 0, size = fontsize/.pt) 
   
return(plot.a)}


plot.a <- plot.simu("comb1", "A", c("100k", "100k", "100k", "100k"))
plot.b <- plot.simu("comb2", "B", c("40k", "40k", "100k", "40k"))
plot.c <- plot.simu("comb3", "C", c("0", "40k", "100k", "0"))
plot.d <- plot.simu("comb4.ctn", "D", c("100k", "100k", "100k","100k"), "Varying") + 
annotate("text", x = 4.2, y=-0.006, label= "NS: non-significant", hjust = 1, size =fontsize/.pt) +
annotate("text", x = 4.2, y=-0.010, label= "*: P<0.05" , hjust = 1, size =fontsize/.pt) + 
annotate("text", x = 4.2, y=-0.014, label= "**: P<0.01" , hjust = 1, size =fontsize/.pt) +
annotate("text", x = 4.2, y=-0.018, label= "***: P<0.001", hjust = 1, size =fontsize/.pt) +
annotate("text", x = 2.1, y=-0.014, label= "+: increase", hjust = 1, size =fontsize/.pt) + 
annotate("text", x = 2.1, y=-0.018, label= "-: decrease" , hjust = 1, size =fontsize/.pt)

plot <- plot_grid(plot.a, plot.b, plot.c, plot.d, nrow=1, rel_widths = c(4, 4, 4, 4))
show(plot)
pdf(paste0("Fig2.pdf"),  width=8, height=4.5)
print(plot)
dev.off()

# Figure 3
comp <- data.frame(PRS = character(), pop=character(), PHENO=character(), 
BETA=double(), SE=double(), P=double(), Adj.R2=double(), N=integer())
for (tag in c("BMI", "SBP", "DBP", "hdladj", "ldladj", "choladj", "logtg")) {
dat <- fread(paste0("../data/", tag,".2.result.tsv"))
comp <- rbind(comp, dat)
}

mix <- comp[PRS=="DD", ]
match <- comp[PRS==pop | (PRS=="tbd" & pop=="AMR"), ]
other <- comp[PRS!="DD", ]
max <- comp[PRS!="DD" & PRS!="tbd" , ] %>%
    group_by(PHENO, pop) %>%
    filter(Adj.R2==max(Adj.R2)) %>%
    ungroup


comp1 <- merge(mix, match, by =c("pop", "PHENO"))
comp1$increase <- (comp1$Adj.R2.x - comp1$Adj.R2.y)/ comp1$Adj.R2.y

comp1$sig <- "NS"
for (POP in (unique(comp1$pop))){
old <- comp1[pop==POP, ]$Adj.R2.y
new <- comp1[pop==POP, ]$Adj.R2.x
p = t.test(old, new, paired = TRUE)$p.value       
print(paste(POP, p, mean(comp1[pop==POP, ]$increase)))
comp1[pop==POP, ]$sig <- ifelse(p<0.05 & mean(old) < mean(new), "+*", comp1[pop==POP, ]$sig)
comp1[pop==POP, ]$sig <- ifelse(p<0.05 & mean(old) > mean(new), "-*", comp1[pop==POP, ]$sig)
comp1[pop==POP, ]$sig <- ifelse(p<0.01 & mean(old) < mean(new), "+**", comp1[pop==POP, ]$sig)
comp1[pop==POP, ]$sig <- ifelse(p<0.01 & mean(old) > mean(new), "-**", comp1[pop==POP, ]$sig)
comp1[pop==POP, ]$sig <- ifelse(p<0.001 & mean(old) < mean(new), "+***", comp1[pop==POP, ]$sig)
comp1[pop==POP, ]$sig <- ifelse(p<0.001 & mean(old) > mean(new), "-***", comp1[pop==POP, ]$sig)
}

fontsize = 10

comp1$pop <- ifelse(comp1$pop=="tbd", "OTH", comp1$pop)
comp1$pop <- factor(comp1$pop, level = c("AFR", "EAS",  "SAS", "EUR", "AMR", "OTH"))
phenolabel = c("BMI"="BMI", "SBP"="SBP", "DBP"="DBP", "choladj"="TC", "hdladj"="HDL", "ldladj"="LDL", "logtg"="lg(TG)")
colorvalue = c("BMI"="#a6cee3",  "DBP"="#1f78b4", "SBP"="#b2df8a", "choladj"="#33a02c", "hdladj"="#fb9a99", "ldladj"="#e31a1c", "logtg"="#fdbf6f")

comp1$increase2 <- ifelse(comp1$increase > 0.8, comp1$increase-0.4, comp1$increase)

comp2 <- comp1 %>%
group_by(pop) %>%
summarise(
    y_pos = max(increase2, na.rm = TRUE) +0.01,
    sig = unique(sig),
    .groups = "drop"
  ) %>%
ungroup()

xpos=5.54
ypos=0.37
ystep=0.031
fontsize2=9

plot1 <- ggplot(comp1, aes(x=pop, y=increase2)) + 
    geom_violin(linewidth=0.15, position=position_dodge(), width=0.3, scale="width") +
    geom_abline(linewidth=0.5, alpha=0.3, slope =0, intercept=0) +
    geom_beeswarm(aes(color = PHENO), cex=2, size=2, 
                priority="ascending", dodge.width = 0, groupOnX = T)+
    stat_summary(fun = mean, geom = "crossbar", width = 0.7, size=0.3) +
    expand_limits(x = 0, y=0.02) +
    geom_text(data=comp2, aes(x = pop, y = y_pos, label = sig),
              inherit.aes = FALSE, vjust = 0, size = 4 ) +
    theme_classic() +
    theme(text = element_text(size=fontsize), 
          title = element_text(size=fontsize),
          legend.title=element_blank(),
          legend.position=c(.87,.9)) +
    scale_color_manual(values = colorvalue, labels = phenolabel) +
    scale_y_continuous(labels = scales::percent) + scale_y_continuous(
    breaks = c(-0.1, 0, 0.1, 0.2, 0.3, 0.5, 0.55, 0.6),  # show "before and after" values
    labels = c("-10%", "0%", "10%","20%", "30%", "90%", "95%", "100%")
  ) +
    xlab("Target population") + 
    ylab(bquote('Relative '~R^2~'increase')) +
    annotate("text", x = xpos, y=ypos, label= "+    increase", hjust = 0, size =fontsize2/.pt) + 
    annotate("text", x = xpos, y=ypos-ystep, label= "-    decrease" , hjust = 0, size =fontsize2/.pt) +
    annotate("text", x = xpos, y=ypos-2*ystep, label= "NS non-significant", hjust = 0, size =fontsize2/.pt) +
    annotate("text", x = xpos, y=ypos-3*ystep, label= "*    P<0.05" , hjust = 0, size =fontsize2/.pt) + 
    annotate("text", x = xpos, y=ypos-4*ystep, label= "**   P<0.01" , hjust = 0, size =fontsize2/.pt) 
# annotate("text", x = xpos, y=ypos-5*ystep, label= "***: P<0.001", hjust = 1, size =fontsize/.pt) +
show(plot1)

pdf(paste0("/medpop/esp2/yruan/projects/disco.plot/update2.fig3.pdf"),  width=8, height=6)
print(plot1)
adj = 0.005
grid.lines(y = unit(c(0.756, 0.766), "npc"),x = unit(c(0.075 + adj, 0.075+ adj), "npc"), gp = gpar(lwd =2, col = "white"))
grid.lines(y = unit(c(0.760, 0.775), "npc"),x = unit(c(0.065+ adj, 0.085+ adj), "npc"), gp = gpar(lwd =1.5))
grid.lines(y = unit(c(0.748, 0.764), "npc"),x = unit(c(0.065+ adj, 0.085+ adj), "npc"), gp = gpar(lwd =1.5))

dev.off()

# Figure 4
comp <- fread("fig4.tsv")
plot.fig4 <- function(tag, comp){
comp0 <- comp[PHENO == tag, ]
print(dim(comp0))
plot.a<-ggplot(comp0, aes(x= PRS, y=BETA, color=PRS)) + 
  geom_point(size=1, alpha=0.7)+
  geom_errorbar(aes(ymax=up, ymin=low), width=0.5) +
  ylab("BETA")+
  # scale_x_discrete(labels= phenolabel) + 
  # scale_color_discrete(labels= phenolabel, guide="none") +
  # facet_grid(PHENO ~ pop, labeller = as_labeller(phenolabel))+
  facet_grid(cohort ~ pop)+
  theme_classic() +
  theme(text = element_text(size=fontsize), 
        title = element_text(size=fontsize),
        axis.ticks.x = element_line(),
        legend.position="none",
        strip.background = element_rect( fill="lightgrey", linewidth=1.5, linetype = "blank"),
        axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, size=fontsize),
        # axis.text.x = element_blank(),
        panel.spacing = unit(0.5, "lines")
  ) +
  geom_hline(yintercept = -Inf, linewidth = 0.6)
  ggtitle(tag)
return(plot.a)
}
plot.cad <- plot.fig4("CAD", comp)
plot.cad <- plot.fig4("DM2", comp)
plot <- plot_grid(plot.cad, plot.dm2, nrow=2)
show(plot)
pdf(paste0("fig4.pdf"),  width=8, height=8)
print(plot)
dev.off()

