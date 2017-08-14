setwd("D:/Box Sync/Data/")

library(tidyverse)
library(ggthemes)
theme_set(theme_few(base_size = 16))

filename <- "compiledLabeled.csv"

if (!file.exists("compiledLabeled.csv")){
        url <- "https://raw.githubusercontent.com/JamesHWade/XenograftProteinProfiling/master/compiledLabeled.csv"
        filename <- basename(url)
        download.file(url, filename)
}


# Read in data
dat <- read_csv("compiledLabeled.csv")
# dat$Treatment <- factor(dat$Treatment, levels = c("(+)-Serum", "(-)-Serum", 
#                 "DMSO", "Erlotinib", "GNE-317", "Apitolisib", "Palbociclib"))

# Save current wd to return to later and setwd to plots folder
directory <- getwd()
setwd("../XPP_Plots/")

## Plot all data combined
g <- ggplot(dat, aes(x = Target, y = `Net Shift`))

all.point <- g + 
        geom_point(aes(color = Treatment), position = "jitter", alpha = 0.7) +
        facet_grid(`Cell Line` ~ `Time Point`) + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Full Dataset")
        

ggsave(all.point, filename = "everything_point.png", 
       width = 12, height = 8)

all.boxplot <- g + 
        geom_boxplot(aes(fill = Treatment)) +
        facet_grid(`Cell Line` ~ `Time Point`) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle("Full Dataset")

ggsave(all.boxplot, filename = "everything_boxplot.png",
       width = 12, height = 8)

## Write loop to go through each target and plot the result

targetList <- unique(dat$Target)

for(i in targetList) {
        
        dat.tar <- filter(dat, Target == i & Treatment != "(+)-Serum" & 
                        Treatment != "(-)-Serum")
        
        g <- ggplot(dat.tar, aes(x = Treatment, y = `Net Shift`))
        
        plotName <- unlist(strsplit(i, "/"))[1]
        
        target.point <- g + 
                geom_point(aes(color = Treatment), position = "jitter") + 
                facet_grid(`Cell Line`~`Time Point`) +
                ggtitle(paste0("Target: ", i)) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        target.box <- g + 
                geom_boxplot(aes(fill = Treatment)) + 
                facet_grid(`Cell Line`~`Time Point`) +
                ggtitle(paste0("Target: ", i)) + 
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        dir.create("Target Plots")
        
        ggsave(target.point, 
               filename = paste0("Target Plots/", plotName, "_point.png"), 
               width = 12, height = 8)
        
        ggsave(target.box, 
               filename = paste0("Target Plots/", plotName, "_boxplot.png"),
               width = 12, heigh = 8)
}

## Write loop to go through each treatment and plot the result

treatmentList <- unique(dat$Treatment)

for (i in treatmentList) {
        
        dat.rx <- filter(dat, Treatment == i)
        
        g <- ggplot(dat.rx, aes(Target, `Net Shift`))
        
        fig4 <- g + 
                geom_point(aes(color = Target), position = "jitter") + 
                facet_grid(`Cell Line`~`Time Point`) +
                ggtitle(paste0("Treatment: ", i)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        fig5 <- g + 
                geom_boxplot(aes(fill = Target)) + 
                facet_grid(`Cell Line`~`Time Point`) +
                ggtitle(paste0("Treatment: ", i)) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        dir.create("Treatment Plots")
        
        ggsave(fig4, filename = paste0("Treatment Plots/", i, "_point.png"),
               width = 12, height = 8)
        
        ggsave(fig5, filename = paste0("Treatment Plots/", i, "_boxplot.png"),
               width = 12, heigh = 8)
}

## Plot Treatments

# Function to plot individual treatments
PlotTreatment <- function(control, treatment, targets){
        dat.cntl <- filter(dat, Treatment == control & `Time Point` == 1)
        dat.all <- rbind(filter(dat, Treatment == treatment), dat.cntl)
        
        dat.all$`Cell Line` <- factor(dat.all$`Cell Line`, 
                                      labels = c("GBM 6", "GBM 26"))
        
        dat.all$`Time Point` <- factor(dat.all$`Time Point`)
        
        g.all <- ggplot(dat.all, 
                        aes(x = interaction(`Time Point`, Treatment, Target),
                            y = `Net Shift`, 
                            fill = Target)) +
                geom_boxplot() + 
                facet_grid(`Cell Line`~.) +
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(unique(dat.all$Target)))) +
                xlab("Treatment Time") + 
                ggtitle(paste0("Treatment: ", treatment)) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(g.all, 
               filename = paste0("Treatment Plots/", treatment, ".png"), 
               width = 12, height = 6)
        
        dat.rx.6 <- filter(dat.all, Target %in% targets & 
                                   `Cell Line` == "GBM 6")
        
        g <- ggplot(dat.rx.6, 
                    aes(x = interaction(`Time Point`, Treatment, 
                                        Target, `Cell Line`), 
                        y = `Net Shift`, 
                        fill = Target))
        
        txt_6 <- g +
                geom_boxplot() + 
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(targets))) +
                xlab("Treatment Time") + 
                ggtitle(paste0("Treatment: ", treatment, ", GBM 6")) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(txt_6, 
               filename = paste0("Treatment Plots/", treatment, "_GBM6.png"),
               width = 8, height = 6)
        
        dat.rx.26 <- filter(dat.all, Target %in% targets & 
                                    `Cell Line` == "GBM 26")
        
        g <- ggplot(dat.rx.26, 
                    aes(x = interaction(`Time Point`, Treatment, 
                                        Target, `Cell Line`), 
                        y = `Net Shift`, 
                        fill = Target))
        
        txt_26 <- g + 
                geom_boxplot() + 
                labs(fill = "") + 
                scale_x_discrete(labels = rep(c("0 h", "1 h", "24 h"), 
                                              length(targets))) +
                xlab("Treatment Time") + 
                ggtitle(paste0("Treatment: ", treatment, ", GBM 26")) +
                theme(axis.text.x = 
                              element_text(angle = 90, hjust = 1, vjust = 0.5))
        
        ggsave(txt_26, 
               filename = paste0("Treatment Plots", treatment, "_GBM26.png"), 
               width = 8, height = 6)
}

PlotTreatment(control = "DMSO",
              treatment = "Palbociclib",
              targets = c("pRb Ser780", "pRb Ser807/11", 
                          "pGSK3b Ser9", "pPDK1 Ser341"))
PlotTreatment(control = "DMSO",
              treatment = "Erlotinib",
              targets = c("pMAPK Thr202/Tyr204","pAkt Thr308", "pAkt Ser473",
                           "pS6 Ser235/6", "pS6 Ser240/4", "pmTOR Ser2448"))
PlotTreatment(control = "DMSO",
              treatment = "GNE-317",
              targets = c("pAkt Ser473", "pS6 Ser235/6", "pS6 Ser240/4", 
                           "pmTOR Ser2448", "pGSK3b Ser9", "pp70S6K Thr389"))
PlotTreatment(control = "DMSO",
              treatment = "Apitolisib",
              targets = c("pAkt Ser473", "pS6 Ser235/6", "pS6 Ser240/4",
                          "pmTOR Ser2448", "pGSK3b Ser9", "pp70S6K Thr389"))

## Pairwise Treatment Comparisons

# Target list
compTargets <- c("pAkt Ser473", "pS6 Ser235/6", "pS6 Ser240/4",
                 "pp70S6K Thr389", "pRb Ser780", "pRb Ser807/11")

# Pairwise treatment list
txtPairs <- combn(unique(dat$Treatment), 2, simplify = FALSE)

# Function to plot treatment comparisons
TreatmentComp <- function(treatments, targets, cellLine){
        dat.rx <- filter(dat, Treatment == treatments & Target %in% targets &
                                 `Cell Line` == cellLine)
        
        g <- ggplot(dat.rx, aes(x = Target, 
                                group = interaction(Treatment, 
                                                    `Time Point`,
                                                    Target),
                                y = `Net Shift`))
        
        fig_comp <- g + 
                geom_boxplot(aes(fill = interaction(Treatment, `Time Point`))) + 
                labs(fill = "") + 
                xlab("Target") + 
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = 
                              element_text(angle = 45, hjust = 1))
        
        fig_comp.2 <- g + 
                geom_boxplot(aes(fill = Treatment)) + 
                facet_wrap(~`Time Point`) + 
                labs(fill = "") + 
                xlab("Target") + 
                ggtitle(paste0(treatments[1], " vs ", treatments[2])) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        dir.create(path = "Treatment Comparisons")
        
        ggsave(fig_comp,
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, ".png"), 
               width = 8, height = 6)
        ggsave(fig_comp.2, 
               filename = paste0("Treatment Comparisons/", treatments[1], "_", 
                                 treatments[2], "_", cellLine, "_2.png"), 
               width = 8, height = 6)
}

# Run through pair-wise list to plot treatment comparisons
lapply(txtPairs, function(i){
        TreatmentComp(treatments = as.vector(i), targets = compTargets,
                      cellLine = 6)
        TreatmentComp(treatments = as.vector(i), targets = compTargets,
                      cellLine = 26)
})
