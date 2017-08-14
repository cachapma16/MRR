setwd("D:/Box Sync/Data/")
library(tidyverse)
library(reshape2)
library(Hmisc)
library(corrplot)
library(igraph)


dat <- read_csv("compiledLabeled.csv")
dat <- select(dat, -Shift.1, -Shift.2, - Group)
dat <- filter(dat, !grepl("p53|Abl", Target))
dat$Id <- dat$Ring %% 4

PlotTargets <- function(cellLine, treatment, timePoint){
        
        
        # Filter data based on cell line, treatment, & time point
        dat.tar <- filter(dat, Treatment %in% treatment & `Cell Line` %in% 
                                cellLine & `Time Point` %in% timePoint)
        
        # Recast data as wide format from dataframe
        dat.dcast <<- dcast(dat.tar,
                            Treatment ~ Target + `Cell Line` + `Time Point`,
                            value.var = "Net Shift",
                            fun.aggregate = mean)

        # Remove label colums and convert to matrix
        dat.alt <- dat.dcast
        rownames(dat.alt) <- dat.alt$Target
        dat.alt$Treatment <- NULL
        m.alt <- as.matrix(dat.alt)
        
        holder.r <<- rcorr(m.alt, type = "pearson")
        alt.r <- rcorr(m.alt, type = "spearman")
        cov.r <<- cov(m.alt)
        
        g <<- graph.adjacency(holder.r$r, mode = "undirected", 
                             weighted = TRUE, diag = FALSE)
        g.spear <<- graph.adjacency(holder.r$r, mode = "undirected", 
                                     weighted = TRUE, diag = FALSE)
        g.cov <<-graph.adjacency(cov.r, mode = "undirected", weighted = TRUE,
                                 diag = FALSE)
}

PlotTreatments <- function(cellLine, treatment, timePoint){
        
        
        # Filter data based on cell line, treatment, & time point
        dat.txt <- filter(dat, Treatment %in% treatment & `Cell Line` %in% 
                              cellLine & `Time Point` %in% timePoint)
        
        dat.avg <- dat.txt %>% 
                group_by(Target, `Time Point`, `Cell Line`, Treatment) %>%
                summarise_each(funs(mean), c(`Net Shift`))
        
        names(dat.avg)[names(dat.avg) == 'mean'] <- 'Net Shift'
        
        # Recast data as wide format from dataframe
        dat.dcast <<- dcast(dat.avg,
                            Target ~ `Cell Line` + `Time Point` + Treatment)
        
        # Remove label colums and convert to matrix
        dat.alt <- dat.dcast
        rownames(dat.alt) <- dat.alt$Target
        m.alt <<- as.matrix(dat.alt[,-1])
        
        holder.r <<- rcorr(m.alt)
        cov.r <<- cov(m.alt)
        
        g <<- graph.adjacency(holder.r$r, mode = "undirected", 
                              weighted = TRUE, diag = FALSE)
        
        g.cov <<-graph.adjacency(cov.r, mode = "undirected", weighted = TRUE,
                                 diag = FALSE)
}

PlotTargets(cellLine = c(6, 26), 
            treatment = c("Apitolisib", "GNE-317", "DMSO", 
                          "Erlotinib", "Palbociclib"), 
            timePoint = c(1, 24))

# Format edges
CustomGraph <- function(graph) {
        E(graph)$cor <- E(graph)$weight
        E(graph)$weight <- abs(E(graph)$cor)
        E(graph)$color <- ifelse(E(graph)$cor < 0, "darkblue", "grey50")
        E(graph)$width <- 3 * atanh(E(graph)$weight)
        V(graph)$label <- NA
        
}

graph <- g

CustomGraph(graph = graph)

V(graph)$label <- NA
colorList <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
               "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", 
               "grey50")

targetList <- c("pPDK1 Ser341", "pGSK3β Ser9", "pp70S6K Thr389", "pS6 Ser235/6", 
                "pRb Ser780", "pAkt Thr308", "pS6 Ser240/4", "pRb Ser807/11",
                "pmTOR Ser2448", "pMAPK Thr202/Tyr204", "pAkt Ser473", 
                "pSrc Tyr416", "hydroxy-HIF Pro564")

colors <- rep(colorList, each = 1)

V(graph)$color <- colors

graph <- delete.edges(graph, E(graph)[ abs(weight) < 0.8 ])


plot(graph, layout = layout.fruchterman.reingold)
plot(graph, layout = layout_as_star)
plot(graph, layout = layout_with_graphopt)
plot(graph, layout = layout_as_tree)

legend(x = -2, y = 1, targetList, pch=21,
       col="grey50", pt.bg= colorList ,bty="n", ncol=1)

corrplot(holder.r$r, method = "color")

hmcol <- colorRampPalette(c("darkblue", "white", "darkred"))(256)

heatmap(holder.r$r, col = hmcol)

PlotTreatments(cellLine = c(6, 26), 
             treatment = c("Erlotinib", "Palbociclib", "DMSO", 
                           "Apitolisib", "GNE-317"), 
             timePoint = c(1, 24))
CustomGraph(graph = g)

# Format edges
E(g)$cor <- E(g)$weight
E(g)$weight <- abs(E(g)$cor)
E(g)$color <- ifelse(E(g)$cor < 0, "darkblue", "grey")
E(g)$width <- atanh(E(g)$weight)

# Format vertices
#V(g)$label.color <- "black"
V(g)$label <- NA
colors <- c("orange", "blue")
V(g)$color <- ifelse(grepl("26", V(g)$name), "#ffffbf", "#91bfdb")


g <- delete.edges(g, E(g)[ abs(weight) < 0.6 ])

plot(g, layout = layout.fruchterman.reingold)
plot(g, layout = layout_with_kk)

legend(x=-1.75, y=-0.1, c("GBM 26","GBM 6"), pch=21,
       col="#777777", pt.bg=c('#ffffbf', "#91bfdb"), pt.cex=3, cex=3, bty="n")


corrplot(holder.r$r, method = "color")

hmcol <- colorRampPalette(c("darkblue", "white", "darkred"))(256)

heatmap(holder.r$r, col = hmcol)

# I need to find a way to get unique identifiers before dcast so that the 
# values don't get combined
