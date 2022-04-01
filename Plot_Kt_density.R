setwd("~/Documents/Scolaire/Stage_M1/Cetaceans/euler-2021-08-23-reduced/")

devtools::load_all("~/Documents/Scolaire/Stage_M1/RevGadgets")
library(deeptime)

Kt_mean <- readOBDP( start_time_trace_file="start_time_trace_N85_S200.txt", 
                     popSize_distribution_matrices_file="Kt_trace_N85_S200.txt", 
                     trees_trace_file="mcmc_OBDP_Cetaceans_wellMixedTrees.trees" )

p <- plotDiversityOBDP( Kt_mean,
                        xlab="Time (My)",
                        ylab="Number of lineages",
                        xticks.n.breaks=21,
                        col.Hidden="dodgerblue3",
                        col.LTT="gray25",
                        col.Total="forestgreen",
                        col.Hidden.interval="dodgerblue2",
                        col.Total.interval="darkolivegreen4",
                        palette.Hidden=c("transparent", "dodgerblue2", "dodgerblue3", "dodgerblue4", "black"),
                        palette.Total=c("transparent", "green4", "forestgreen", "black"),
                        line.size=0.7,
                        interval.line.size=0.5,
                        show.Hidden=F,
                        show.LTT=TRUE,
                        show.Total=TRUE,
                        show.intervals=TRUE,
                        show.densities=TRUE,
                        show.expectations=TRUE,
                        use.interpolate=TRUE )
p

q <- gggeo_scale(p, dat="periods", height=unit(1.3, "line"), abbrv=F, size=4.5, neg=T)
r <- gggeo_scale(q, dat="epochs", height=unit(1.1, "line"), abbrv=F, size=3.5, neg=T, skip=c("Paleocene", "Pliocene", "Pleistocene", "Holocene"))
s <- gggeo_scale(r, dat="stages", height=unit(1, "line"), abbrv=T, size=2.5, neg=T)
s
ggsave(filename="nbLineages_Cetacea_genera_Total.svg", plot=s, device="svg", width=12, height=6)
ggsave(filename="nbLineages_Cetacea_genera_Total.png", plot=s, device="png", width=12, height=6)

library("ggtree")

library("treeio")
tree <- read.beast(beast_file)
tree

tree <- read.mrbayes("mcmc_OBDP_Cetaceans.tre")
mit_rates <- read.table("mcmc_OBDP_Cetaceans_mit_rates.out", header=T, row.names=1)

ggtree(tree, aes(color=c(colMeans(mit_rates), rep(.015,16)))) +
  scale_color_continuous(low='darkgreen', high='red', name="Mit. Rates") +
  theme(legend.position="right")
