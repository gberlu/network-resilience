# getting started ----

# load packages
library(network) # network analysis
library(sna)
library(ergm)
library(haven) # data management
library(tidyverse)
library(kableExtra) # print tables
library(texreg) # convert statistical model outputs to tables
library(wesanderson) # color palette

# create function to check network
# from https://www.mjdenny.com/Preparing_Network_Data_In_R.html
Check_Network_Integrity <- function(network_object, # the network object you created
                                    n_edges_to_check = 10 # defaults to 10 edges
){
  # make sure you are providing a network object to the function
  if(class(network_object) != "network"){
    stop("You must provide this function with an object of class network!")
  }
  # get network information and edges
  num_nodes <- length(network_object$val)
  edgelist <- as.matrix(network_object, matrix.type = "edgelist")
  names <- get.vertex.attribute(network_object, "vertex.names")
  # make sure you do not ask for more edges than are in your network
  if(n_edges_to_check > length(edgelist[, 1])){
    n_edges_to_check <- length(edgelist[, 1])
  }
  # select a sample of edges
  edges_to_check <- sample(x = 1:length(edgelist[, 1]),
                           size = n_edges_to_check,
                           replace = FALSE)
  # print out the edges to check
  cat("Check source data file to make sure that the following edges exist:\n\n")
  for(i in 1:n_edges_to_check){
    cat(names[edgelist[edges_to_check[i], 1]], "<-->", names[edgelist[edges_to_check[i], 2]], "\n")
  }
  cat("\nIf any of these edgest do not match your source data files, there was likely a problem reading in your data, and you should revisit the process.\n")
}

# read in co-offending data, and store them in a list
mat_list <- list(read.csv("data/co-offending-1.csv", header = FALSE, sep = ","), # edgelist for network at time 1
                 read.csv("data/co-offending-2.csv", header = FALSE, sep = ","), # edgelist for network at time 2
                 read.csv("data/co-offending-3.csv", header = FALSE, sep = ",")) # edgelist for network at time 3
head(mat_list[[1]]) # check data

# read in attribute data, and store them in a list
att_list <- list(read.csv("data/attributes-1.csv", header = TRUE, sep = ","), # attributes for network at time 1
                 read.csv("data/attributes-2.csv", header = TRUE, sep = ","), # attributes for network at time 2
                 read.csv("data/attributes-3.csv", header = TRUE, sep = ",")) # attributes for network at time 3
head(att_list[[1]]) # check data

# create network objects, store them in a list, and add data
# net_list <- ls() # create empty list
# for (i in 1:3){ 
#   net_list[i] <- list(network(mat_list[[i]], # create network objects from edgelists in mat_list
#                               directed = FALSE, # networks are undirected
#                               matrix.type = "e", # data are stored in edgelists
#                               vertex.attr = att_list[[i]][, 2:11], # add attributes from att_list
#                               vertex.attrnames = colnames(att_list[[i]][, 2:11]))) # and attributes' names
# }
# old code above doesn't seem to work anymore - nodes are not listed in alphabetical order...
# ...so new code is below - it ensures nodes are listed in alphabetical order

# create empty network objects
net_list <- list(network.initialize(length(att_list[[1]]$id), directed = FALSE), # initialize network with right no. nodes (length of list of node labels at time 1)
                 network.initialize(length(att_list[[2]]$id), directed = FALSE), # repeat for time 2
                 network.initialize(length(att_list[[3]]$id), directed = FALSE)) # and then again for time 3

# add names to nodes and values to edges to get co-offending networks
for (i in 1:3){
  network.vertex.names(net_list[[i]]) <- as.vector(att_list[[i]]$id) # assign names to nodes
  net_list[[i]][as.matrix(mat_list[[i]][, 1:2])] <- 1 # add in edges from the edgelist to the network object
  set.vertex.attribute(net_list[[i]], attrname = colnames(att_list[[i]][2:11]), value = att_list[[i]][, 2:11]) # add attributes and their names from att_list
}
Check_Network_Integrity(net_list[[1]]) # check new network

# descriptive statistics ----

# create empty matrix
stats <- matrix(ncol = 3, nrow = 28, # matrix has 3 columns and 28 rows
                dimnames = list(c("Network size", "Edge count", "Density", "Mean degree", # add row names
                                  "Maximum degree", "Degree centralisation", "Clustering coefficient", 
                                  "Nationality", "Italian", "Moroccan", "Romanian", "Other/unknown",
                                  "Task", "Supplier", "Trafficker", "Courier", "Buyer", "Support", "Unknown", 
                                  "Role in the drug supply chain", "Supply", "Importation", "Distribution", "Unknown",
                                  "Status", "High", "Medium", "Low"), 
                                c("Phase 1", "Phase 2", "Phase 3"))) # add column names

# calculate descriptive statistics, and store them in the matrix
for (i in 1:3){
  stats[, i] <- c(network.size(net_list[[i]]), # network size
                  network.edgecount(net_list[[i]]), # network edge count
                  round(gden(net_list[[i]]), 2), # density
                  round(mean(degree(net_list[[i]], gmode = "graph")),2), # mean degree
                  max(degree(net_list[[i]], gmode = "graph")), # max. degree
                  round(centralization(net_list[[i]], degree, mode = "graph", normalize = TRUE), 2), # degree centralization
                  round(gtrans(net_list[[i]], mode = "graph", measure = "weak"), 2), # clustering coefficient
                  "",
                  sum(att_list[[i]]$nationality == "Italian"), # actors' nationality
                  sum(att_list[[i]]$nationality == "Moroccan"),
                  sum(att_list[[i]]$nationality == "Romanian"),
                  sum(att_list[[i]]$nationality == "Colombian" | att_list[[i]]$nationality == "Spanish" | att_list[[i]]$nationality == "Unknown"),
                  "",
                  sum(att_list[[i]]$task == "Supplier"), # actors' task
                  sum(att_list[[i]]$task == "Trafficker"),
                  sum(att_list[[i]]$task == "Courier"),
                  sum(att_list[[i]]$task == "Buyer"),
                  sum(att_list[[i]]$task == "Support"),
                  sum(att_list[[i]]$task == "Unknown"),
                  "",
                  sum(att_list[[i]]$role == "Supply"), # actors' role in the drug supply chain
                  sum(att_list[[i]]$role == "Importation"),
                  sum(att_list[[i]]$role == "Distribution"),
                  sum(att_list[[i]]$role == "Unknown/Peripheral involvement"),
                  "",
                  sum(att_list[[i]]$status == "High"), # actors' status
                  sum(att_list[[i]]$status == "Medium"),
                  sum(att_list[[i]]$status == "Low"))
}

# print and save descriptive statistics
stats # print table
View(stats) # another way to print table

write.table(stats, "output/table1_stats.txt", sep="\t") # save table in a .txt file with values separated by tabs

# network visualization ----

# list all nodes' IDs (n = 128)
ids_all <- sort(unique(c(network.vertex.names(net_list[[1]]), network.vertex.names(net_list[[2]]), network.vertex.names(net_list[[3]])))) 

# create empty networks with isolates (n = 128)
net_all_list <- list(network.initialize(length(ids_all)), network.initialize(length(ids_all)), network.initialize(length(ids_all))) # initialize a network with all nodes (n = 128)

# assign names to nodes and values to edges
for (i in 1:3){
  network.vertex.names(net_all_list[[i]]) <- ids_all # assign names to nodes
  net_all_list[[i]][as.matrix(mat_list[[i]][, 1:2])] <- 1
}

# create network with all nodes (n = 128)
net_all <- network.initialize(length(ids_all)) # initialize a network with all nodes (n = 128)
network.vertex.names(net_all) <- ids_all # assign names to nodes
net_all[as.matrix(mat_list[[1]][, 1:2])] <- 1 # assign values to edges in network
net_all[as.matrix(mat_list[[2]][, 1:2])] <- 1
net_all[as.matrix(mat_list[[3]][, 1:2])] <- 1
Check_Network_Integrity(net_all) # check new network

# save coordinates for all nodes
coords_all <- gplot(net_all, usearrows = FALSE)

# set color palette
names(wes_palettes) # check color palettes
pal <- wes_palette("Royal1", 4, type = "discrete") 

# use different color for N3
vect_ids <- rep(FALSE, times = length(network.vertex.names(net_all))) # create attribute vector
lab_ids <- rep("", times = length(network.vertex.names(net_all))) # create label vector
for (i in 3) {
  vect_ids[i] <- TRUE # make N3' attribute TRUE
  lab_ids[i] <- network.vertex.names(net_all)[i] # add labels
}

# set node colors
netcolours <- ifelse(vect_ids == TRUE, pal[2], pal[1]) # red for N3, grey for all other nodes

# print and save graphs
png (filename = "output/figure1_plot.png", height = 1000, width = 2500, bg = "white") # save figure in a .png file
par(mfrow=c(1,3)) # print graphs in one row and three columns

gplot(net_all_list[[1]], # plot network at time 1
      usearrows = FALSE, # don't use arrow heads
      label = lab_ids, # use label only for N3
      displayisolates = FALSE, # don't display isolates
      usecurve = TRUE, # edges are slightly curved
      edge.curve = .01, 
      vertex.col = netcolours, # color N3 in red, other nodes in grey
      vertex.cex = 2, # set nodes' size
      label.cex = 4.5, # set labels' size
      edge.lwd = 1, # set edge width
      edge.col = pal[1], # use grey edges
      coord = coords_all, # use same coordinates for all three networks
      main = "Phase 1", # include figure title
      cex.main = 6) # set title's size

gplot(net_all_list[[2]], usearrows = FALSE, label = lab_ids, displayisolates = FALSE, # repeat for network at time 2
      usecurve = TRUE, edge.curve = .01, vertex.col = netcolours, vertex.cex = 2, label.cex = 4.5, 
      edge.lwd = 1, edge.col = pal[1], coord = coords_all, main = "Phase 2", cex.main = 6)

gplot(net_all_list[[3]], usearrows = FALSE, label = lab_ids, displayisolates = FALSE, # and then again for network at time 3
      usecurve = TRUE, edge.curve = .01, vertex.col = netcolours, vertex.cex = 2, label.cex = 4.5, 
      edge.lwd = 1, edge.col = pal[1], coord = coords_all, main = "Phase 3", cex.main = 6)

par(mfrow=c(1,1)) # go back to default setting
dev.off()

# changes across phases ----

# create empty matrix
overlap <- matrix(ncol = 2, nrow = 11, # matrix has 2 columns and 11 rows
                  dimnames = list(c("No. observed actors", "Actors at t", "Actors at t + 1", # add row names
                                    "Size difference", "Joint", "Combined", "Jaccard index",
                                    "No. observed ties", "Ties at t", "Ties at t + 1", "Jaccard index"),
                                  c("Phase 1 to 2", "Phase 2 to 3"))) # add column names

# combine edgelists and identify edges in each of the three networks
mat_all <- bind_rows(mat_list[[1]], mat_list[[2]], mat_list[[3]]) %>% # combine the three edgelists
  distinct(V1, V2, .keep_all = TRUE) %>% # remove duplicates based on first and second column
  select(-3) %>% # remove third column
  left_join(mat_list[[1]], by = c("V1", "V2")) %>% # create new dummy variable where 1 == edge is present in network at time 1
  left_join(mat_list[[2]], by = c("V1", "V2")) %>%
  left_join(mat_list[[3]], by = c("V1", "V2")) %>%
  rename(time1 = V3.x, time2 = V3.y, time3 = V3) %>% # rename new columns
  replace_na(list(time1 = 0, time2 = 0, time3 = 0)) # replace NAs with zeros

# calculate change statistics
for (i in 1:2){
  joint_nodes <- length(intersect(network.vertex.names(net_list[[i]]), network.vertex.names(net_list[[i+1]]))) # no. nodes shared by two networks
  combined_nodes <- length(union(network.vertex.names(net_list[[i]]), network.vertex.names(net_list[[i+1]]))) # total no. nodes in both networks (shared and not shared)
  joint_edges <- dim(filter(mat_all, mat_all[,i+2] == 1 & mat_all[,i+3] == 1))[1] # no. edges shared by two networks
  combined_edges <- dim(filter(mat_all, mat_all[,i+2] == 1 | mat_all[,i+3] == 1))[1] # total no. edges in both networks (shared and not shared)
  overlap[, i] <- c("", network.size(net_list[[i]]), network.size(net_list[[i+1]]), # size of each networks
                    network.size(net_list[[i+1]]) - network.size(net_list[[i]]), # size difference
                    joint_nodes, combined_nodes, # from above - nodes shared by / in both networks
                    round(joint_nodes / combined_nodes * 100, 1), # jaccard index for node overlap
                    "", sum(mat_all[,i+2]), sum(mat_all[,i+3]), # no. edges in each networks
                    round(joint_edges / combined_edges * 100, 1)) # jaccard index for edge overlap
}

# print and save descriptive statistics
overlap # print table
View(overlap) # another way to print table

write.table(overlap, "output/table2_changes.txt", sep="\t") # save table in a .txt file with values separated by tabs

# centrality ----

# check N3's degree centrality score
ord_by_deg <- order(degree(net_list[[1]], gmode = "graph"), decreasing = TRUE) # order nodes at time 1 by degree centrality
cbind(network.vertex.names(net_list[[1]]), degree(net_list[[1]], gmode = "graph")) [ord_by_deg, ] # print node labels and degree scores

# repeat for betweenness centrality
ord_by_bet <- order(betweenness(net_list[[1]], gmode = "graph"), decreasing = TRUE) # order nodes at time 1 by betweenness centrality
cbind(network.vertex.names(net_list[[1]]), betweenness(net_list[[1]], gmode = "graph")) [ord_by_bet, ]

# ergm ----

# read in kinship data, and store them in a list
mat_kin_list <- list(read.csv("data/kinship-1.csv", header = FALSE, sep = ","), read.csv("data/kinship-2.csv", header = FALSE, sep = ","), read.csv("data/kinship-3.csv", header = FALSE, sep = ","))
head(mat_kin_list[[1]]) # check data

# create empty network objects
kin_list <- list(network.initialize(network.size(net_list[[1]]), directed = FALSE), network.initialize(network.size(net_list[[2]]), directed = FALSE), network.initialize(network.size(net_list[[3]]), directed = FALSE)) # initialize network with right no. nodes

# add names to nodes and values to edges to get kinship networks
for (i in 1:3){
  network.vertex.names(kin_list[[i]]) <- as.vector(att_list[[i]]$id) # assign names to nodes
  kin_list[[i]][as.matrix(mat_kin_list[[i]][, 1:2])] <- 1 # add in edges from the edgelist to the network object
}
Check_Network_Integrity(kin_list[[1]]) # check new network

# read in affiliation data, and store them in a list
mat_clan_list <- list(read.csv("data/clan-1.csv", header = FALSE, sep = ","), read.csv("data/clan-2.csv", header = FALSE, sep = ","), read.csv("data/clan-3.csv", header = FALSE, sep = ","))
head(mat_clan_list[[1]]) # check data

# create empty network objects
clan_list <- list(network.initialize(network.size(net_list[[1]]), directed = FALSE), network.initialize(network.size(net_list[[2]]), directed = FALSE), network.initialize(network.size(net_list[[3]]), directed = FALSE)) # initialize network with right no. nodes

# add names to nodes and values to edges to get affiliation networks
for (i in 1:3){
  network.vertex.names(clan_list[[i]]) <- as.vector(att_list[[i]]$id) # assign names to nodes
  clan_list[[i]][as.matrix(mat_clan_list[[i]][, 1:2])] <- 1 # add in edges from the edgelist to the network object
}
Check_Network_Integrity(clan_list[[1]]) # check new network

# set seed
set.seed(8)

# check data
ls()

# check node attributes' names
list.vertex.attributes(net_list[[1]])

# check node attributes' categories
table(net_list[[1]]%v%"role") # [1] Distribution, [2] Importation, [3] Supply, [4] Unknown/Peripheral involvement
table(net_list[[1]]%v%"task") # [1] Buyer, [2] Courier, [3] Supplier, [4] Support, [5] Trafficker, [6] Unknown
table(net_list[[1]]%v%"status") # [1] High, [2] Low, [3] Medium
table(net_list[[1]]%v%"affiliation.2") # [1] Moroccan clan, [2] Ndrangheta, [3] No affiliation, [4], Romanian clan
table(net_list[[1]]%v%"nationality") # [1] Italian, [2] Moroccan, [3] Romanian, [4] Spanish, [5] Unknown
table(net_list[[1]]%v%"targeted") # [1] 0 (No), [2] 1 (Yes)

# check homophily
summary(net_list[[1]] ~ nodematch("task"))
summary(net_list[[1]] ~ nodematch("role"))

# run ergms of the three networks
fit_list <- list() # create empty list
alpha_list <- c(0.05, 0.05, 0.15)

for (i in 1:3){
  fit_list[[i]] <- ergm(net_list[[i]] ~ edges
                        + nodefactor("status", levels = -(2:3)) # individual attribute: high status
                        + nodefactor("affiliation.2", levels = 2) # individual attribute: 'Ndrangheta member
                        + nodefactor("task", levels = 5) # individual attribute: trafficker
                        
                        + nodematch("task") # uniform homophily: task
                        + nodematch("role") # uniform homophily: role
                        
                        + edgecov(kin_list[[i]]) # multiplexity: kinship ties
                        + edgecov(clan_list[[i]]) # multiplexity: shared criminal affiliation
                        
                        + gwdegree(0.25, fixed = TRUE) # preferential attachment
                        + gwesp(alpha_list[[i]], fixed = TRUE) # triadic closure
                        + gwdsp(alpha_list[[i]], fixed = TRUE) # indirect ties
                        
                        + nodematch("nationality") # uniform homophily: nationality
                        + nodefactor("targeted") # individual attribute: targeted
                        
                        , control = control.ergm(MCMC.samplesize = 4096, MCMC.interval = 8192)
                        , verbose = TRUE
                        )
  print(summary(fit_list[[i]]))
}

# assess model convergence of the three models
for (i in 1:3){
  pdf(paste("output/diag-", i, ".pdf", sep = ""))
  mcmc.diagnostics(fit_list[[i]])
  dev.off()
}

# convert ergm outputs to table
models <- screenreg(fit_list, # ergms of the three networks
                    single.row = TRUE, # place coefficients and standard errors in the same line
                    custom.model.names = c("Phase 1", "Phase 2", "Phase 3"), # set model names
                    custom.coef.map = list("edges" = "Edges", # reorder, group, and rename coefficients
                                           "nodefactor.status.High" = "Activity spread (high status)",
                                           "nodefactor.task.Trafficker" = "Activity spread (trafficker)",
                                           "nodefactor.affiliation.2.Ndrangheta" = "Activity spread ('Ndrangheta member)",
                                           "nodematch.task" = "Homophily by task",
                                           "nodematch.role" = "Homophily by role",
                                           "edgecov.kin_list[[i]]" = "Multiplex ties (kinship)",
                                           "edgecov.clan_list[[i]]" = "Multiplex ties (formal org.)",
                                           "gwdeg.fixed.0.25" = "Preferential attachment",
                                           "gwesp.fixed.0.05" = "Triadic closure",
                                           "gwesp.fixed.0.15" = "Triadic closure",
                                           "gwdsp.fixed.0.05" = "Indirect connections",
                                           "gwdsp.fixed.0.15" = "Indirect connections",
                                           "nodematch.nationality" = "Homophily by nationality",
                                           "nodefactor.targeted.1" = "Activity spread (targeted)")
                    )

# print and save descriptive statistics
models # print table

write.table(models, "output/table4_models.txt", sep="\t") # save table in a .txt file with values separated by tabs

# goodness of fit ----

# test whether the models fit the data
gof_list <- list() # create empty list
for (i in 1:3){
  gof_list[[i]] <- gof(fit_list[[i]], verbose = TRUE, interval = 5e+4) # check gof
  png(paste("output/gof-", i, ".png", sep = ""), height = 400, width = 1600, bg = "white") # save gof as png file
  par(mfrow = c(1,4)) 
  plot(gof_list[[i]], cex.lab = 1.6, cex.axis = 1.6, plotlogodds = TRUE) # plot gof results
  par(mfrow = c(1,1)) 
  dev.off()
}

# save data ----
save.image("output/all_results.RData")
