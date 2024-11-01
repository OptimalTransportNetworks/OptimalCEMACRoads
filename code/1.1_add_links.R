#############################################
# Adding High-Value Links to Existing Network
#############################################

# Note: This is a self-contained section due to the use of specific indices for manual adjustments. 
# It uses a previous version of the network identical to the current/generated version but where the nodes and edges are in a different (random) order. 
# The output is a spatial data frame 'add_links' with the final proposed links. These can be added to the current/generated version.

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, sfnetworks, stplanr, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

net <- qread("data/transport_network/net_discrete_final.qs")
dist_ttime_mats <- qread("data/transport_network/net_dist_ttime_mats.qs")
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL)

nodes_dmat <- st_distance(nodes) |> set_units("km")
diag(nodes_dmat) <- NA
nodes_dmat[upper.tri(nodes_dmat)] <- NA

# Routes to be calculated
routes_ind <- which(!is.na(nodes_dmat), arr.ind = TRUE)
nrow(routes_ind)

nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

keep_routes <- !intercepted_routes(routes_ind, nodes_coord, sym_dist_mat, feasible = TRUE, alpha = 45, mr = 1/0.843, fmr = 1.3) # US: 0.843 EU: 0.767
sum(keep_routes)

add_links <- with(nodes_coord, lapply(mrtl(routes_ind[keep_routes, ]), function(r) st_linestring(cbind(lon[r], lat[r])))) |> 
              st_sfc(crs = 4326) |> st_as_sf()
add_links_df <- line2df(add_links)[-1] %*=% 1000 %>% dapply(as.integer)
edges_df <- line2df(edges)[-1] %*=% 1000 %>% dapply(as.integer)
m <- add_links_df %in% edges_df | add_links_df %in% gv(edges_df, c(3:4, 1:2))
nrow(add_links) - sum(m)
add_links <- add_links[!m, ]
rm(add_links_df, edges_df)

add_links |> qsave("data/transport_network/add_links_network_30km_alpha45_mrUS_fmr13.qs")
add_links <- qread("data/transport_network/add_links_network_30km_alpha45_mrUS_fmr13.qs")

# Manual Adjustments: Removing links crossing waterbodies or protected areas
# mapview(qread("data/transport_network/segments.qs"), map.types = c(mapviewGetOption("basemaps"), "Esri.WorldStreetMap", "Esri.WorldTerrain")) + mapview(nodes) + mapview(add_links[-rem, ][-rem2, ][-rem3, ][-rem4, ][-rem5, ][-rem6, ], color = "green") # 
rem <- c(57, 88, 89, 140, 144, 145, 146, 137, 130, 142, 73, 112, 102, 100, 90, 81, 71, 125, 34, 4, 18, 2, 8, 17, 
            75, 74, 61, 65, 51, 56, 38, 33, 39, 49, 10, 13, 14)
rem2 <- c(52, 33, 35, 16, 15, 91, 83, 85, 75, 79, 100)
rem3 <- c(38, 71, 75, 76, 77, 83)
rem4 <- c(7, 26, 27, 12, 80)
rem5 <- c(47, 76, 47, 55, 57, 72)
rem6 <- c(86, 82, 23, 41)
add <- list(c(12, 57), c(57, 94), c(117, 156), c(154, 166))
add_links <- rbind(add_links[-rem, ][-rem2, ][-rem3, ][-rem4, ][-rem5, ][-rem6, ],
                    with(nodes_coord, lapply(add, function(r) st_linestring(cbind(lon[r], lat[r])))) |> st_sfc(crs = 4326) |> st_as_sf())

nrow(add_links)
if(nrow(add_links) != 86) stop("File should have 481 rows")
add_links |> qsave("data/transport_network/add_links_network_30km_alpha45_mrUS_fmr13_adjusted.qs")
# mapview(qread("data/transport_network/segments.qs"), map.types = c(mapviewGetOption("basemaps"), "Esri.WorldStreetMap", "Esri.WorldTerrain")) + mapview(nodes) + mapview(add_links, color = "green") 
