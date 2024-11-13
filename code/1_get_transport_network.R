#####################################################################
# Trans-African Network Between Cities > 100k and International Ports
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, s2, units, stplanr, sfnetworks, osrm, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

####################################
# Part 1: Important Cities and Ports
####################################

# CD13: CORRIDOR POINTE NOIRE – BRAZZAVILLE – OUESSO – BANGUI – N’DJAMENA (CD13) (CONGO-RCA-TCHAD) PHASE 2 -------------------------------

# Phase 1
# Mbaiki-Bangui (would need extra node for Mbaiki), 
# Ouesso-Pokola (Pokola is a small town to the east, not connected)
# Pokola-Gouga (CAR) (Where is this? -> minor upgrading)

# Phase 2
# Bossembélé-Baoro (3 segments, included)
# Mbaiki-Boda-Yaloke (Mbaiki-Boda included, Boda-Yaloke not yet, need to add Yaloke as waypoint)
# Bossembélé-Bossangoa (1 segment, included)
# Pokola-Gouga (CAR) Border road  (full asphalting)

# ROUTE KELO-PALA-LERE-FRONTIERE DU CAMEROUN (TCHAD-CAMEROUN) -------------

# Kélo-Pala-Léré-border Road (226km) (Kelo-Pala, 1 link, and Para-Sorawel (across the border, 1 link))

# KOUGOULEU-MEDOUNEU-AKURENAM (GABON-GUINEE-EQUATORIALE) Road development and Asphalting ----------------

# Kougouleu-Mela-Medounou-Akurenam (168 km) Road Asphalting (Present, but subset of a larger segment -> add Medounou-Akurenam as node?
# -> Is there a road connection between Medounou and Akurenam?

# ROUTE GAROUA BOULAI-BABOUA SUR LE CORRIDOR 2 (CAMEROUN-RCA) ------------------

# GAROUA-BOULAI-> BABOUA (54km) (a short piece part of the Gallo->GAROUA-BOULAI segment. Would need extra waypoint at Baboua)

# NDENDE-DOLISIE Road -------------------------

# Multiple segments (274km in total), but with nodes at Ndende and Dolisie
# -> 290.264 million € transport costs

# Also add points at Lastoursville
# Perhaps need AP24 -> has Yaloke and Pokola but not Gouga

CEMAC <- africamonitor::am_countries$ISO3[africamonitor::am_entities$CEMAC]

AP24 <- qread("/Users/sebastiankrantz/Documents/Data/Africapolis/2024/Africapolis_agglomeration_2024_geonames_centroids.qs") |>
        rm_stub("agglomeration_") |> rm_stub("_geo", FALSE) |>
        select(id, name, name_ascii, iso3c, lon, lat, 
               population = population_2025, built_up = built_up_2025, metropole,
               closest_metro = closest_metro_2025, dist_to_metro = dist_to_metro_2025)

AP24_CEMAC <- AP24 |> subset(iso3c %in% CEMAC) 
AP24_CEMAC |> count(iso3c)
# AP24_CEMAC |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population, name) |> mapview::mapview()

# First Deduplication
AP24_CEMAC_red <- largest_within_radius(AP24_CEMAC, c("lon", "lat"), "population", radius_km = 30) 
AP24_CEMAC_red <- largest_within_radius(AP24_CEMAC_red, c("lon", "lat"), "population", radius_km = 60)  

# Plotting
# AP24_CEMAC_red |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# AP24_CEMAC_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population, name) |> mapview::mapview()

# Now a second deduplication using the 50k threshold
AP24_CEMAC_g50k_red <- largest_within_radius(subset(AP24_CEMAC_red, population >= 5e4), c("lon", "lat"), "population", radius_km = 100) 

# Now ensuring strategically important places are present
waypoints <- AP24_CEMAC |> subset(name_ascii %ilike% "m'baiki|ouesso|pokola|gouga|boda|yaloke|bossembele|bossangoa|kelo|pala|figuil centre|garoua-boulai|baboua|ndende|dolisie|lastoursville")
# Need to add nodes at: kougouleu|medounou|akurenam (kougouleu) will be added
waypoints %<>% rowbind(rbind(data.frame(id = 1e6+1e5+1, name = "Medounou", name_ascii = "Medounou", iso3c = "GAB", lon = 10.796814, lat = 1.014809, population = 500),
                             data.frame(id = 1e6+2e5+2, name = "Akurenam", name_ascii = "Akurenam", iso3c = "GNQ", lon = 10.670471, lat = 1.036778, population = 500)), fill = TRUE)
AP24_CEMAC_g50k_red %<>% rowbind(subset(waypoints, id %!in% AP24_CEMAC_g50k_red$id))  
# Remove desert, island, and unnavigable towns
AP24_CEMAC_g50k_red %<>% subset(!name_ascii %ilike% "faya|arada|iriba|malabo|bredjing|impfondo")

# Now Adding populations within 30 km
AP24_CEMAC_g50k_red <- join_within_radius(AP24_CEMAC_g50k_red, AP24_CEMAC, c("lon", "lat"), radius_km = 30)
# AP24_CEMAC_g50k_red |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# AP24_CEMAC_g50k_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(name, population) |> mapview::mapview()


# Load International Ports from The World Bank: https://datacatalog.worldbank.org/search/dataset/0038118/Global---International-Ports -----
WBP_CEMAC <- fread("data/other_inputs/continental_african_ports.csv") |>
             subset(iso3c %in% CEMAC & name != "Malabo")
WBP_CEMAC_g50k <- WBP_CEMAC |> subset(outflows > 1e5)

## Visualization
# WBP_CEMAC_g50k |> with(plot(lon, lat, cex = sqrt(outflows)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WBP_CEMAC_g50k |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(outflows) |> mapview::mapview()

# Deduplication
WBP_CEMAC_g50k_red <- largest_within_radius(WBP_CEMAC_g50k, size = "outflows", radius_km = 100)
WBP_CEMAC_g50k_red <- join_within_radius(WBP_CEMAC_g50k_red, WBP_CEMAC, size = "outflows", radius_km = 99.9)

## Visualization
# WBP_CEMAC_g50k_red |> with(plot(lon, lat, cex = sqrt(outflows)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# WBP_CEMAC_g50k_red |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(outflows) |> mapview::mapview()

# Now Joining Cities and Ports ---------------------------------

cities_ports_dmat <- st_distance(with(AP24_CEMAC_g50k_red, s2_lnglat(lon, lat)), 
                                 with(WBP_CEMAC_g50k_red, s2_lnglat(lon, lat)))
fsum(cities_ports_dmat < as_units(30, "km"))

m_ind <- dapply(cities_ports_dmat < as_units(30, "km"), function(x) if (any(x)) which.max(x) else NA_integer_)

settransform(WBP_CEMAC_g50k_red, row = m_ind)
settransform(AP24_CEMAC_g50k_red, row = seq_along(name))

cities_ports <- AP24_CEMAC_g50k_red |> 
  join(subset(WBP_CEMAC_g50k_red, !is.na(row), row, port_locode = locode, port_name = name, port_status = status, outflows), on = "row") |> 
  select(-row)
rm(m_ind)

## Visualization
# cities_ports |> with(plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5)))
# cities_ports |> st_as_sf(coords = c("lon", "lat"), crs = st_crs(4326)) |> select(population) |> mapview::mapview()

# Raw matrix of routes:
dist_ttime_mats <- split_large_dist_matrix(select(cities_ports, lon, lat))

# Checks
all.equal(dist_ttime_mats$sources, dist_ttime_mats$destinations)
allv(diag(dist_ttime_mats$durations), 0)
allv(diag(dist_ttime_mats$distances), 0)
diag(cor(dist_ttime_mats$sources, select(cities_ports, lon, lat)))
s2_distance(with(dist_ttime_mats$sources, s2_lnglat(lon, lat)),
            with(select(cities_ports, lon, lat), s2_lnglat(lon, lat))) |> descr()
pwcor(unattrib(dist_ttime_mats$distances), unattrib(dist_ttime_mats$durations))

# Creating unique city identifier
settfm(cities_ports, city_country = stringi::stri_trans_general(paste(name_ascii, iso3c, sep = " - "), "latin-ascii"))
fndistinct(cities_ports)

# Saving
cities_ports |> fwrite("data/transport_network/cities_ports.csv")
dist_ttime_mats |> qsave("data/transport_network/cities_ports_dist_ttime_mats.qs")

# Cleanup
rm(list = ls())
source("code/helpers/helpers.R")
fastverse_conflicts()

####################################
# Part 2: Transport Network
####################################

# Reading again
cities_ports <- fread("data/transport_network/cities_ports.csv")
dist_ttime_mats <- qread("data/transport_network/cities_ports_dist_ttime_mats.qs")

# As spatial data frame
cities_ports_sf <- cities_ports |> st_as_sf(coords = c("lon", "lat"), crs = 4326)
# Distance between city centroids and route start/end points
descr(diag(st_distance(cities_ports_sf, st_as_sf(dist_ttime_mats$sources, coords = c("lon", "lat"), crs = 4326))))
# Create route-start-point version of the dataset
cities_ports_rsp_sf <- cities_ports_sf |> mutate(geometry = st_as_sf(dist_ttime_mats$sources, coords = c("lon", "lat"), crs = 4326)$geometry)

# Routes between all connections ------------------------------------------

cities_ports_dmat <- dist_ttime_mats$sources |> with(s2_lnglat(lon, lat)) |> st_distance() |> set_units("km")
diag(cities_ports_dmat) <- NA
cities_ports_dmat[upper.tri(cities_ports_dmat)] <- NA

# Routes to be calculated
routes_ind <- which(!is.na(cities_ports_dmat), arr.ind = TRUE)
nrow(routes_ind)


# Determining Ideal Hypothetical (US-Grade) Network 
# See: https://christopherwolfram.com/projects/efficiency-of-road-networks/

# US Route efficeincy = 0.843, thus mr = 1/0.843 = 1.186 
keep_routes <- !intercepted_routes(routes_ind, dist_ttime_mats$sources, NULL, alpha = 20, mr = 1/0.843) 
sum(keep_routes)

# mapview::mapview(subset(routes, keep_routes), zcol = "distance") +
#   mapview::mapview(st_as_sf(cities_ports, coords = c("lon", "lat"), crs = st_crs(4326)))

# Plot ideal Network
# <Figure 14: RHS>
with(cities_ports, {
  oldpar <- par(mar = c(0,0,0,0))
  plot(lon, lat, cex = sqrt(population)/1e3, pch = 16, col = rgb(0, 0, 0, alpha=0.5), 
       axes = FALSE, xlab = NA, ylab = NA, asp=1)
  for (r in mrtl(routes_ind[keep_routes, ])) { # Comment loop to get LHS of Figure 14
    lines(lon[r], lat[r])
  }
  par(oldpar)
})

dev.copy(pdf, "figures/trans_CEMAC_network_US_20deg.pdf", width = 5, height = 5)
dev.off()

# EU Route Efficiency  
keep_routes <- !intercepted_routes(routes_ind, dist_ttime_mats$sources, NULL, alpha = 45, mr = 1/0.767) 
sum(keep_routes)

# Plot ideal Network
with(cities_ports, {
  oldpar <- par(mar = c(0,0,0,0))
  plot(lon, lat, cex = 1.3, pch = 16, axes = FALSE, xlab = NA, ylab = NA, asp=1)
  for (r in mrtl(routes_ind[keep_routes, ])) { # & intersects == 2
    lines(lon[r], lat[r])
  }
  par(oldpar)
})

dev.copy(pdf, "figures/trans_CEMAC_network_EU_45deg.pdf", width = 10, height = 10)
dev.off()


# Fetching all Routes and Simplifying -------------------------------------------

# -> Better skip and load finished segments below

routes <- data.table(from_city_country = cities_ports$city_country[routes_ind[, 1]], 
                     to_city_country = cities_ports$city_country[routes_ind[, 2]], 
                     duration = NA_real_, 
                     distance = NA_real_, 
                     geometry = list())
# Fetch Routes
graphhopper::gh_set_api_url("https://graphhopper.com/api/1") # /1/route
for (i in seq_row(routes_ind)) {
  cat(i, " ")
  Sys.sleep(0.5)
  # route <- osrmRoute(ss(cities_ports, routes_ind[i, 1L], c("lon", "lat")),
  #                    ss(cities_ports, routes_ind[i, 2L], c("lon", "lat")), overview = "full") |>
  route <- graphhopper::gh_get_route(list(as.numeric(ss(cities_ports, routes_ind[i, 1L], c("lat", "lon"))),
                                          as.numeric(ss(cities_ports, routes_ind[i, 2L], c("lat", "lon")))),
                                     key="30e0a17e-c0ee-4929-b539-2ebec283da2e") |> graphhopper::gh_as_sf() |>
            tryCatch(error = function(e) NULL)
  if(is.null(route)) {
    cat(sprintf("\nroute %d from %s to %s could not be calculated\n", i, routes$from_city_country[i], routes$to_city_country[i]))
    next    
  }
  set(routes, i, 3:5, select(route, duration = time, distance, geometry))
}
routes <- routes |> subset(is.finite(duration)) |> st_as_sf(crs = st_crs(route))

# Saving
routes |> qsave("data/transport_network/routes_raw_gh.qs")

# Adding Gravity to Routes (https://en.wikipedia.org/wiki/Newton%27s_law_of_universal_gravitation)
routes <- qread("data/transport_network/routes_raw_gh.qs")
dmat <- st_distance(cities_ports_rsp_sf)
diag(dmat) <- NA
frange(dmat) 

colnames(dmat) <- rownames(dmat) <- cities_ports_rsp_sf$city_country
ddf <- pivot(qDF(dmat, "from"), "from", names = list("to", "sp_dist"), na.rm = TRUE) |> 
  join(select(cities_ports, from = city_country, from_pop = population)) |> 
  join(select(cities_ports, to = city_country, to_pop = population)) |> 
  mutate(gravity = as.double(from_pop) * to_pop / sp_dist / 1e6) # Figure in TEU for port cities??

routes <- routes |> 
  join(ddf, on = c(from_city_country = "from", to_city_country = "to"), drop = "x") |>
  mutate(gravity_rd = as.double(from_pop) * to_pop / distance / 1e6, 
         gravity_dur = as.double(from_pop) * to_pop / duration / 1e6)

# Intersecting Routes
segments <- overline2(routes |> mutate(passes = 1L), 
                      attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
rm(routes); gc()
segments %<>% ss(!is_linepoint(.)) %>% st_make_valid()

# Saving
segments |> qsave("data/transport_network/segments.qs")


# Loading Segments --------------------------------------------------------------------

segments <- qread("data/transport_network/segments.qs") 

# First Round of subdivision
segments <- rmapshaper::ms_simplify(segments, keep = 0.2, snap_interval = deg_m(500)) |> 
            subset(vlengths(geometry) >= 4) # |> st_cast("LINESTRING") |> st_make_valid() |> st_cast("LINESTRING")
segments <- overline2(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
if(any(is_linepoint(segments))) segments %<>% ss(!is_linepoint(.))
# plot(segments)
# mapview(segments) + mapview(rnet_get_nodes(segments))

# -> Feel free to skip and load final network below

# Creating Network
net <- as_sfnetwork(segments, directed = FALSE)
# plot(net)
summ_fun <- function(fun) list(passes = fun, gravity = fun, gravity_rd = fun, gravity_dur = fun, "ignore")
filter_smooth <- function(net) {
  net |> 
    tidygraph::activate("edges") |> 
    dplyr::filter(!tidygraph::edge_is_multiple()) |> 
    dplyr::filter(!tidygraph::edge_is_loop()) |> 
    tidygraph::convert(to_spatial_smooth, 
                       protect = cities_ports_rsp_sf$geometry,
                       summarise_attributes = summ_fun("mean"))
}


## Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)
# plot(net)
# mapview(st_geometry(net, "edges")) + mapview(st_geometry(net, "nodes"))

# Saving Smoothed Version (pre contraction)
net |> qsave("data/transport_network/net_smoothed.qs")
net <- qread("data/transport_network/net_smoothed.qs")

## Contracting network: Manually 
segments <- net |> activate("edges") |> tidygraph::as_tibble() |> mutate(.tidygraph_edge_index = NULL)
nodes <- nodes_max_passes(segments, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
nodes_clustered <- cluster_nodes_by_cities(nodes, cities_ports_rsp_sf, 
                                           city_radius_km = 20, 
                                           cluster_radius_km = 10, 
                                           algo.dbscan = FALSE,
                                           weight = "gravity_rd")
segments_contracted <- contract_segments(segments, nodes_clustered, attrib = c("passes", "gravity", "gravity_rd", "gravity_dur"))
# plot(segments_contracted)
# mapview(segments_contracted) + mapview(rnet_get_nodes(segments_contracted))

net <- as_sfnetwork(segments_contracted, directed = FALSE)

# Filtering, smoothing and spatial subdivision 
net <- filter_smooth(net)
net <- tidygraph::convert(net, to_spatial_subdivision)
net <- filter_smooth(net)
# plot(net)
# mapview(tidygraph::as_tibble(activate(net, "edges"))) + mapview(st_geometry(net, "nodes"))

## Saving 
net |> qsave("data/transport_network/net_discrete_final.qs")

## Loading Final Network -----------------------------------------

net <- qread("data/transport_network/net_discrete_final.qs")

## Plotting
plot(net)
nodes <- net |> activate("nodes") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_node_index = NULL)
edges <- net |> activate("edges") |> tidygraph::as_tibble() |> 
  mutate(.tidygraph_edge_index = NULL, log10_gravity = log10(gravity))

nodes$city_port <- rowSums(st_distance(nodes, cities_ports_rsp_sf) < as_units(20, "km")) > 0
sum(nodes$city_port)
descr(edges$gravity_rd)

# Needed throughout 
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

tmap_mode("plot")
tmap_options(raster.max.cells = 1e7)

# First a plot of just the routes
# <Figure 12: LHS>
pdf("figures/trans_CEMAC_network_actual.pdf", width = 7.5, height = 10)
tm_basemap("Esri.WorldGrayCanvas", zoom = 6) +
  tm_shape(segments) + tm_lines(col = "black") +
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2, fill = "orange2") +
  tm_layout(frame = FALSE)
dev.off()

# Now the plot with the discretized representation
# <Figure 12: RHS>
pdf("figures/trans_CEMAC_network_actual_discretized_gravity_plus_orig.pdf", width = 7.5, height = 10)
tm_basemap("Esri.WorldGrayCanvas", zoom = 6) +
  tm_shape(segments) + tm_lines(col = "black") +
  tm_shape(edges) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !city_port)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# Recalculating Routes Matrix for Network ----------------------------------------------------------

# -> Feel free to skip and load matrix below

nrow(nodes_coord)
all.equal(unattrib(nodes_coord), mctl(st_coordinates(nodes)))
# This can take a few minutes: now generating distance matrix of all nodes
dist_ttime_mats <- split_large_dist_matrix_gh(nodes_coord, verbose = TRUE)

# Checks
all.equal(dist_ttime_mats$sources, dist_ttime_mats$destinations)
allv(diag(dist_ttime_mats$durations), 0)
allv(diag(dist_ttime_mats$distances), 0)
diag(cor(dist_ttime_mats$sources, nodes_coord))
s2_distance(with(dist_ttime_mats$sources, s2_lnglat(lon, lat)),
            with(nodes_coord, s2_lnglat(lon, lat))) |> descr()
pwcor(unattrib(dist_ttime_mats$distances), unattrib(dist_ttime_mats$durations))
# Now finding places that are on islands (e.g. Zanzibar City): should not exist here
if(any(which(fnobs(dist_ttime_mats$durations) < 20))) stop("Found Islands")

dist_ttime_mats |> qsave("data/transport_network/net_dist_ttime_mats.qs")

# Loading again -----------------------------------

dist_ttime_mats <- qread("data/transport_network/net_dist_ttime_mats.qs")

# Check
all.equal(st_as_sf(net, "nodes")$geometry, nodes$geometry)

# Making symmetric
sym_dist_mat <- (dist_ttime_mats$distances + t(dist_ttime_mats$distances)) / 2
sym_time_mat <- (dist_ttime_mats$durations + t(dist_ttime_mats$durations)) / 2

# Add average distance and travel time to edges
edges_ind <- edges |> qDF() |> select(from, to) |> qM()
edges$sp_distance <- st_length(edges)
edges$distance <- sym_dist_mat[edges_ind]
edges$duration <- sym_time_mat[edges_ind]

# Loading, Integrating and Plotting Additional Connections ---------------------------------------------------

# -> To generate the 'add_links' file, execute '7.1_add_links.R' in a clean R session

add_links <- qread("data/transport_network/add_links_network_30km_alpha45_mrEU_fmr15.qs")
add_links_df <- line2points(add_links)
dmat <- st_distance(nodes$geometry, add_links_df$geometry)
add_links_df$node <- dapply(dmat, which.min)
add_links_df$geometry <- nodes$geometry[add_links_df$node]
add_links <- add_links_df |> group_by(id) |> 
              summarise(from = ffirst(node),
                        to = flast(node),
                        geometry = st_combine(geometry)) |> st_cast("LINESTRING")
# Checks
all(line2df(add_links) %>% select(fx, fy) %in% qDF(st_coordinates(nodes)))
all(line2df(add_links) %>% select(tx, ty) %in% qDF(st_coordinates(nodes)))
all(select(line2df(add_links), fx:ty) %!in% select(line2df(edges), fx:ty))
rm(dmat, add_links_df)

tmap_mode("plot")

# Same as Figure 15 but with discrete edges (Figure 15 with real roads is obtained below)
pdf("figures/trans_CEMAC_network_actual_discretized_gravity_new_roads_Esri.WorldStreetMap.pdf", width = 7.5, height = 10)
tm_basemap("Esri.WorldStreetMap", zoom = 7) +
  tm_shape(segments) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 1) + 
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !city_port)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()



# Now Adding Populations ----------------------------------------------------------

# Distance Matrix
dmat <- nodes |> st_distance(st_as_sf(AP24_CEMAC, coords = c("lon", "lat"), crs = 4326)) 
# Remove cities part of port cities
used <- dapply(dmat[nodes$city_port, ] < as_units(30, "km"), any)
table(used)
dmat <- dmat[!nodes$city_port, !used]
dim(dmat)
dmat[dmat >= as_units(30, "km")] <- NA
col <- numeric(nrow(dmat))
m <- dapply(dmat, function(x) if (allNA(x)) col else replace(col, which.min(x), 1)) # Assigning each city to closest node
nodes$population <- NA
nodes$population[!nodes$city_port] <- m %*% AP24_CEMAC$population[!used] # Summing populations of closest cities
nodes$city_country <- NA_character_
nodes$city_country[!nodes$city_port] <- qDT(AP24_CEMAC) |> 
  extract(!used, paste(name_ascii, iso3c, sep = " - ")) |> 
  extract(dapply(m %r*% AP24_CEMAC$population[!used], 
                 function(x) if(any(x > 0.1)) which.max(x) else NA_integer_, MARGIN = 1))
# Now adding port cities (closest matches)
ind <- dapply(st_distance(cities_ports_rsp_sf, nodes), function(x) if (any(x < 20e3)) which.min(x) else NA_integer_)
nodes$population[!is.na(ind)] <- cities_ports_rsp_sf$population[na_rm(ind)]
nodes$city_country[!is.na(ind)] <- cities_ports_rsp_sf$city_country[na_rm(ind)]
# Cleanup
rm(col, m, ind, used, dmat)

# Ratios
sum(nodes$population) / sum(cities_ports_rsp_sf$population)
sum(nodes$population) / sum(AP24_CEMAC$population)
(sum(nodes$population > 0) - nrow(cities_ports_rsp_sf)) / (nrow(nodes) - nrow(cities_ports_rsp_sf))

# Saving all necessary objects in an RData file ---------------------------------------------------

# First adding info to edges
tmp <- net |> st_as_sf("edges") |> atomic_elem() |> qDT() |> 
  join(atomic_elem(edges), overid = 2L) |> 
  select(-from, -to, -.tidygraph_edge_index)

net %<>% activate("edges") %>% dplyr::mutate(tmp)
rm(tmp)

save(nodes, edges, edges_ind, nodes_coord, net, add_links, AP24_CEMAC, 
     cities_ports, cities_ports_sf, cities_ports_rsp_sf, 
     dist_ttime_mats, sym_dist_mat, sym_time_mat, 
     file = "data/transport_network/trans_CEMAC_network.RData")


##########################################################
# Fetching Simplified Routes (Edges) for Visual Inspection
##########################################################

load("data/transport_network/trans_CEMAC_network.RData")

# This is the previous Plot
tm_basemap("Esri.WorldGrayCanvas", zoom = 6) +
  tm_shape(edges) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 1) + 
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, !city_port)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

# dev.copy(pdf, "figures/trans_CEMAC_network_actual_discretized_gravity_new_roads.pdf", 
#          width = 10, height = 10)
# dev.off()

# Fetching Simplified Routes

# -> Feel free to skip and load result below

edges_real <- edges |> qDT() |> 
  transform(geometry = list(NULL), distance = NA_real_, duration = NA_real_)
nodes_coord <- st_coordinates(nodes) |> qDF() |> setNames(c("lon", "lat"))

# Fetch Routes
graphhopper::gh_set_api_url("https://graphhopper.com/api/1") # /1/route
for (i in seq_row(edges_ind)) {
  cat(i, " ")
  Sys.sleep(0.5)
  # route <- osrmRoute(ss(nodes_coord, edges_ind[i, 1]),
  #                    ss(nodes_coord, edges_ind[i, 2]), overview = "simplified")
  route <- graphhopper::gh_get_route(list(as.numeric(ss(nodes_coord, edges_ind[i, "from"], c("lat", "lon"))), 
                                          as.numeric(ss(nodes_coord, edges_ind[i, "to"], c("lat", "lon")))), 
                                     key="30e0a17e-c0ee-4929-b539-2ebec283da2e") |> graphhopper::gh_as_sf() |>
    tryCatch(error = function(e) NULL)
  if(is.null(route)) {
    cat(sprintf("\nroute %d from %s to %s could not be calculated", i, nodes$city_country[edges$from[i]], nodes$city_country[edges$to[i]]))
    next    
  }
  set(edges_real, i, c("duration", "distance", "geometry"), 
      select(route, duration = time, distance, geometry))
}
edges_real <- edges_real |> st_as_sf(crs = st_crs(route)) |> st_make_valid()
edges_real |> qsave("data/transport_network/edges_real.qs")
rm(route, i)

edges_real %<>% rmapshaper::ms_simplify(keep = 0.1)
edges_real %<>% st_make_valid()
edges_real |> qsave("data/transport_network/edges_real_simplified.qs")


# Draw the updated plot
edges_real <- qread("data/transport_network/edges_real_simplified.qs")

# <Figure 15>
pdf("figures/trans_CEMAC_network_actual_discretized_gravity_new_roads_real_edges.pdf", width = 7.5, height = 10)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(add_links) + tm_lines(col = "green4", lwd = 1) + # limegreen # , lty = "twodash"
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.15) +
  tm_shape(subset(nodes, !city_port & population > 0)) + tm_dots(size = 0.1, fill = "grey20") +
  tm_shape(subset(nodes, !city_port & population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) #, inner.margins = c(0.1, 0.1, 0.1, 0.1))
dev.off()

