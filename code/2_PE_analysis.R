#####################################################################
# African Transport Network: Partial Equilibrium Analysis
#####################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"), nthreads = 4)
fastverse_extend(qs, sf, units, sfnetworks, tmap, install = TRUE)
source("code/helpers/helpers.R")
fastverse_conflicts()

GADM0_africa <- qread("data/other_inputs/GADM0_africa_simplified.qs")
CEMAC <- africamonitor::am_countries$ISO3[africamonitor::am_entities$CEMAC]

# All essential objects from previous sections
# Need only load("data/transport_network/trans_CEMAC_network.RData"), the result of 1_get_transport_network.R
# The .RData file with the _param suffix is created at the end of this file and adds a 'parameterized' version 
# of the network using results computed in this file. It thus just has additional objects and also works as input. 
load("data/transport_network/trans_CEMAC_network_google.RData")
edges$duration %/=% 60 # From seconds to minutes -> also with google
dist_ttime_mats$durations %/=% 60
edges_real <- qread("data/transport_network/edges_real_simplified.qs") |>
  rmapshaper::ms_simplify(keep = 0.1) |> st_make_valid()

# Average route efficiency per edge
aere <- unattrib(mean(edges$sp_distance / edges$distance)) 
print(aere)

# Adding distance to new edges based on average edge route efficiency
settfm(add_links, sp_distance = unattrib(st_length(geometry)))
settfm(add_links, distance = sp_distance / aere)


# Shortest Paths ----------------------------------------------------------
if(!net |> activate("edges") |> tidygraph::as_tibble() |> select(from, to) |> 
   atomic_elem() |> all.equal(atomic_elem(select(edges, from, to)))) stop("Mismatch") # Check

distances <- st_network_cost(net, weights = edges$distance) # distance in m
n_links <- igraph::distances(net)
range(n_links)

# Checks 
identical(st_geometry(net, "nodes"), nodes$geometry)
dist(unattrib(st_coordinates(nodes$geometry)), unattrib(qM(dist_ttime_mats$sources)))
sp_distances <- st_distance(nodes)

# Network Route Efficiency
nre <- mean(sp_distances / distances, na.rm = TRUE)
rnre <- mean(sp_distances / dist_ttime_mats$distances, na.rm = TRUE) # Real NRE

# Now adding edges
identical(st_geometry(net, "edges"), edges$geometry)
net_ext_data <- rbind(select(edges, distance, geometry), select(add_links, distance, geometry))
net_ext <- as_sfnetwork(net_ext_data, directed = FALSE)
plot(net_ext)
identical(st_geometry(net_ext, "nodes"), nodes$geometry) # Not the case, thus need to recalculate spherical distance as well
ind_ext <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext, "nodes"))))
sp_distances_ext <- st_distance(st_geometry(net_ext, "nodes"))[ind_ext, ind_ext]
identical(sp_distances_ext, sp_distances)

distances_ext <- st_network_cost(net_ext, weights = "distance")[ind_ext, ind_ext]
all.equal(sp_distances_ext, sp_distances)
nre_ext <- mean(sp_distances_ext / distances_ext, na.rm = TRUE)

nre_ext / nre
rnre * (nre_ext / nre) # Reported increase
mean(distances / distances_ext, na.rm = TRUE) 

# Per link gain in NRE: takes a few mins
add_links$nre_per_link <- sapply(seq_row(add_links), function(i) {
  net_extd = as_sfnetwork(rbind(select(edges, distance), 
                                subset(add_links, i, distance)), directed = FALSE)
  distances_extd = st_network_cost(net_extd, weights = "distance")
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_extd, "nodes"))))
  mean(sp_distances / distances_extd[ind, ind], na.rm = TRUE)
})
# Percent increase
add_links$nre_gain_perc <- (unattrib(add_links$nre_per_link / nre) - 1) * 100
descr(add_links$nre_gain_perc)

# Plot percent increase
pdf("figures/PE/trans_CEMAC_network_NRE_gain_perc.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "nre_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),   
           col.legend = tm_legend(expression(Delta~"%"~"NRE"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()


# Gravity weighted versions
gravity <- replace_inf(tcrossprod(nodes$population) / sp_distances, set = TRUE) |> unclass()
nre_wtd <- fmean(unattrib(sp_distances / distances), w = gravity)
rnre_wtd <- fmean(unattrib(sp_distances / dist_ttime_mats$distances), w = gravity)
all.equal(sp_distances_ext, sp_distances)
nre_ext_wtd <- fmean(unattrib(sp_distances_ext / distances_ext), w = gravity)

rnre_wtd
nre_ext_wtd / nre_wtd
rnre_wtd * (nre_ext_wtd / nre_wtd) # Reported increase

# Per link gain in NRE: weighted: takes a few mins
add_links$nre_wtd_per_link <- sapply(seq_row(add_links), function(i) {
  net_extd = as_sfnetwork(rbind(select(edges, distance), 
                                subset(add_links, i, distance)), directed = FALSE)
  distances_extd = st_network_cost(net_extd, weights = "distance")
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_extd, "nodes"))))
  fmean(unattrib(sp_distances / distances_extd[ind, ind]), w = gravity)
})
# Percent increase
add_links$nre_wtd_gain_perc <- (unattrib(add_links$nre_wtd_per_link / nre_wtd) - 1) * 100
descr(add_links$nre_wtd_gain_perc)

# Plot percent increase
# <Figure 19: RHS>
pdf("figures/PE/trans_CEMAC_network_NRE_wtd_gain_perc.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "nre_wtd_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),   
           col.legend = tm_legend(expression(Delta~"%"~"NRE WTD"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()


# (Optional) loading saved results ------------------------------------------------------------

qreadm("results/transport_network/PE/PE_results_MACR_90kmh_google.qs")

# Market Access -------------------------------------------------------------------------------

# Matching to IWI, and choosing distance-weighted nearest cell value within 10km
IWI <- qread("/Users/sebastiankrantz/Documents/IFW Kiel/Africa-Infrastructure/data/intermediate/other/IWI.qs")

IWI_grid <- IWI |>
  # subset(countrycode::countryname(country_name, "iso3c") %iin% c(CEMAC, "NGA", "COD")) |>
  transform(round_to_kms_fast(lon, lat, 10)) |>
  group_by(lon, lat) |>
  select(IWI) |> fmean() |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) 

# Only in CEMAC
IWI_grid %<>% subset(vlengths(st_within(IWI_grid, base::subset(GADM0_africa, GID_0 %in% CEMAC))) > 0)

pdf("figures/IWI_CEMAC_GRID.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(IWI_grid) + 
  tm_symbols(fill = "IWI", size = 0.235, shape = 15, fill_alpha = 0.6, # ?points
             fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
             fill.legend = tm_legend("IWI", position = c("right", "bottom"), frame = FALSE, 
                                     text.size = 1.4, title.size = 1.4)) + 
  tm_layout(frame = FALSE)
dev.off()

outcomes_coords <- IWI |> select(lon, lat) |> qM()
nodes_coord_mat <- st_coordinates(nodes) # graph_nodes |> select(lon, lat) |> qM() 
nodes$N_IWI <- nodes$IWI <- NA_real_
for (i in seq_row(nodes)) {
  d = geodist::geodist(nodes_coord_mat[i, , drop = FALSE], outcomes_coords, measure = "haversine") 
  ind = which(d < 10e3)
  w = 1 / d[ind]
  nodes$IWI[i] = fmean.default(IWI$IWI[ind], w = w)
  nodes$N_IWI[i] = length(ind)
}
settfm(nodes, wealth = IWI * population)
rm(d, ind, w, nodes_coord_mat, outcomes_coords, IWI); gc()
fndistinct(atomic_elem(nodes))

pdf("figures/trans_CEMAC_network_IWI.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(nodes) + 
  tm_dots(fill = "IWI", size = 0.4,
          fill.scale = tm_scale_intervals(7, style = "pretty", values = "inferno"),
          fill.legend = tm_legend("Productivity (IWI)", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.6)) + 
  tm_layout(frame = FALSE)
dev.off()

# Creating GDP variable based on wealth index 
CEMAC_GDP <- africamonitor::am_data(ctry = CEMAC, series = c("NY_GDP_MKTP_KD", "SP_POP_TOTL"), 
                                    expand.date = TRUE, gen = "Year", keep.date = FALSE) |>
             group_by(year = Year) |> 
             summarise(gdp = fsum(NY_GDP_MKTP_KD), 
                       population = fsum(SP_POP_TOTL))
# GDP growth
G(CEMAC_GDP, t = ~year) |> subset(year >= 2000) |> fmedian()
# 80% of CEMAC city population
fsum(nodes$population) / fsum(AP24_CEMAC$population)
# 45% of total CEMAC population
fsum(nodes$population) / CEMAC_GDP[year == max(year), population]
# Now scaling
nodes$gdp_other = proportions(nodes$wealth) * 0.8 * 90889016753 # CEMAC_GDP[year == max(year), gdp] # Assume 80% of GDP (big cities are more productive)

pdf("figures/trans_CEMAC_network_GDP_other.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(mutate(nodes, gdp_b = gdp_other / 1e9)) + 
  tm_dots(fill = "gdp_b", size = 0.4,
          fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
          fill.legend = tm_legend("Other Prod. ($'15B)", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.6)) + 
  tm_layout(frame = FALSE)
dev.off()

# Using Crop Production Instead ----------------------------------------------------

# load("data/transport_network/trans_CEMAC_network_param_google.RData")
# nodes <- nodes_param
# edges <- edges_param
# add_links <- add_links_param
# rm(nodes_param, edges_param, add_links_param)

SPAM <- fread("/Users/sebastiankrantz/Documents/Data/SPAM2020/Global_CSV/spam2020V1r0_global_production/spam2020V1r0_global_P_TA.csv")
SPAM[["TOTAL"]] <- SPAM |> select(BANA_A:YAMS_A) |> psum()
fndistinct(SPAM)
SPAM_CEMAC <- SPAM |> subset(countrycode::countryname(ADM0_NAME, "iso3c") %iin% CEMAC)
SPAM %<>% select(-(BANA_A:YAMS_A))
# Focus on top 11 crops (80% of production) -> Improve crop production through transportation
TOP80_crops <- SPAM_CEMAC |> select(BANA_A:YAMS_A) |> fsum() |> proportions() |> 
  sort(decreasing = TRUE) |> cumsum() |> is_less_than(0.8) |> which() |> names()
SPAM_CEMAC[["TOP80"]] <- SPAM_CEMAC |> get_vars(TOP80_crops) |> psum()
# SPAM_CEMAC %<>% select(-(BANA_A:YAMS_A))
gc()


tfm(edges_real) <- atomic_elem(edges)
for (v in c(TOP80_crops, "TOTAL", "TOP80")) {
  print(v)
  Sys.sleep(1)
  pdf(sprintf("figures/crops/SPAM_CEMAC_%s.pdf", v), width = 6.5, height = 8)
  print(tm_basemap("CartoDB.Positron", zoom = 6) +
    tm_shape(edges_real) + tm_lines(col = "grey30") + 
    tm_shape(SPAM_CEMAC |>
      st_as_sf(coords = c("x", "y"), crs = 4326) |>
      tfmv(v, `/`, 1e3)) +
    tm_symbols(fill = v, size = 0.217, shape = 15, fill_alpha = 0.7, # ?points
            fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
            fill.legend = tm_legend("Crop Production (GT)", position = c("right", "bottom"), frame = FALSE, 
                                    text.size = 1.4, title.size = 1.6)) + 
    tm_layout(frame = FALSE))
  dev.off()
}

st_distance(nodes) %>% extract(lower.tri(.)) %>% quantile(c(0, 0.001, 0.01, 0.025))

outcomes_coords <- SPAM_CEMAC |> select(x, y) |> qM()
nodes_coord_mat <- st_coordinates(nodes) # graph_nodes |> select(lon, lat) |> qM() 
tfm(nodes) <- c("N_SPAM", "SPAM_TOTAL", "SPAM_TOP80", paste0("SPAM_", TOP80_crops)) |> 
  sapply(\(x) NA_real_, simplify = FALSE)
for (i in seq_row(nodes)) {
  d = geodist::geodist(nodes_coord_mat[i, , drop = FALSE], outcomes_coords, measure = "haversine") 
  ind = which(d < 50e3)
  w = 1 - d[ind] / 50e3
  nodes$SPAM_TOTAL[i] = if(length(ind)) fsum.default(SPAM_CEMAC$TOTAL[ind], w = w) else 0
  nodes$SPAM_TOP80[i] = if(length(ind)) fsum.default(SPAM_CEMAC$TOP80[ind], w = w) else 0
  for (c in TOP80_crops) {
    nodes[[paste0("SPAM_", c)]][i] = if(length(ind)) fsum.default(SPAM_CEMAC[[c]][ind], w = w) else 0
  }
  nodes$N_SPAM[i] = length(ind)
}
rm(d, ind, w, nodes_coord_mat, outcomes_coords, SPAM); gc()
fndistinct(atomic_elem(nodes))

for (v in c("SPAM_TOTAL", "SPAM_TOP80")) {
  pdf(sprintf("figures/trans_CEMAC_network_%s.pdf", v), width = 6.5, height = 8)
  print(tm_basemap("CartoDB.Positron", zoom = 6) +
    tm_shape(edges_real) + tm_lines(col = "grey30") + 
    tm_shape(tfmv(nodes, v, `/`, 1e3)) + 
    tm_dots(fill = v, size = 0.4,
            fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
            fill.legend = tm_legend("Crop Production (GT)", position = c("right", "bottom"), frame = FALSE, 
                                    text.size = 1.4, title.size = 1.6)) + 
    tm_layout(frame = FALSE))
  dev.off()
}

# Conflict (ACLED) -----------------------------------------------------------------

ACLED <- fread("/Users/sebastiankrantz/Documents/Data/ACLED/Africa_1997-2024_Nov22.csv")
settfm(ACLED, iso3c = countrycode::countryname(country, "iso3c"))
ACLED %<>% subset(iso3c %in% c(CEMAC, "NGA", "COD") & year >= 2018) 

ACLED_grid <- ACLED |>
  transform(round_to_kms_fast(longitude, latitude, 10)) |>
  group_by(lon, lat) |>
  compute(sqrt_fatalities = sqrt(fatalities), keep = "fatalities") |> fsum() |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

pdf("figures/ACLED_CEMAC_FATALITIES.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(ACLED_grid) + 
  tm_symbols(fill = "sqrt_fatalities", size = 0.235, shape = 15, fill_alpha = 0.7, # ?points
             fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
             fill.legend = tm_legend(expression(sum(sqrt("ACLED Fatalities"))), position = c("right", "bottom"), frame = FALSE, 
                                     text.size = 1.4, title.size = 1.4)) + 
  tm_layout(frame = FALSE)
dev.off()

outcomes_coords <- ACLED |> select(longitude, latitude) |> qM()
nodes_coord_mat <- st_coordinates(nodes) # graph_nodes |> select(lon, lat) |> qM() 
nodes$N_ACLED <- nodes$ACLED_fatalities_sqrt <- NA_real_
for (i in seq_row(nodes)) {
  d = geodist::geodist(nodes_coord_mat[i, , drop = FALSE], outcomes_coords, measure = "haversine") 
  ind = which(d < 50e3)
  w = 1 - d[ind] / 50e3
  nodes$ACLED_fatalities_sqrt[i] = if(length(ind)) fsum.default(sqrt(ACLED$fatalities[ind]), w = w) else 0
  nodes$N_ACLED[i] = length(ind)
}
fndistinct(atomic_elem(nodes))

pdf("figures/trans_CEMAC_network_ACLED.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(nodes) + 
  tm_dots(fill = "ACLED_fatalities_sqrt", size = 0.4,
          fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
          fill.legend = tm_legend(expression(sum(sqrt("ACLED Fatalities"))), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.6)) + 
  tm_layout(frame = FALSE)
dev.off()

# Putting conflict into edges as cost estimates...
edges_real$N_ACLED <- edges_real$ACLED_fatalities_sqrt <- NA_real_
dmat <- unclass(st_distance(edges_real, st_as_sf(ACLED, coords = c("longitude", "latitude"), crs = 4326)))
for (i in seq_row(edges_real)) {
  ind = which(dmat[i, ] < 50e3)
  w = 1 - dmat[i, ind] / 50e3
  edges_real$ACLED_fatalities_sqrt[i] = if(length(ind)) fsum.default(sqrt(ACLED$fatalities[ind]), w = w) else 0
  edges_real$N_ACLED[i] = length(ind)
}
tfm(edges) <- edges_real |> qDF() |> gvr("ACLED")
fndistinct(atomic_elem(edges))

pdf("figures/trans_CEMAC_network_ACLED_edges.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + 
  tm_lines(col = "ACLED_fatalities_sqrt", lwd = 2,
           col.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
           col.legend = tm_legend(expression(sum(sqrt("ACLED Fatalities"))), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.6)) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

add_links$N_ACLED <- add_links$ACLED_fatalities_sqrt <- NA_real_
dmat <- unclass(st_distance(add_links, st_as_sf(ACLED, coords = c("longitude", "latitude"), crs = 4326)))
for (i in seq_row(add_links)) {
  ind = which(dmat[i, ] < 50e3)
  w = 1 - dmat[i, ind] / 50e3
  add_links$ACLED_fatalities_sqrt[i] = if(length(ind)) fsum.default(sqrt(ACLED$fatalities[ind]), w = w) else 0
  add_links$N_ACLED[i] = length(ind)
}
fndistinct(atomic_elem(add_links))

pdf("figures/trans_CEMAC_network_ACLED_edges_add.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(rowbind(select(edges_real, ACLED_fatalities_sqrt), select(add_links, ACLED_fatalities_sqrt))) + 
  tm_lines(col = "ACLED_fatalities_sqrt", lwd = 2,
           col.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
           col.legend = tm_legend(expression(sum(sqrt("ACLED Fatalities"))), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.6)) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

rm(d, dmat, ind, w, nodes_coord_mat, outcomes_coords, ACLED_grid); gc()

# Also: Think of roads reducing conflict: Conflict as a kind of market access measures..

# Combine Statistics --------------------------------------------------------------

nodes |>
  select(gdp_other, SPAM_TOP80) |>
  descr()

nodes$gdp = nodes$gdp_other + proportions(nodes$SPAM_TOP80) * 0.2 * 90889016753
all.equal(sum(nodes$gdp), 90889016753) # Check

pdf("figures/trans_CEMAC_network_GDP.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + tm_lines(col = "grey30") + 
  tm_shape(mutate(nodes, gdp_b = gdp / 1e9)) + 
  tm_dots(fill = "gdp_b", size = 0.4,
          fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
          fill.legend = tm_legend("GDP ($'15B)", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.6)) + 
  tm_layout(frame = FALSE)
dev.off()

# # MISC NOTES
# # HDI-like measure: GDP should take 75% of the measure, Crop Production 25%
# 
# # For each scenario: These are the critical links. 
# # And also: These are the critical borders: 
# # Cost associated of those
# 
# # -> Additional simulations with CROPs/Conflict and borders
# 
# # -> Create a few slides!!! !!
# # -> Need to start operationalization: use PE to select links (set budget), and then GE to compute the program
# 
# edges_real$from_city <- nodes$city_country[edges_real$from]
# edges_real$to_city <- nodes$city_country[edges_real$to]
# 
# edges_real |>
#   # subset(from_ctry != to_ctry, from_city, to_city) |>
#   mapview::mapview()
# 
# # Which cross-border link are more critical to improve for crop production and GDP
# # 2 types of Program: corss-border cost + phyical impromvement (separate simulations)
# # Also give names of borders ...
# 
# nodes |> qDT() |>
#   select(gdp, SPAM_TOTAL, ACLED_fatalities_sqrt) |>
#   pivot() |>
#   ggplot(aes(fill = value)) + 
#     geom_histogram() +
#     facet_wrap(~variable)

# Computing total market access ----------------------------------------------------

(MA_real <- total_MA(dist_ttime_mats$distances, nodes$gdp))
(MA <- total_MA(distances, nodes$gdp)) # distances^3.8

# Total gain
(MA_ext <- total_MA(distances_ext, nodes$gdp)) # distances^3.8 

MA_ext / MA

# Needed for later
ma_gain_per_km <- (MA_ext - MA) * 1000

# Compute change in MA from each link
add_links$MA_per_link <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, distance), 
                            subset(add_links, i, distance)), directed = FALSE)
  inv_dist = 1 / unclass(st_network_cost(nete, weights = "distance"))
  diag(inv_dist) = 0
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  sum(inv_dist %*% nodes$gdp[ind])
})
# Percent increase
add_links$MA_gain_perc <- (add_links$MA_per_link / MA - 1) * 100
descr(add_links$MA_gain_perc)

# <Figure 20: LHS> (Use distances^3.8 above to generate RHS)
pdf("figures/PE/trans_CEMAC_network_MACR_gain_perc.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Adding trade costs ----------------------------------------------------------------------

border_time <- fread("data/QSE/model_border_time_mat.csv") |> qM(1)
border_time_sym <- (border_time + t(border_time)) / 2
border_dist <- fread("data/QSE/model_border_dist_mat.csv") |> qM(1)
border_dist_sym <- (border_dist + t(border_dist)) / 2
border_dist_transit <- fread("data/QSE/model_border_dist_mat_transit.csv") |> qM(1)
border_time_transit <- fread("data/QSE/model_border_time_mat_transit.csv") |> qM(1)

# Adding Country Classification
# GADM0_africa <- st_read("/Users/sebastiankrantz/Documents/Data/GADM/gadm_410-levels.gpkg", layer = "ADM_0") %>%
#   subset(GID_0 %in% africamonitor::am_countries$ISO3) %>% st_make_valid()
# GADM0_africa %<>% rmapshaper::ms_simplify(keep = 0.2) %>% st_make_valid()
# mapview::mapView(GADM0_africa)
ctry <- st_within(nodes, GADM0_africa)
table(vlengths(ctry))
ctry[vlengths(ctry) == 0] <- NA
ctry <- as.integer(ctry)
# mapview::mapview(nodes[is.na(ctry), "geometry"])
nodes$iso3c <- GADM0_africa$GID_0[ctry]
anyNA(nodes$iso3c)
rm(ctry); gc()

# Adding trade costs (symmetric since graph is undirected)
edges$from_ctry <- nodes$iso3c[edges$from]
edges$to_ctry <- nodes$iso3c[edges$to]
edges$border_dist <- sapply(seq_row(edges), function(i) border_dist_sym[edges$from_ctry[i], edges$to_ctry[i]])
edges$border_time <- sapply(seq_row(edges), function(i) border_time_sym[edges$from_ctry[i], edges$to_ctry[i]])

# Checks: pretty large costs, but less for distance...
edges |> qDT() |> extract(border_dist > 0, distance / border_dist) |> descr()
edges |> qDT() |> extract(border_time > 0, duration / border_time) |> descr()

# Same for Additional Routes
add_links$from_ctry <- nodes$iso3c[add_links$from]
add_links$to_ctry <- nodes$iso3c[add_links$to]
add_links$border_dist <- sapply(seq_row(add_links), function(i) border_dist_sym[add_links$from_ctry[i], add_links$to_ctry[i]])
add_links$border_time <- sapply(seq_row(add_links), function(i) border_time_sym[add_links$from_ctry[i], add_links$to_ctry[i]])
# mapview::mapview(select(add_links, from_ctry, to_ctry))

# Now: Repeat Market Access Simulations

# Computing total real market access
bdt_nodes <- border_dist_transit[nodes$iso3c, nodes$iso3c]
MA_bc <- total_MA(distances + bdt_nodes, nodes$gdp) # dist_ttime_mats$distances
(MA_bc_real <- total_MA(dist_ttime_mats$distances + bdt_nodes, nodes$gdp))
# Total gain
MA_ext_bc <- total_MA(distances_ext + bdt_nodes, nodes$gdp) 

MA_ext_bc / MA_bc

# Needed for later
ma_gain_per_km_bc <- (MA_ext_bc - MA_bc) * 1000

# Compute change in MA from each link, with border costs
add_links$MA_per_link_bc <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, distance), 
                            subset(add_links, i, distance)), directed = FALSE)
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(nete, "nodes"))))
  inv_dist = 1 / (unclass(st_network_cost(nete, weights = "distance"))[ind, ind] + bdt_nodes)
  diag(inv_dist) = 0
  sum(inv_dist %*% nodes$gdp)
})
# Percent increase
add_links$MA_gain_perc_bc <- (add_links$MA_per_link_bc / MA_bc - 1) * 100
descr(add_links$MA_gain_perc_bc)

# <Figure 21: LHS>
pdf("figures/PE/trans_CEMAC_network_MACR_gain_perc_bc.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_perc_bc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.025, 0.1, 0.25, 0.5, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Compute Ratios
settfm(add_links, 
       MA_gain_bc_ratio = perch_to_diff(MA_per_link_bc, MA_gain_perc_bc) / perch_to_diff(MA_per_link, MA_gain_perc), 
       MA_gain_perc_bc_ratio = MA_gain_perc_bc / MA_gain_perc)

add_links |> gvr("ratio") |> descr()
add_links$MA_gain_bc_ratio |> replace_outliers(c(0, 1), "clip", set = TRUE)

# <Figure 21: RHS>
pdf("figures/PE/trans_CEMAC_network_MACR_gain_bc_ratio.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_bc_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)),
           col.legend = tm_legend(expression(Delta~"MA Ratio"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()


# Estimating Network Building Costs -------------------------------------------

# Adding border costs
settfm(add_links, total_dist = distance + border_dist)
settfm(edges, total_dist = distance + border_dist, total_time = duration + border_time)

# This commented code shows how the 3000m buffer around links is computed and applied
add_links_buff_3km <- st_buffer(add_links, as_units(3000, "m"))
edges_buff_3km <- st_buffer(edges_real, as_units(3000, "m"))

# Adding Ruggedness: https://diegopuga.org/data/rugged/
rugg <- terra::rast("/Users/sebastiankrantz/Documents/Data/Ruggedness/tri.txt")
# max(rugg)
add_links$rugg <- exactextractr::exact_extract(rugg, add_links_buff_3km, fun = "mean")
edges$rugg <- exactextractr::exact_extract(rugg, edges_buff_3km, fun = "mean")
# Adding Population (WorldPop 2020 1km2 global)
pop_wpop <- terra::rast("/Users/sebastiankrantz/Documents/Data/WorldPop/africa_pop_2020_1km.tif")
# max(pop_wpop)
add_links$pop_wpop <- exactextractr::exact_extract(pop_wpop, add_links_buff_3km, fun = "sum")
add_links$pop_wpop_km2 <- unattrib(add_links$pop_wpop / (st_area(add_links_buff_3km) / 1e6))
edges$pop_wpop <- exactextractr::exact_extract(pop_wpop, edges_buff_3km, fun = "sum")
edges$pop_wpop_km2 <- unattrib(edges$pop_wpop / (st_area(edges_buff_3km) / 1e6))
# Cleanup
rm(rugg, add_links_buff_3km, edges_buff_3km, pop_wpop); gc()
all.equal(select(qDF(edges_real), from, to), select(qDF(edges), from, to))
tfm(edges_real) <- atomic_elem(edges)
  
# Plot Ruggedness
# <Figure 23: LHS>
pdf("figures/trans_CEMAC_network_rugg.pdf", width = 6.5, height = 8)
tm_basemap("Esri.WorldTopoMap", zoom = 6) +
  tm_shape(mutate(rbind(select(edges_real, rugg), select(add_links, rugg)), rugg = rugg / 1000)) +
  tm_lines(col = "rugg",
           col.scale = tm_scale_continuous_log1p(7, values = "turbo"),
           col.legend = tm_legend("Ruggedness", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

add_links |> gvr("_km2") |> descr()
edges |> gvr("_km2") |> descr()

# Plot Population Density
# <Figure 23: RHS>
pdf("figures/trans_CEMAC_network_pop_wpop_km2.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(rbind(select(edges_real, pop_wpop_km2), select(add_links, pop_wpop_km2))) +
  tm_lines(col = "pop_wpop_km2",
           col.scale = tm_scale_continuous_log1p(7, values = "turbo"),
           col.legend = tm_legend(expression("Population/km"^2), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# # Plot conflict (ACLED) (OLD): Theo requested
# ACLED <- fread("/Users/sebastiankrantz/Documents/Data/ACLED/Africa_1997-2024_Nov22.csv")
# settfm(ACLED, iso3c = countrycode::countryname(country, "iso3c"))
# ACLED %<>% subset(iso3c %in% CEMAC & year >= 2018) %>% 
#   st_as_sf(coords = c("longitude", "latitude"), crs = 4326)
# 
# pdf("figures/trans_CEMAC_network_ACLED.pdf", width = 6.5, height = 8)
# tm_basemap("CartoDB.Positron", zoom = 6) +
#   tm_shape(rbind(select(edges_real, pop_wpop_km2), 
#                  select(add_links, pop_wpop_km2))) + tm_lines(col = "grey30") + 
#   tm_shape(ACLED) + 
#   tm_dots(fill = "fatalities", size = 0.2,
#           fill.scale = tm_scale_intervals(7, style = "fisher", values = "inferno"),
#           fill.legend = tm_legend("ACLED Fatalities", position = c("right", "bottom"), frame = FALSE, 
#                                   text.size = 1.4, title.size = 1.7)) + 
#   tm_shape(nodes) + tm_dots(size = 0.1, fill = "grey70") +
#   tm_layout(frame = FALSE)
# dev.off()
# 
# rm(ACLED); gc()

descr(edges$ACLED_fatalities_sqrt)
fmedian(log(edges$ACLED_fatalities_sqrt[edges$ACLED_fatalities_sqrt > 0.1]))
hist(log(edges$ACLED_fatalities_sqrt), breaks = 100)

# --> Median approximately at 3, and collier coefficient is at 0.3, thus multiply log by 0.1

# Estimating Road Costs following my Optimal Roads Paper: ----------------------------

# Load Project data from the region
CDB <- readxl::read_xlsx("data/other_inputs/CEMAC_Road_Costs/Road investment costs data IAWT3&4_April18.xlsx") |>
       janitor::clean_names()
CDB %<>% 
  rename(type_of_works = type_of_works_new_construction_rehabilitation_reconstruction_periodic_maintenance_reseals_patching_resurfacing) %>%
  mutate(length = plast(initial_length_km, as.double(final_length_km_if_different)), 
         type_of_works = trimws(tolower(stringi::stri_trans_general(type_of_works, "latin-ascii"))) |> rm_stub("new "), 
         cost_per_km = final_cost_of_works_us_million / length)
descr(CDB$length)
descr(CDB$type_of_works)
descr(CDB$cost_per_km)
CDB %<>% subset(type_of_works %in% c("construction", "rehabilitation", "reconstruction") & number_of_lanes_per_direction == 1)
descr(CDB, cost_per_km ~ type_of_works) |> print(digits = 3) # Pretty much the same...

# Estimates from Collier, Kirchberger & Söderbom (2016)
# Pop Coef
mean(c(0.11, 0.088, 0.082,   # Table 4
       0.077, 0.083, 0.074)) # Taable 5

# Calibrating cost eauation to match median 2L Highway construction cost in Africa (611 million/km)
mean(with(edges, exp(log(120e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1)))) / 1000 

# Applying Equation
settfm(edges, cost_km = exp(log(120e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1))) # 0.00085 * pop_wpop_km2
descr(edges$cost_km)

settfm(add_links, cost_km = exp(log(120e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1))) # 0.00085 * pop_wpop_km2
descr(add_links$cost_km)

descr(rbind(select(edges, cost_km), select(add_links, cost_km)))

# Total Network Length and Cost
sum(edges$distance / 1000) / 1e3
sum(edges$cost_km * edges$distance / 1000) / 1e9
sum(add_links$distance / 1000) / 1e3
sum(add_links$cost_km * add_links$distance / 1000) / 1e9

edges_real$cost_km = edges$cost_km

# Plots
# <Figure A12: LHS>
pdf("figures/trans_CEMAC_network_cost_km_conflict_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(rbind(select(edges_real, cost_km), select(add_links, cost_km))) +
  tm_lines(col = "cost_km",
           col.scale = tm_scale_continuous(7, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# Cost-Benefit Analysis: Adding New Links -------------------------------------------

# Total Cost-Benefit Ratios

# No Frictions
# CEMACGDP22 <- 87520769488 # CEMAC GDP 2022 in constant 2015 USD (3.11% of African GDP in 2022)
ma_gain_per_km / 1e9 # MA gain in billions
ma_gain_per_km / sum(with(add_links, cost_km * distance / 1000)) # MA gain per investment
# With Frictions
ma_gain_per_km_bc / 1e9 # MA gain in billions
ma_gain_per_km_bc / sum(with(add_links, cost_km * distance / 1000)) # MA gain per investment

# MA Gain per Dollar
settfm(add_links, MA_gain_pusd = perch_to_diff(MA_per_link, MA_gain_perc) * 1000 / (cost_km * distance / 1000)) # * 1216
descr(add_links$MA_gain_pusd)
proportions(table(add_links$MA_gain_pusd < 1))

# <Figure 24: LHS (Top)>
pdf("figures/PE/trans_CEMAC_network_MACR_gain_pusd_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Under Frictions: Static
settfm(add_links, MA_gain_pusd_bc = perch_to_diff(MA_per_link_bc, MA_gain_perc_bc) * 1000 / (cost_km * distance / 1000)) # * 1216
descr(add_links$MA_gain_pusd_bc)

# <Figure 24: RHS (Top)>
pdf("figures/PE/trans_CEMAC_network_MACR_gain_pusd_bc_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(add_links) + 
  tm_lines(col = "MA_gain_pusd_bc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Consensus Package
settfm(add_links, 
       consensus = MA_gain_pusd > 1 & MA_gain_pusd_bc > 1, # (| MA_gain_pusd_bc_opt > 1),
       MA_gain_pusd_cons = pmean(MA_gain_pusd, MA_gain_pusd_bc)) # , MA_gain_pusd_bc_opt

# <Figure 24: RHS (Bottom)>
pdf("figures/PE/trans_CEMAC_network_MACR_gain_pusd_cons_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(add_links, consensus, MA_gain_pusd_cons)) + 
  tm_lines(col = "MA_gain_pusd_cons", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Consensus Gains
nrow(subset(add_links, consensus)) / nrow(add_links)
subset(add_links, consensus) |> with(sum(cost_km * distance / 1000)) |> divide_by(1e9)

net_ext_cons <- as_sfnetwork(rbind(select(edges, distance, total_dist), 
                                   subset(add_links, consensus, 
                                          distance, total_dist)), directed = FALSE)
plot(net_ext_cons)
ind_ext_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext_cons, "nodes"))))
identical(st_distance(st_geometry(net_ext_cons, "nodes"))[ind_ext_cons, ind_ext_cons], sp_distances)

distances_ext_cons <- st_network_cost(net_ext_cons, weights = "distance")[ind_ext_cons, ind_ext_cons]
# distances_ext_bc_cons <- st_network_cost(net_ext_cons, weights = "total_dist")[ind_ext_cons, ind_ext_cons]
distances_ext_bc_cons <- distances_ext_cons + bdt_nodes

# Here need to choose which matrices to use: with/without internal/added border costs: 
sum(distances_ext_bc_cons) / sum(distances_ext_cons)
mean(distances_ext_bc_cons / distances_ext_cons, na.rm = TRUE)

# Total gain
(MA_ext_cons_bc <- total_MA(distances_ext_bc_cons, nodes$gdp)) 

MA_ext_cons_bc / MA_bc

ma_gain_per_km_cons <- (MA_ext_cons_bc - MA_bc) * 1000

ma_gain_per_km_cons / 1e9 # MA gain in billions
ma_gain_per_km_cons / sum(with(subset(add_links, consensus), cost_km * distance / 1000)) # MA gain per investment


# Now: Improving Existing Links ------------------------------------------------------------------------

settfm(edges, speed_kmh = (distance / 1000) / (duration / 60))
descr(edges$speed_kmh)

# <Figure 25: LHS>
hist(edges$speed_kmh, breaks = 80, xlab = "Average Link Speed in km/h", main = NULL)
dev.copy(pdf, "figures/trans_CEMAC_network_average_link_speed_hist_google.pdf", width = 6.5, height = 8)
dev.off()

# Inspect
# edges |> select(speed_kmh) |> mapview::mapview() 

tfm(edges_real) <- atomic_elem(edges)

# Plot
# <Figure 25: RHS>
pdf("figures/trans_CEMAC_network_average_link_speed_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "speed_kmh", 
           col.scale = tm_scale_continuous(5, values = "turbo"),
           col.legend = tm_legend("Speed in km/h", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Computing Times 
if(!net |> activate("edges") |> tidygraph::as_tibble() |> select(from, to) |> 
   atomic_elem() |> all.equal(atomic_elem(select(edges, from, to)))) stop("Mismatch") # Check

times <- st_network_cost(net, weights = edges$duration)
edges %<>% mutate(speed_kmh_imp = iif(speed_kmh < 90, 90, speed_kmh),
                  duration_imp = duration * speed_kmh / speed_kmh_imp)
descr(edges, cols = .c(duration, duration_imp))
times_imp <- st_network_cost(net, weights = edges$duration_imp)

# Computing total real market access
(MA_real <- total_MA(dist_ttime_mats$durations, nodes$gdp)) # Original: 1748.128 billion USD/min
(MA <- total_MA(times, nodes$gdp)) #  

# Total gain
(MA_imp <- total_MA(times_imp, nodes$gdp))

MA_imp / MA
# Gain from original: 
(MA_imp / MA) * MA_real

# Needed for later
ma_gain_per_min <- MA_imp - MA

edges$MA_90_min_speed <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_90_min_speed_perc <- (edges$MA_90_min_speed / MA - 1) * 100
descr(edges$MA_90_min_speed_perc)

tfm(edges_real) <- atomic_elem(edges)

# <Figure 26: A>
pdf("figures/PE/trans_CEMAC_network_MACR_90_min_speed_perc_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_90_min_speed_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.01, 0.025, 0.1, 0.25, 0.5, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Considering the addition of proposed links under 90km/h or 65km/h assumption
settfm(add_links, 
       duration_90kmh = distance / kmh_to_mmin(90), 
       duration_65kmh = distance / kmh_to_mmin(65))

# Temporary networks as needed
net_ext_tmp <- as_sfnetwork(rbind(select(edges, duration = duration_imp), 
                                  select(add_links, duration = duration_90kmh)), directed = FALSE)
ind_ext_tmp <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_ext_tmp, "nodes"))))
times_ext_tmp <- st_network_cost(net_ext_tmp, weights = "duration")[ind_ext_tmp, ind_ext_tmp]

# Total gain
(MA_tmp <- total_MA(times_ext_tmp, nodes$gdp)) / 1e9

MA_tmp / MA # * MA_real
(MA_tmp - MA) / 1e9
rm(list = ls()[endsWith(ls(), "_tmp")]); gc()

# Computing total Crop Access ----------------------------------------------------

for (c in c("TOTAL", "TOP80", TOP80_crops)) {
  
  cat(" ", c)
  cvar <- paste0("SPAM_", c)
  print(MA_crop <- total_MA_wtd(times, nodes[[cvar]], nodes$population))
  
  edges[[sprintf("MA_%s_90_min_speed", cvar)]] <- sapply(seq_row(edges), function(i) {
    w = copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
    total_MA_wtd(st_network_cost(net, weights = w), nodes[[cvar]], nodes$population)
  })
  
  # Percent increase
  edges[[sprintf("MA_%s_90_min_speed_perc", cvar)]] <- (edges[[sprintf("MA_%s_90_min_speed", cvar)]] / MA_crop - 1) * 100
  descr(edges[[sprintf("MA_%s_90_min_speed_perc", cvar)]])
  
  tfm(edges_real) <- atomic_elem(edges)
  
  Sys.sleep(1)
  pdf(sprintf("figures/PE/CA/trans_CEMAC_network_MA_%s_90_min_speed_perc_google.pdf", cvar), width = 6.5, height = 8)
  print(tm_basemap("CartoDB.Positron", zoom = 6) +
    tm_shape(edges_real) +
    tm_lines(col = sprintf("MA_%s_90_min_speed_perc", cvar), 
             col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.01, 0.025, 0.1, 0.25, 0.5, Inf)),
             col.legend = tm_legend(expression(Delta~"%"~"CA [MTP/min]"), 
                                    position = c("right", "bottom"), frame = FALSE, 
                                    text.size = 1.5, title.size = 1.7), lwd = 2) + 
    tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE))
  dev.off()
}

# Joint Plot
tmap_options(raster.max.cells = 1e7)

TOP80_crops_perc <- SPAM_CEMAC |> select(BANA_A:YAMS_A) |> fsum() |> 
  proportions() |> extract(TOP80_crops) |> multiply_by(100)
sum(TOP80_crops_perc)

pdf("figures/PE/CA/trans_CEMAC_network_MA_SPAM_JOINT_90_min_speed_perc_google.pdf", width = 11, height = 11)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real |> gvr("SPAM_") |> gvr("_perc") |> 
             gvr("TOTAL", invert = TRUE) |> rm_stub("MA_SPAM_") |>
             rm_stub("_90_min_speed_perc", FALSE) |> 
             rename(function(x) paste0(x, " (", round(TOP80_crops_perc[x], 1), "%)"), 
                    cols = TOP80_crops) |> rm_stub("_A", regex = TRUE) |>
             rename(TOP80 = "TOP80 (80%)") |>
             pivot("geometry")) +
  tm_lines(col = "value", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.01, 0.025, 0.1, 0.25, 0.5, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"CA [MTP/min]"), 
                                  position = tm_pos_in(0.56, 0.4), frame = FALSE, 
                                  text.size = 0.55, title.size = 0.75), lwd = 1.5) + 
  tm_facets_wrap("variable", ncols = 4) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE,
            panel.label.bg.color = "white",
            panel.label.frame = FALSE)
dev.off()


# Adding Border Frictions and Repeating ------------------------------------------

btt_nodes <- border_time_transit[nodes$iso3c, nodes$iso3c]

# Computing total real market access
(MA_real_bt <- total_MA(dist_ttime_mats$durations + btt_nodes, nodes$gdp)) 
(MA_bt <- total_MA(times + btt_nodes, nodes$gdp)) / 1e9 # dist_ttime_mats$durations 

MA_bt / MA * MA_real # Reported increase

# Total gain
(MA_imp_bt <- total_MA(times_imp + btt_nodes, nodes$gdp)) / 1e9
MA_imp_bt / MA * MA_real # Reported 
MA_imp_bt / MA_bt 

# Needed for later
ma_gain_per_min_bt <- MA_imp_bt - MA_bt 

# Per link
edges$MA_90_min_speed_bt <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$duration, i, edges$duration_imp, vind1 = TRUE)
  inv_dur = 1 / (unclass(st_network_cost(net, weights = w)) + btt_nodes)
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_90_min_speed_bt_perc <- (edges$MA_90_min_speed_bt / MA_bt - 1) * 100
descr(edges$MA_90_min_speed_bt_perc)

tfm(edges_real) <- atomic_elem(edges)

# <Figure 27: LHS>
pdf("figures/PE/trans_CEMAC_network_MACR_90_min_speed_bt_perc_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_90_min_speed_bt_perc", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.01, 0.025, 0.1, 0.25, 0.5, Inf)),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# Compute Ratio
settfm(edges, 
       MA_90_min_speed_bt_ratio = replace_na(perch_to_diff(MA_90_min_speed_bt, MA_90_min_speed_bt_perc) / perch_to_diff(MA_90_min_speed, MA_90_min_speed_perc), 0), 
       MA_90_min_speed_bt_perc_ratio = MA_90_min_speed_bt_perc / MA_90_min_speed_perc)
descr(edges$MA_90_min_speed_bt_ratio)
edges$MA_90_min_speed_bt_ratio |> replace_outliers(c(0, 1), "clip", set = TRUE)

# <Figure 27: RHS>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_CEMAC_network_MACR_90_min_speed_bt_ratio_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_90_min_speed_bt_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)), # breaks = c(0, 0.25, 0.5, 1, 1.5, 2),  # For Ratio
           col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Adding border costs
settfm(edges, total_time = duration + border_time, total_time_imp = duration_imp + border_time)

# Link Upgrading Costs ---------------------------------------------------

settfm(edges, upgrade_cat = nif(speed_kmh < 60, "Upgrade: <60km/h", speed_kmh >= 60 & speed_kmh < 80, "Mixed Works: 60-80km/h", 
                                speed_kmh >= 80 & speed_kmh < 90, "Resurfacing: 80-90km/h", speed_kmh >= 90, "Nothing: >90km/h") |> 
         factor(levels = c("Nothing: >90km/h", "Resurfacing: 80-90km/h", "Mixed Works: 60-80km/h", "Upgrade: <60km/h")))
table(edges$upgrade_cat)
anyNA(edges$upgrade_cat)

# <Figure 29: LHS>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/trans_CEMAC_network_type_of_work_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "upgrade_cat", 
           col.scale = tm_scale_categorical(values = "turbo"),
           col.legend = tm_legend("Type of Work", 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.2, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Now costing the categories
descr(with(edges, # subset(edges, upgrade_cat == "Resurfacing"), 
           exp(log(28.4e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1)))) 
descr(with(edges, # subset(edges, upgrade_cat == "Mixed Works"), 
           exp(log(64.6e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1)))) 
descr(with(edges, # subset(edges, upgrade_cat == "Upgrade"), 
           exp(log(101.6e3) - 0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1)))) 

edges %<>%
  mutate(ug_cost_km = -0.11 * (distance > 50e3) + 0.12 * log(rugg) + 0.085 * log(pop_wpop_km2+1) + 0.1 * log(ACLED_fatalities_sqrt+1), 
         ug_cost_km = nif(upgrade_cat == "Resurfacing: 80-90km/h", clip5perc(exp(ug_cost_km + log(28.4e3))),
                          upgrade_cat == "Mixed Works: 60-80km/h", clip5perc(exp(ug_cost_km + log(64.6e3))),
                          upgrade_cat == "Upgrade: <60km/h", clip5perc(exp(ug_cost_km + log(101.6e3))), 
                          upgrade_cat == "Nothing: >90km/h", 0)) 

descr(edges$ug_cost_km)
descr(edges, ug_cost_km ~ upgrade_cat)
hist(edges$ug_cost_km, breaks = 80)

# <Figure 29: RHS>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/trans_CEMAC_network_upgrading_costs_conflict_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "ug_cost_km", 
           col.scale = tm_scale_continuous(7, values = "yl_or_rd"), 
           col.legend = tm_legend("USD'15/km", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.4, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Total Costs and Breakdown
sum(edges$ug_cost_km * edges$distance / 1000) / 1e9
fsum(edges$ug_cost_km * edges$distance / 1000, edges$upgrade_cat) / 1e9


# Cost-Benefit Analysis ---------------------------------------------------

# No Frictions
ma_gain_per_min / 1e9 # MA gain in billions
ma_gain_per_min / sum(with(edges, ug_cost_km * distance / 1000)) # MA gain per investment

# With Frictions
ma_gain_per_min_bt / 1e9 # MA gain in billions
ma_gain_per_min_bt / sum(with(edges, ug_cost_km * distance / 1000)) # MA gain per investment

# MA Gain per Dollar
settfm(edges, MA_gain_pusd = perch_to_diff(MA_90_min_speed, MA_90_min_speed_perc) / (ug_cost_km * distance / 1000)) # * 1216
edges$MA_gain_pusd |> replace_inf(set = TRUE)
# edges$MA_gain_pusd |> replace_na(0, set = TRUE)
descr(edges$MA_gain_pusd)
proportions(table(edges$MA_gain_pusd < 1))
proportions(table(edges$MA_gain_pusd < 2))

# <Figure 30: LHS (Top)>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_CEMAC_network_MACR_gain_90_min_speed_pusd_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) + 
  tm_lines(col = "MA_gain_pusd", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.25, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Under Frictions: Static
settfm(edges, MA_gain_pusd_bt = perch_to_diff(MA_90_min_speed_bt, MA_90_min_speed_bt_perc) / (ug_cost_km * distance / 1000)) # * 1216
edges$MA_gain_pusd_bt |> replace_inf(set = TRUE)
# edges$MA_gain_pusd_bt |> replace_na(0, set = TRUE)
descr(edges$MA_gain_pusd_bt)
proportions(table(edges$MA_gain_pusd_bt < 1))
proportions(table(edges$MA_gain_pusd_bt < 2))


# <Figure 30: RHS (Top)>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_CEMAC_network_MACR_gain_90_min_speed_pusd_bt_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_gain_pusd_bt", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.25, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Consensus Package
settfm(edges, 
       consensus = is.finite(MA_gain_pusd) & (MA_gain_pusd > 1 & MA_gain_pusd_bt > 1), #  | MA_gain_pusd_bt_opt > 1
       MA_gain_pusd_cons = pmean(MA_gain_pusd, MA_gain_pusd_bt)) # , MA_gain_pusd_bt_opt

# <Figure 30: RHS (Bottom)>
tfm(edges_real) <- atomic_elem(edges)
pdf("figures/PE/trans_CEMAC_network_MACR_gain_90_min_speed_pusd_cons_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(subset(edges_real, !consensus)) + tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(edges_real, consensus, MA_gain_pusd_cons)) +
  tm_lines(col = "MA_gain_pusd_cons", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"),
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.25, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Consensus Gains
nrow(subset(edges, consensus)) / nrow(edges)
subset(edges, consensus) |> with(sum(ug_cost_km * distance / 1000)) |> divide_by(1e9)

net_imp_cons <- as_sfnetwork(rbind(subset(edges, !consensus, duration, total_time), 
                                   subset(edges, consensus, duration = duration_imp, total_time = total_time_imp)), 
                             directed = FALSE)
# plot(net_imp_cons)
ind_imp_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_cons, "nodes"))))
identical(st_distance(st_geometry(net_imp_cons, "nodes"))[ind_imp_cons, ind_imp_cons], sp_distances)

times_imp_cons <- st_network_cost(net_imp_cons, weights = "duration")[ind_imp_cons, ind_imp_cons]
times_imp_bt_cons <- st_network_cost(net_imp_cons, weights = "total_time")[ind_imp_cons, ind_imp_cons]
# times_imp_bt_cons <- times_imp_cons + bdt_nodes
sum(times_imp_bt_cons) / sum(times_imp_cons)
mean(times_imp_bt_cons / times_imp_cons, na.rm = TRUE)

# Total gain
MA_imp_cons <- total_MA(times_imp_cons, nodes$gdp)
MA_imp_cons / MA

ma_gain_per_min_cons <- MA_imp_cons - MA

ma_gain_per_min_cons / 1e9 # MA gain in billions
ma_gain_per_min_cons / sum(with(subset(edges, consensus), ug_cost_km * distance / 1000)) # MA gain per investment


# Cost-Benefit Analysis: Joint Scenarios ---------------------------------------------------

settfm(add_links, total_time_90kmh = duration_90kmh + border_time, 
                  total_time_65kmh = duration_65kmh + border_time)

# Plot All Costs
all_costs <- rowbind(existing = select(edges_real, cost_km = ug_cost_km, distance, duration, duration_imp, border_time, total_time, total_time_imp), 
                     new = select(add_links, cost_km = cost_km, distance, duration_imp = duration_90kmh, 
                                  border_time, total_time_imp = total_time_90kmh) |> 
                           transform(duration = duration_imp, total_time = total_time_imp), 
                     idcol = "type")

descr(all_costs$cost_km)

# <Figure 31: LHS>
hist(all_costs$cost_km / 1000, breaks = 80, xlab = "Cost per Km in Thousands of 2015 USD", main = NULL)
dev.copy(pdf, "figures/trans_CEMAC_network_all_costs_conflict_hist_google.pdf", width = 6.5, height = 8)
dev.off()

# <Figure 31: RHS>
pdf("figures/trans_CEMAC_network_all_costs_conflict_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(all_costs) +
  tm_lines(col = "cost_km", 
           col.scale = tm_scale_continuous(7, values = "turbo"), 
           col.legend = tm_legend("USD'15/km", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.25, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# Total Costs
sum(all_costs$cost_km * all_costs$distance / 1000) / 1e9

# Need to repeat simulations for new links with MA denominated in time 

# Baseline Scenario
add_links$MA_per_link_90kmh <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration), 
                            subset(add_links, i, duration = duration_90kmh)), directed = FALSE)
  ind = ckmatch(mctl(st_coordinates(st_geometry(nete, "nodes"))), nodes_coord)
  inv_dur = 1 / unclass(st_network_cost(nete, weights = "duration"))
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$gdp[ind])
})
add_links$MA_per_link_90kmh_perc <- (add_links$MA_per_link_90kmh / MA - 1) * 100
descr(add_links$MA_per_link_90kmh_perc)

# Added Cost Scenario
add_links$MA_per_link_90kmh_bt <- sapply(seq_row(add_links), function(i) {
  nete = as_sfnetwork(rbind(select(edges, duration), 
                            subset(add_links, i, duration = duration_90kmh)), directed = FALSE)
  ind = ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(nete, "nodes"))))
  inv_dur = 1 / (unclass(st_network_cost(nete, weights = "duration"))[ind, ind] + btt_nodes)
  diag(inv_dur) = 0
  sum(inv_dur %*% nodes$gdp)
})
add_links$MA_per_link_90kmh_bt_perc <- (add_links$MA_per_link_90kmh_bt / MA_bt - 1) * 100
descr(add_links$MA_per_link_90kmh_bt_perc)

# Computing Cost-Benefit Ratios
settfm(add_links, 
   MA_gain_90kmh_pusd = perch_to_diff(MA_per_link_90kmh, MA_per_link_90kmh_perc) / (cost_km * distance / 1000) |> replace_inf(),
   MA_gain_90kmh_pusd_bt = perch_to_diff(MA_per_link_90kmh_bt, MA_per_link_90kmh_bt_perc) / (cost_km * distance / 1000) |> replace_inf()
)

# Combining Datasets
all_cb_ratios <- rbind(edges_real |> select(MA_gain_pusd, MA_gain_pusd_bt), # MA_gain_pusd_bt_opt
                       add_links |> select(MA_gain_90kmh_pusd, MA_gain_90kmh_pusd_bt) |> rm_stub("90kmh_", regex = TRUE)) # MA_gain_90kmh_pusd_bt_opt
tfm(all_cb_ratios) <- all_costs |> atomic_elem()
settfm(all_cb_ratios, MA_gain_pusd_bt_ratio = MA_gain_pusd_bt / MA_gain_pusd)
descr(all_cb_ratios)

for (v in .c(pusd, pusd_bt)) { # pusd_bt_opt
  print(v)
  # <Figure 32: First 3 Plots>
  pdf(sprintf("figures/PE/trans_CEMAC_network_MACR_gain_all_90kmh_%s_google.pdf", v), width = 6.5, height = 8)
  print(tm_basemap("CartoDB.Positron", zoom = 6) +
    tm_shape(all_cb_ratios) + 
    tm_lines(col = paste0("MA_gain_", v), 
            col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)),
            col.legend = tm_legend(expression(Delta~"MA/USD"),
                                    position = c("right", "bottom"), frame = FALSE, 
                                    text.size = 1.25, title.size = 1.7), lwd = 2) + 
    tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE))
  dev.off()
}; rm(v)

# Ratio: 
pdf("figures/PE/trans_CEMAC_network_MACR_gain_all_90kmh_pusd_bt_ratio_google.pdf", width = 6.5, height = 8)
print(tm_basemap("CartoDB.Positron", zoom = 6) +
        tm_shape(all_cb_ratios) + 
        tm_lines(col = "MA_gain_pusd_bt_ratio", 
                 col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)), # breaks = c(0, 0.25, 0.5, 1, 1.5, 2),  # For Ratio
                 col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                        position = c("right", "bottom"), frame = FALSE, 
                                        text.size = 1.5, title.size = 2), lwd = 2) + 
        tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
        tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
        tm_layout(frame = FALSE))
dev.off()


# Need again without real edges (for total gain evaluation):
all_cb_ratios_se <- rbind(edges |> select(MA_gain_pusd, MA_gain_pusd_bt), # MA_gain_pusd_bt_opt
                          add_links |> select(MA_gain_90kmh_pusd, MA_gain_90kmh_pusd_bt) |> rm_stub("90kmh_", regex = TRUE)) # MA_gain_90kmh_pusd_bt_opt
tfm(all_cb_ratios_se) <- all_costs |> atomic_elem()
descr(all_cb_ratios_se)


for (i in c(0, 0.5, 1, 2)) {
cat("MA Gain greatern than ", i, fill = TRUE)

# Consensus Package
# settfm(all_cb_ratios, 
#        consensus = is.finite(MA_gain_pusd) & (MA_gain_pusd > i & MA_gain_pusd_bt > i), # MA_gain_pusd_bt_opt > i # 1, 2, or 4
#        MA_gain_pusd_cons = pmean(MA_gain_pusd, MA_gain_pusd_bt)) # MA_gain_pusd_bt_opt
settfm(all_cb_ratios, 
       consensus = is.finite(MA_gain_pusd) & MA_gain_pusd > i & type == "existing", # MA_gain_pusd_bt_opt > i # 1, 2, or 4
       MA_gain_pusd_cons = MA_gain_pusd)
tfm(all_cb_ratios_se) <- atomic_elem(all_cb_ratios)

# <Figure 32: Last 3 Plots>
pdf(sprintf("figures/PE/trans_CEMAC_network_MACR_gain_all_90kmh_pusd_cons_nofr_noadd_MAg%g_google.pdf", i), width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(subset(all_cb_ratios, !consensus)) + tm_lines(lwd = 2, col = "grey70") +
  tm_shape(subset(all_cb_ratios, consensus, MA_gain_pusd_cons)) +
  tm_lines(col = "MA_gain_pusd_cons", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)),
           col.legend = tm_legend(expression(Delta~"MA/USD"), position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.25, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 
dev.off()

# Consensus Gains
all_cb_ratios %$% table(type, consensus) |> proportions(1)
all_cb_ratios %$% table(type, consensus, w = distance) |> # addmargins() |> 
  t() |> fsum.matrix(TRA = "%")
sum(subset(all_cb_ratios, consensus)$distance) / sum(all_cb_ratios$distance)

# Cost
subset(all_cb_ratios, consensus) |> with(sum(cost_km * distance / 1000)) |> divide_by(1e9) |> print()

net_imp_cons <- as_sfnetwork(rbind(
  subset(all_cb_ratios_se, consensus, duration = duration_imp, total_time = total_time_imp),
  subset(all_cb_ratios_se, !consensus & type == "existing", duration, total_time)), directed = FALSE)

# plot(net_imp_cons)
ind_imp_cons <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_cons, "nodes"))))
identical(st_distance(st_geometry(net_imp_cons, "nodes"))[ind_imp_cons, ind_imp_cons], sp_distances)

times_imp_cons <- st_network_cost(net_imp_cons, weights = "duration")[ind_imp_cons, ind_imp_cons]
times_imp_bt_cons <- st_network_cost(net_imp_cons, weights = "total_time")[ind_imp_cons, ind_imp_cons]
# Again adjust frictions scenario. Default: cumulative frictions.
# times_imp_bt_cons <- times_imp_cons + btt_nodes
sum(times_imp_bt_cons) / sum(times_imp_cons)
mean(times_imp_bt_cons / times_imp_cons, na.rm = TRUE)

# Total gain: No Frictions
MA_imp_cons <- total_MA(times_imp_cons, nodes$gdp) 
print(MA_imp_cons / MA) 
ma_gain_per_min_cons <- MA_imp_cons - MA 
print(ma_gain_per_min_cons / 1e9) # MA gain in billions
# MA gain per investment:
print(ma_gain_per_min_cons / sum(with(subset(all_cb_ratios, consensus), cost_km * distance / 1000)))

# Total gain: Frictions
MA_imp_bt_cons <- total_MA(times_imp_bt_cons, nodes$gdp) 
print(MA_imp_bt_cons / MA_bt) 
ma_gain_per_min_cons <- MA_imp_bt_cons - MA_bt
print(ma_gain_per_min_cons / 1e9) # MA gain in billions
# MA gain per investment:
print(ma_gain_per_min_cons / sum(with(subset(all_cb_ratios, consensus), cost_km * distance / 1000)))
}

# Optimizing Border Posts ------------------------------------------------------------------------------------------

# Two Scenarios: 50% or 100% Reduction

edges_real |>
  subset(from_ctry != to_ctry) |> mapview::mapview()

# Optimizing Agents
times_bt <- st_network_cost(net, weights = edges$total_time)
(MA_bt_opt <- total_MA(times_bt, nodes$gdp)) / 1e9

# 50% reduction ++++++++++++++++++++++++++++++++++++++++++++++
edges$MA_bt_50perc_red <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$total_time, i, edges$duration[i] + 0.5 * edges$border_time[i], vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_bt_50perc_red_perc <- (edges$MA_bt_50perc_red / MA_bt_opt - 1) * 100
edges$MA_bt_50perc_red_perc[edges$MA_bt_50perc_red_perc < 1e-5] <- NA
descr(edges$MA_bt_50perc_red_perc)

tfm(edges_real) <- atomic_elem(edges)

pdf("figures/PE/trans_CEMAC_network_MACR_bt_50perc_red_perc_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_bt_50perc_red_perc", 
           col.scale = tm_scale_continuous(6, values = "yl_or_rd", trans = "log1p"),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# 90% reduction to 12 hours ++++++++++++++++++++++++++++++++++++++++++++++

edges$border_time |> unique(sort = TRUE) |> extract(-1) |> median() |> divide_by(60)

edges$MA_bt_90percto12h_red <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$total_time, i, edges$duration[i] + 12*60, vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_bt_90percto12h_red_perc <- (edges$MA_bt_90percto12h_red / MA_bt_opt - 1) * 100
edges$MA_bt_90percto12h_red_perc[edges$MA_bt_90percto12h_red_perc < 1e-5] <- NA
descr(edges$MA_bt_90percto12h_red_perc)

tfm(edges_real) <- atomic_elem(edges)

pdf("figures/PE/trans_CEMAC_network_MACR_bt_90percto12h_red_perc_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_bt_90percto12h_red_perc", 
           col.scale = tm_scale_continuous(6, values = "yl_or_rd", trans = "log1p"),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

# 100% reduction ++++++++++++++++++++++++++++++++++++++++++++++

edges$MA_bt_100perc_red <- sapply(seq_row(edges), function(i) {
  w = copyv(edges$total_time, i, edges$duration[i], vind1 = TRUE)
  inv_dur = 1 / unclass(st_network_cost(net, weights = w))
  diag(inv_dur) = 0 
  sum(inv_dur %*% nodes$gdp) 
})
# Percent increase
edges$MA_bt_100perc_red_perc <- (edges$MA_bt_100perc_red / MA_bt_opt - 1) * 100
edges$MA_bt_100perc_red_perc[edges$MA_bt_100perc_red_perc < 1e-5] <- NA
descr(edges$MA_bt_100perc_red_perc)

tfm(edges_real) <- atomic_elem(edges)

pdf("figures/PE/trans_CEMAC_network_MACR_bt_100perc_red_perc_google.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(col = "MA_bt_100perc_red_perc", 
           col.scale = tm_scale_continuous(6, values = "yl_or_rd", trans = "log1p"),
           col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 1.7), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()


# Macroeconomic Cost-Benefit Analysis (Minimal Required Growth Returns) --------------------------------------------

CEMACGDP23 <- 90889016753 # * 0.8 # CEMAC GDP 2023 in constant 2015 USD

my_PDV <- function(gdp = CEMACGDP23, df = 1.1, gr_old_perc = 3.1, gr_new_perc = 4.11, max_years = 30) {
  gr_old = 1 + gr_old_perc / 100
  gr_new = 1 + gr_new_perc / 100
  years = 1:max_years
  FV = gdp * (gr_new^years - gr_old^years)
  PDV = sum(FV / df^years)
  return(PDV)
}

# What minimum new growth required 
inv_PDV <- function(PDV = 40e9, ...) {
  objective <- function(x) abs(PDV - my_PDV(gr_new_perc = x, ...))
  optimize(objective, c(0, 100), tol = .Machine$double.eps)$minimum
}

# Test
inv_PDV(my_PDV())

# +++ This builds the components of <Table 7> +++
packages <- c(
  # "All Links MA > 0.5" = 3.305958e9, 
  # "All Links MA > 1" = 1.856438e9, 
  # "All Links MA > 2" = 0.8424899e9
  "All Works" = 27e9,
  "All Upgrades" = 19e9,
  "All Links MA > 0.5 Avg" = 4.3e9, 
  "All Links MA > 1 NoFR" = 5.4e9, 
  "All Links MA > 1 NFNN" = 4e9,
  "Regional Projects" = 3.4e9
)

calc_rates <- function(x, bgr) {
  gr_new = inv_PDV(x, gr_old_perc = bgr)
  c("Rate" = gr_new, "Growth of Rate" = (gr_new / bgr - 1) * 100)
}

# Optimistic Growth Scenario
sapply(packages, calc_rates, 3.5) |> t() |> round(3)

# Pessimistic Growth Scenario
sapply(packages, calc_rates, 2) |> t() |> round(3)





#############################################
# Saving PE Results
#############################################

nodes_tmp <- nodes |>
  join(compute(cities_ports, city_port = TRUE,
               keep = .c(city_country, port_locode, port_name, port_status, outflows)), 
       on = c("city_country", "city_port"))

list(nodes = nodes_tmp,
     edges = edges, 
     add_links = add_links) |>
  qsave("results/transport_network/PE/PE_results_MACR_90kmh_google.qs")

nodes_tmp |> transform(set_names(mctl(st_coordinates(geometry)), c("lon", "lat"))) |> 
  atomic_elem() |> qDT() |> fwrite("results/transport_network/PE/csv/PE_results_nodes_MACR_90kmh_google.csv")
edges |> atomic_elem() |> qDT() |> fwrite("results/transport_network/PE/csv/PE_results_edges_MACR_90kmh_google.csv")
add_links |> atomic_elem() |> qDT() |> fwrite("results/transport_network/PE/csv/PE_results_add_links_MACR_90kmh_google.csv")

rm(nodes_tmp)

# Also saving Simplified real edges
edges_real <- qread("data/transport_network/edges_real_simplified.qs")
edges_real %<>% join(atomic_elem(qread("results/transport_network/PE/PE_results_MACR_90kmh_google.qs")$edges), on = c("from", "to"), drop = "x")

# mapview::mapview(edges_real, zcol = "MA_gain_pusd_cons")
edges_real %>% qsave("results/transport_network/PE/PE_edges_real_MACR_90kmh_google.qs")

#############################################
# Saving Graphs for OTN (GE Analysis)
#############################################

qsu(edges)

graph_orig <- edges |> qDT() |> 
  select(from, from_ctry, to, to_ctry, sp_distance, distance, duration, speed_kmh, 
         speed_kmh_imp, duration_imp, border_dist, border_time, total_time, total_time_imp, 
         rugg, pop_wpop, pop_wpop_km2, cost_km, upgrade_cat, ug_cost_km)


settfm(add_links, total_time_90kmh = duration_90kmh + border_time, total_time_65kmh = duration_65kmh + border_time)

graph_add <- add_links |> qDT() |> 
  select(from, from_ctry, to, to_ctry, sp_distance, distance, duration_90kmh, 
         duration_65kmh, border_dist, border_time, total_dist, total_time_90kmh, total_time_65kmh,
         rugg, pop_wpop, pop_wpop_km2, cost_km)

anyNA(cities_ports$city_country)
any_duplicated(na_rm(nodes$city_country))
sum(nodes$city_port)

graph_nodes <- nodes |> qDT() |> 
  transform(set_names(mctl(st_coordinates(geometry)), c("lon", "lat"))) |> 
  select(lon, lat, iso3c, city_country, city_port, population,  IWI, gdp, wealth, gdp_other, SPAM_TOTAL, SPAM_TOP80) |> 
  join(compute(cities_ports, city_port = TRUE,
               keep = .c(city_country, port_locode, port_name, port_status, outflows)), 
       on = c("city_country", "city_port")) |> 
  mutate(outflows = replace_na(outflows))

# Consistency Checks
identical(graph_orig$from_ctry, graph_nodes$iso3c[graph_orig$from])
identical(graph_orig$to_ctry, graph_nodes$iso3c[graph_orig$to])
identical(graph_add$from_ctry, graph_nodes$iso3c[graph_add$from])
identical(graph_add$to_ctry, graph_nodes$iso3c[graph_add$to])

# Saving
for (name in .c(graph_orig, graph_add, graph_nodes)) {
  sprintf("data/transport_network/csv/%s_MACR_90kmh_google.csv", name) |> 
    fwrite(x = get(name))
}

# Also Adding the information to the RData file
load("data/transport_network/trans_CEMAC_network_google.RData")

# Load previous saved graphs
graphs <- sapply(.c(graph_orig, graph_add, graph_nodes), function(name)
  sprintf("data/transport_network/csv/%s_MACR_90kmh_google.csv", name) |> fread())
graphs$graph_orig$add <- FALSE
graphs$graph_add$add <- TRUE

# Joining
nodes %<>% transform(qDF(round(st_coordinates(.), 4))) %>% 
  join(tfm(graphs$graph_nodes, X = round(lon, 4), Y = round(lat, 4)), 
       on = c("X", "Y", "population", "city_port"), drop = "x", overid = 0) %>% select(-X, -Y)
edges %<>% join(graphs$graph_orig, on = c("from", "to"), drop = "x", overid = 2)
add_links %<>% join(graphs$graph_add, on = c("from", "to"), drop = "x", overid = 2)

# Check that network aligns with nodes
allv(st_distance(st_geometry(net, "nodes"), st_geometry(nodes), by_element = TRUE), 0)
identical(st_geometry(net, "edges"), st_geometry(edges))

# Add to network
net %<>% activate("nodes") %>% dplyr::mutate(nodes |> atomic_elem() |> qDF())
net %<>% activate("edges") %>% dplyr::mutate(select(edges, -(from:gravity_dur)) |> atomic_elem() |> qDF())

# Save
TAN_env <- new.env()
load("data/transport_network/trans_CEMAC_network_google.RData", envir = TAN_env)
TAN_env$nodes_param <- nodes
TAN_env$edges_param <- edges
TAN_env$add_links_param <- add_links
TAN_env$net_param <- net
save(list = ls(TAN_env), file = "data/transport_network/trans_CEMAC_network_param_MACR_90kmh_google.RData", envir = TAN_env)



#############################################
# Evaluation Regional Road Projects
#############################################

# load("data/transport_network/trans_CEMAC_network_param_google.RData")
# # ne <- new.env()
# # load("data/transport_network/trans_CEMAC_network.RData", envir = ne)
# # identical(edges$from, ne$edges$from)
# # identical(edges$to, ne$edges$to)
# edges$ug_cost_km <- edges_param$ug_cost_km
# edges$duration %/=% 60
# edges_res <- fread("results/transport_network/PE/csv/PE_results_edges.csv")
# edges %<>% join(edges_res, on = c("from", "to", "distance"), drop =TRUE)
# add_links_res <- fread("results/transport_network/PE/csv/PE_results_add_links.csv")
# add_links %<>% join(add_links_res, on = c("from", "to"), drop =TRUE)
 
edges_real <- qread("data/transport_network/edges_real_simplified.qs") |>
  rmapshaper::ms_simplify(keep = 0.1) |> st_make_valid()

list2env(qread("results/transport_network/PE/PE_results_MACR_90kmh_google.qs"), globalenv())
tfm(edges_real) <- atomic_elem(edges)
# add_links <- add_links_param

edges_real$from_city <- edges$from_city <- nodes$city_country[edges$from] |> setv(NA, edges$from)
edges_real$to_city <- edges$to_city <- nodes$city_country[edges$to] |> setv(NA, edges$to)
add_links$from_city <- nodes$city_country[add_links$from] |> setv(NA, add_links$from)
add_links$to_city <- nodes$city_country[add_links$to] |> setv(NA, add_links$to)

edges_real |>
  select(from_city, to_city) |>
  mapview::mapview()

planned_segments <- matrix(c(
  c("m'baiki - CAF", "bangui - CAF"),
  c("ouesso - COG", "pokola - COG"), # Better use the bypass. 
  # c("bossembele - CAF", ""), # ?
  c("bouar - CAF", "baoro - CAF"), # ?
  c("yaloke - CAF", "174"), # ?
  c("boda - CAF", "174"),
  c("yaloke - CAF", "bossembele - CAF"),
  c("160", "yaloke - CAF"),
  c("binon - CAF", "160"),
  c("baoro - CAF", "binon - CAF"),
  c("bossangoa - CAF", "bossembele - CAF"),
  c("bossangoa - CAF", "169"),
  c("165", "169"),
  c("gore - TCD", "165"),
  c("176", "bangui - CAF"),
  c("bossembele - CAF", "176"),
  c("boda - CAF", "m'baiki - CAF"),
  c("126", "kelo - TCD"),
  c("pala - TCD", "126"),
  c("102", "pala - TCD"),
  c("figuil centre - CMR", "102"),
  c("4", "Medounou - GAB"),
  c("garoua-boulai - CMR", "baboua - CAF"),
  c("52", "dolisie - COG"),
  c("ndende - GAB", "52"),
# ), ncol = 2, byrow = TRUE) |> qDF() |> 
#   set_names(c("from_city", "to_city"))
# 
# planned_segments_theo_chat <- matrix(c(
  # Douala -> Bangui
  c("douala - CMR", "edea - CMR"),
  c("edea - CMR", "19"),
  c("19", "eseka - CMR"),
  c("eseka - CMR", "yaounde - CMR"), # c("yaounde - CMR", "mbalmayo - CMR"),
  c("yaounde - CMR", "akonolinga - CMR"),
  c("akonolinga - CMR", "ayos - CMR"),
  c("ayos - CMR", "85"),
  c("85", "bertoua - CMR"), # c("bertoua - CMR", "91"),
  c("bertoua - CMR", "96"),
  c("96", "garoua-boulai - CMR"),
  c("garoua-boulai - CMR", "baboua - CAF"),
  c("baboua - CAF", "125"),
  c("125", "bouar - CAF"), # c("bouar - CAF", "148"),
  c("bouar - CAF", "baoro - CAF"),
  c("baoro - CAF", "binon - CAF"),
  c("binon - CAF", "160"),
  c("160", "yaloke - CAF"),
  c("yaloke - CAF", "bossembele - CAF"),
  c("bossembele - CAF", "176"),
  c("176", "bangui - CAF"),
  
  # Douala -> N'Djamena
  c("loum - CMR", "douala - CMR"),
  c("loum - CMR", "bafoussam - CMR"),
  c("bafoussam - CMR", "foumban - CMR"),
  c("foumban - CMR", "bankim - CMR"),
  c("bankim - CMR", "mayo-darle - CMR"),
  c("mayo-darle - CMR", "tibati - CMR"),
  c("tibati - CMR", "ngaoundal - CMR"),
  c("ngaoundal - CMR", "ngaoundere - CMR"), # c("ngaoundere - CMR", "touboro - CMR"),
  c("ngaoundere - CMR", "88"),
  c("garoua - CMR", "88"),
  c("garoua - CMR", "figuil centre - CMR"),
  c("maroua - CMR", "figuil centre - CMR"),
  c("maroua - CMR", "moutourwa - CMR"),
  c("maroua - CMR", "mora - CMR"),
  c("mora - CMR", "104"),
  c("100", "104"),
  c("100", "ndjamena - TCD"),
  
  # Yaounde -> Libreville (Paved, focus on frictions)
  c("yaounde - CMR", "mbalmayo - CMR"),
  c("42", "mbalmayo - CMR"),
  c("ebolowa i - CMR", "42"),
  c("ebolowa i - CMR", "34"),
  c("34", "bitam - GAB"),
  c("mengomeyen - GNQ", "bitam - GAB"),
  c("mengomeyen - GNQ", "29"),
  c("29", "oyem - GAB"),
  c("48", "oyem - GAB"),
  c("28", "48"),
  c("Medounou - GAB", "29"),
  c("4", "Medounou - GAB"),
  c("libreville - GAB", "4")

  # 3.⁠ ⁠⁠Lebamba-Mbigou : 84Km;
  # 4.⁠ ⁠⁠Mbigou-Malinga-Mollo : 112 km
  # -> Unimportant for general navigation!!
  
), ncol = 2, byrow = TRUE) |> qDF() |> 
  set_names(c("from_city", "to_city"))


edges$planned <- select(c(edges), from_city, to_city) %in% planned_segments

edges_real |>
  subset(ckmatch(planned_segments, list(from_city, to_city))) |>
  mapview::mapview()

mapview::mapview(add_links)

# add_segments <- c("m'baiki - CAF", "ouesso - COG")

costs_million <- c(CD13_P2 = 994.572, KPLB = 110.008, KMA = 426.858,
                   CBB = 79.28, ND = 290.264)
sum(costs_million)


# Plotting gains: 
edges_real_target <- edges_real |> 
  subset(ckmatch(planned_segments, list(from_city, to_city))) |>
  mutate(duration = duration_imp, total_time = total_time_imp)

edges_real_target <- edges_real_target |> rowbind(fill = TRUE,
  subset(add_links, from_city == "m'baiki - CAF" & to_city == "ouesso - COG", 
         MA_90_min_speed = MA_per_link_90kmh, #MA_gain_perc, 
         MA_90_min_speed_bt = MA_per_link_90kmh_bt, 
         MA_90_min_speed_perc = MA_per_link_90kmh_perc, #MA_gain_perc, 
         MA_90_min_speed_bt_perc = MA_per_link_90kmh_bt_perc, 
         MA_gain_pusd = MA_gain_90kmh_pusd,
         MA_gain_pusd_bt = MA_gain_90kmh_pusd_bt,
         ug_cost_km = cost_km,
         distance = distance,
         duration = duration_90kmh,
         total_time = total_time_90kmh,
         geometry))

edges_real_target |>
  # subset(ckmatch(planned_segments, list(from_city, to_city))) |>
  with(ug_cost_km * distance / 1000) |>
  sum() |> divide_by(1e6) # |> multiply_by(1.19)

for (v in .c(perc, bt_perc)) {
  print(v)
  Sys.sleep(1)
  pdf(sprintf("figures/PE/trans_CEMAC_network_MACR_90_min_speed_%s_planned_projects.pdf", v), width = 6.5, height = 8)
  tm_basemap("CartoDB.Positron", zoom = 6) +
    tm_shape(edges_real) +
    tm_lines(lwd = 2, col = "grey70") +
    tm_shape(edges_real_target) + 
    tm_lines(col = paste0("MA_90_min_speed_", v), 
             col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.01, 0.025, 0.1, 0.25, 0.5, Inf)),
             col.legend = tm_legend(expression(Delta~"%"~"MA [GDP/min]"), position = c("right", "bottom"), frame = FALSE, 
                                    text.size = 1.5, title.size = 1.7), lwd = 2) + 
    tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE) 
  dev.off()
}

pdf("figures/PE/trans_CEMAC_network_MACR_90_min_speed_bt_ratio_planned_projects.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(mutate(edges_real_target, MA_90_min_speed_bt_ratio = replace_na(perch_to_diff(MA_90_min_speed_bt, MA_90_min_speed_bt_perc) / perch_to_diff(MA_90_min_speed, MA_90_min_speed_perc), 0))) +
  tm_lines(col = "MA_90_min_speed_bt_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)), # breaks = c(0, 0.25, 0.5, 1, 1.5, 2),  # For Ratio
           col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()

for (v in .c(pusd, pusd_bt)) {
  print(v)
  Sys.sleep(1)
  pdf(sprintf("figures/PE/trans_CEMAC_network_MACR_gain_all_90kmh_%s_planned_projects.pdf", v), width = 6.5, height = 8)
  print(tm_basemap("CartoDB.Positron", zoom = 6) +
    tm_shape(edges_real) +
    tm_lines(lwd = 2, col = "grey70") +
    tm_shape(edges_real_target) + 
    tm_lines(col = paste0("MA_gain_", v), 
             col.scale = tm_scale_intervals(values = "turbo", breaks = c(0, 0.1, 0.2, 0.5, 1, 2, 5, Inf)),
             col.legend = tm_legend(expression(Delta~"MA/USD"),
                                    position = c("right", "bottom"), frame = FALSE, 
                                    text.size = 1.25, title.size = 1.7), lwd = 2) + 
    tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
    tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
    tm_layout(frame = FALSE))
  dev.off()
} # ; rm(v)

pdf("figures/PE/trans_CEMAC_network_MACR_gain_all_90kmh_pusd_bt_ratio_planned_projects.pdf", width = 6.5, height = 8)
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(edges_real) +
  tm_lines(lwd = 2, col = "grey70") +
  tm_shape(mutate(edges_real_target, MA_90_min_speed_pusd_bt_ratio = MA_gain_pusd_bt / MA_gain_pusd)) +
  tm_lines(col = "MA_90_min_speed_pusd_bt_ratio", 
           col.scale = tm_scale_intervals(values = "turbo", breaks = seq(0, 1, 0.2)), # breaks = c(0, 0.25, 0.5, 1, 1.5, 2),  # For Ratio
           col.legend = tm_legend(expression(Delta~"MA Ratio"), 
                                  position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) + 
  tm_shape(subset(nodes, population > 0)) + tm_dots(size = 0.1) +
  tm_shape(subset(nodes, population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE)
dev.off()



## Total Gains
# nrow(subset(edges, planned)) / nrow(edges)
# subset(edges, planned) |> with(sum(ug_cost_km * distance / 1000)) |> divide_by(1e9)

edges_target <- edges |> subset(ckmatch(planned_segments, list(from_city, to_city))) |>
  mutate(duration = duration_imp, total_time = total_time_imp)

edges_target <- edges_target |> rowbind(fill = TRUE,
  subset(add_links, from_city == "m'baiki - CAF" & to_city == "ouesso - COG", 
         MA_90_min_speed = MA_per_link_90kmh, #MA_gain_perc, 
         MA_90_min_speed_bt = MA_per_link_90kmh_bt, 
         MA_90_min_speed_perc = MA_per_link_90kmh_perc, #MA_gain_perc, 
         MA_90_min_speed_bt_perc = MA_per_link_90kmh_bt_perc, 
         MA_gain_pusd = MA_gain_90kmh_pusd,
         MA_gain_pusd_bt = MA_gain_90kmh_pusd_bt,
         ug_cost_km = cost_km,
         distance = distance,
         duration = duration_90kmh,
         total_time = total_time_90kmh,
         geometry))

# Cost
edges_real_target |> with(ug_cost_km * distance / 1000) |> sum() |> divide_by(1e6) # |> multiply_by(1.19)

net_imp_proj <- as_sfnetwork(rbind(subset(edges, !planned, duration, total_time), 
                                   select(edges_target, duration, total_time)), 
                             directed = FALSE)
# plot(net_imp_proj)
ind_imp_proj <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net_imp_proj, "nodes"))))
identical(st_distance(st_geometry(net_imp_proj, "nodes"))[ind_imp_proj, ind_imp_proj], sp_distances)

times_imp_proj <- st_network_cost(net_imp_proj, weights = "duration")[ind_imp_proj, ind_imp_proj]
times_imp_bt_proj <- st_network_cost(net_imp_proj, weights = "total_time")[ind_imp_proj, ind_imp_proj]
# Again adjust frictions scenario. Default: cumulative frictions.
# times_imp_bt_proj <- times_imp_proj + btt_nodes
sum(times_imp_bt_proj) / sum(times_imp_proj)
mean(times_imp_bt_proj / times_imp_proj, na.rm = TRUE)

if(!exists("MA")) {
  net <- as_sfnetwork(edges, directed = FALSE)
  ind <- ckmatch(nodes_coord, mctl(st_coordinates(st_geometry(net, "nodes"))))
  times <- st_network_cost(net, weights = "duration")[ind, ind]
  MA <- total_MA(times, nodes$gdp)
}

# Total gain: No Frictions
MA_imp_proj <- total_MA(times_imp_proj, nodes$gdp) 
print(MA_imp_proj / MA) 
sum(edges_target$MA_90_min_speed_perc)
ma_gain_per_min_proj <- MA_imp_proj - MA 
print(ma_gain_per_min_proj / 1e9) # MA gain in billions
# MA gain per investment:
print(ma_gain_per_min_proj / sum(with(edges_target, ug_cost_km * distance / 1000)))

# Total gain: Frictions
MA_imp_bt_proj <- total_MA(times_imp_bt_proj, nodes$gdp) 
print(MA_imp_bt_proj / MA_bt) 
sum(edges_target$MA_90_min_speed_bt_perc)
ma_gain_per_min_proj <- MA_imp_bt_proj - MA_bt
print(ma_gain_per_min_proj / 1e9) # MA gain in billions
# MA gain per investment:
print(ma_gain_per_min_proj / sum(with(edges_target, ug_cost_km * distance / 1000)))



# Macroeconomic Cost-Benefit Analysis of Projects ----------------------------

edges_target |>
  subset(ckmatch(planned_segments, list(from_city, to_city))) |>
  with(ug_cost_km * distance / 1000) |>
  sum() |> divide_by(1e6) # |> multiply_by(1.19)


sapply(c(Mycost = 3400257362, Pcost = sum(costs_million)*1e6), calc_rates, 3.5) |> t() |> round(3)

sapply(c(Mycost = 3400257362, Pcost = sum(costs_million)*1e6), calc_rates, 2) |> t() |> round(3)


# Test cost
tmp = fread("data/transport_network/csv/graph_orig_MACR_90kmh_google.csv")
tmp %<>% join(select(edges, from, from_city, to, to_city))

tmp |>
  subset(ckmatch(planned_segments, list(from_city, to_city))) |>
  with(ug_cost_km * distance / 1000) |>
  sum() |> divide_by(1e6) # |> multiply_by(1.19)

