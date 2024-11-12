###############################################################################
# Identifying Major Transport Routes (for GE analysis of Trans-African Network)
###############################################################################

library(fastverse)
set_collapse(mask = c("manip", "helper", "special"))
fastverse_extend(qs, sf, sfnetworks, tmap, install = TRUE)
fastverse_conflicts()

load("data/transport_network/trans_CEMAC_network.RData")
edges_real <- qread("data/transport_network/edges_real_simplified.qs")
nodes %<>% join(tfm(atomic_elem(cities_ports_rsp_sf), city_port = TRUE), 
                on = c("city_port", "population", "city_country"), drop = "y", attr = TRUE)
# Mussing Acurenam: Was eliminated in network smoothing step
cities_ports_rsp_sf[-na_rm(attr(nodes, "join")$match), ]

# Plot high gravity roads
tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(subset(edges_real, gravity_rd >= 5)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(edges_real, gravity_rd < 5)) +
  tm_lines(col = "grey50", lwd = 2) +
  tm_shape(subset(nodes, city_port)) + tm_dots(size = 0.15) +
  tm_shape(subset(nodes, !city_port & population > 0)) + tm_dots(size = 0.1, fill = "grey20") +
  tm_shape(subset(nodes, !city_port & population <= 0)) + tm_dots(size = 0.1, fill = "grey70") +
  tm_layout(frame = FALSE) 

# 16 largest port-cities
large_cities <- nodes %$% which(population > 2e5 | outflows > 1e6)
length(large_cities)
largest <- nodes$city_country[large_cities]

# Fastest Routes between them
# igraph::all_shortest_paths(net, large_cities[1], large_cities)
large_city_paths <- lapply(large_cities, function(i) 
  st_network_paths(net, from = i, to = large_cities, weights = "duration", mode = "all")) |>
  rowbind()

large_city_paths <- list(nodes = unique(unlist(large_city_paths$node_paths)), 
                         edges = unique(unlist(large_city_paths$edge_paths)))

tm_basemap("CartoDB.Positron", zoom = 6) +
  tm_shape(subset(edges_real, large_city_paths$edges)) +
  tm_lines(col = "gravity_rd",
           col.scale = tm_scale_intervals(values = "inferno", breaks = c(0, 1, 5, 25, 100, Inf)),
           col.legend = tm_legend("Sum of Gravity", position = c("right", "bottom"), frame = FALSE, 
                                  text.size = 1.5, title.size = 2), lwd = 2) +
  tm_shape(subset(edges_real, -large_city_paths$edges)) +
  tm_lines(col = "grey70", lwd = 2) +
  tm_shape(subset(nodes, large_city_paths$nodes)) + tm_dots(size = 0.2) +
  tm_shape(subset(nodes, -large_city_paths$nodes)) + tm_dots(size = 0.1, fill = "grey50") +
  tm_shape(subset(nodes, large_cities)) + tm_dots(size = 0.2, fill = "red") +
  tm_layout(frame = FALSE) 


# Plotting --------------------------------

# Classification
settfm(edges, speed_kmh = (distance / 1000) / (duration / 60^2))
tfm(edges_real) <- qDT(edges) |> select(-geometry)

settfm(nodes, major_city_port = replace_na(population > 2e5 | outflows > 1e6))
sum(nodes$major_city_port)


settfm(nodes, product = nif(major_city_port, NA_integer_, # Heterogeneous products
                            # population > 1e5 & outflows > 1e6, 5L, # Large Port-City: doesn't exist
                            outflows > 0, 3L,       # Port
                            population > 1e5, 4L,   # Large City
                            population > 5e4, 2L,   # Medium-Sized City
                            default = 1L))          # Town/Node
table(nodes$product, na.exclude = FALSE)
setv(nodes$product, whichNA(nodes$product), seq_along(largest) + 4L)
# Need to write this to ensure product classification is available for GE simulation !!!!
nodes_param <- fread("data/transport_network/csv/graph_nodes.csv")
nodes_param |> select(lon, lat) |> qM() |> subtract(st_coordinates(nodes)) |> descr() |> print(digits = 7)
nodes_param$product <- nodes$product
nodes_param |> atomic_elem() |> qDT() |> fwrite("data/transport_network/csv/graph_nodes.csv")
rm(nodes_param)
attr(nodes$product, "levels") <- c("Small Town", "City > 50K", "Port", "City > 100K", paste("Major City", seq_along(largest)))
class(nodes$product) <- "factor"

# Plotting
# <Figure 41: LHS> (Use nname <- "all_routes" above to generate the RHS)
pdf("figures/GE/trans_africa_network_reduced_20_products_EWT.pdf", width = 7.5, height = 9)
tm_basemap("Esri.WorldTopoMap", zoom = 6) +
  tm_shape(edges_real) + # Use edges_real to generate <Figure A21>
  tm_lines(col = "speed_kmh", 
           col.legend = tm_legend("Speed", position = c("right", "bottom"), stack = "h", 
                                  frame = FALSE, text.size = 1.5, title.size = 1.7),
           col.scale = tm_scale_continuous(6, values = "turbo"), # 7, 
           lwd = 2) + 
  tm_shape(droplevels(mutate(nodes, product = fifelse(unclass(product) > 4L, NA, product)))) + 
  tm_dots(fill = "product", size = 0.25, 
          fill.scale = tm_scale_categorical(values = "turbo", value.na = "purple3", label.na = "Own Product"),
          fill.legend = tm_legend("Product", position = c("right", "bottom"), # stack = "h", 
                                  frame = FALSE, text.size = 1.5, title.size = 1.7)) +
  tm_shape(subset(nodes, unclass(product) > 5L)) + tm_dots(size = 0.5, fill = "purple3") +
  tm_layout(frame = FALSE)
dev.off()         

