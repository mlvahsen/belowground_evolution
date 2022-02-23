library(patchwork);library(raster); library(maps); library(cowplot);
library(ggsn); library(ggmap); library(tidyverse);
library(ggrepel)

# Create data frame with GPS points for where cores were taken
locations <- tibble(id = c("Corn Island 1", "Corn Island 2", "Kirkpatrick Marsh 1", "Kirkpatrick Marsh 2"),
                    longitude = c(-76.54357, -76.54361, -76.54811, -76.54945),
                    latitude = c(38.87653, 38.87566, 38.87406, 38.87601))

colors <- rep("blue", 4)

# Create map with scale-bar and N arrow
sbbox <- c(left = -76.555, bottom = 38.87, right = -76.53, top = 38.898)
map <- ggmap(get_stamenmap(sbbox, zoom = 16, maptype = "watercolor")) +
  geom_point(data = locations, aes(x = longitude, y = latitude, color = id), size = 5) +
  scale_color_manual(values = colors) +
  geom_label_repel(data = locations,
                   aes(longitude, latitude,label = id,group = id, color = id),
                   size  = 4, fontface = "bold",
                   box.padding = 0.7, point.padding = 0.5, segment.size = 1.2) +
  scalebar(x.min = -76.555, x.max = -76.532,
           y.min = 38.871,  y.max = 38.872,
           dist = 0.5, transform = TRUE, 
           model = "WGS84", height = 0.3, 
           st.dist = 0.5, dist_unit = "km", st.size = 4, ) +
  xlab("Longitude") +
  ylab("Latitude") + theme_bw() + theme(legend.position = "none")

# Save figure in folder
png(here("figs_tables", "FigS1_map.png"), height = 6, width = 6, res = 300, units = "in")
north2(map, x = 0.28, y = 0.95, scale = 0.06, symbol = 12)
dev.off()