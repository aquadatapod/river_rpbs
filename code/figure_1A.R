# RPBS conceptual figure
# Load Libraries
library(ggplot2)
library(cowplot)
library(patchwork)
library(grid)
library(svglite)

# helper: simple river cross-section for rapids/pool/bench
river_segment <- function(type="rapid", width=1, height=0.4){
  if(type=="rapid"){
    x <- c(0, 0.33, 0.66, 1, 1, 0)
    y <- c(0, 0.15, 0.05, 0.2, -0.2, -0.2)
  } else if(type=="pool"){
    x <- c(0, 0.25, 0.5, 0.75, 1, 1, 0)
    y <- c(0, 0.25, 0.35, 0.25, 0, -0.25, -0.25)
  } else { # benchland
    x <- c(0, 1, 1, 0)
    y <- c(0.1, 0.1, -0.2, -0.2)
  }
  df <- data.frame(x = x*width, y = y*height)
  return(df)
}

# Panel A: geomorphic template (longitudinal sketch)
df_rapid <- river_segment("rapid", width=1.2, height=0.8)
df_pool  <- river_segment("pool", width=1.2, height=0.6)
df_bench <- river_segment("benchland", width=1.2, height=0.45)

pA <- ggplot() +
  geom_polygon(data=df_rapid, aes(x=x, y=y), fill="#a6cee3", color=NA) +
  geom_polygon(data=transform(df_pool, x = x + 1.25), aes(x=x, y=y), 
               fill="#1f78b4", color=NA) +
  geom_polygon(data=transform(df_bench, x = x + 2.5, y = y + 0.05), 
               aes(x=x, y=y), fill="#b2df8a", color=NA) +
  annotate("text", x=c(0.6, 1.85, 3.1), y=c(0.6,0.6,0.45),
           label=c("Rapids\n(coarse, high shear)","Pool\n(fine sediment, 
                   organic\naccumulation)","Benchland\n(terrestrial inputs,\nseasonal wetting)"),
           size=3.4, hjust=0.5) +
  xlim(-0.1,4) + ylim(-0.4,0.9) +
  theme_void() +
  ggtitle("A. RPBS geomorphic template") +
  theme(plot.title=element_text(size=10, face="bold", hjust=0))

# Panel B: Environmental filters
env_df <- data.frame(x=1:3, y=1,
                     label=c("High O2\nCoarse substrate\nHigh shear",
                             "Nutrient sink\nFine sediments\nLow flow",
                             "Detritus inputs\nIntermittent inundation\nMoist microhabitats"))
pB <- ggplot(env_df, aes(x=x, y=y)) +
  geom_tile(aes(width=0.9, height=0.8), fill="#f0f0f0", color=NA) +
  geom_text(aes(label=label), size=3.4) +
  xlim(0.5,3.5) + ylim(0.5,1.5) +
  theme_void() +
  ggtitle("B. Environmental filters") +
  theme(plot.title=element_text(size=10, face="bold", hjust=0.02))

# Panel C: Community assembly (guilds)
guilds <- data.frame(x=c(1,2,3), y=c(1,1,1),
                     label=c("Scrapers\nFilter-feeders\nRheophilic taxa",
                             "Gatherer-collectors\nPredators\nDeposit feeders",
                             "Climbers\nSemi-aquatic taxa\nMulti-habit taxa"))
pC <- ggplot(guilds, aes(x=x,y=y)) +
  geom_rect(aes(xmin=x-0.45, xmax=x+0.45, ymin=0.6, ymax=1.4),
            fill=c("#fee0d2","#fc8d59","#91bfdb"), color=NA) +
  geom_text(aes(label=label), size=3.4) +
  xlim(0.5,3.5) + ylim(0.5,1.6) +
  theme_void() +
  ggtitle("C. Community assembly (traits & guilds)") +
  theme(plot.title=element_text(size=10, face="bold", hjust=0.02))

# Panel D: Emergent patterns (taxonomic & functional diversity)
pD <- ggplot() +
  geom_point(aes(x=c(1,1.6,2.2), y=c(1,1,1)), size=8, 
             color=c("#1b9e77","#d95f02","#7570b3")) +
  annotate("text", x=c(1,1.6,2.2), y=c(1.4,1.4,1.4),
           label=c("Taxonomic\n(Unique taxa)","Functional\n(Diverse guilds)","Indicators\n(Habitat-specific)"),
           size=3.4) +
  xlim(0.5,2.7) + ylim(0.7,1.6) +
  theme_void() +
  ggtitle("D. Emergent biodiversity patterns") +
  theme(plot.title=element_text(size=10, face="bold", hjust=0.02))

# Panel E: Implications (theory & application)
imp_df <- data.frame(x=1, y=1,
                     label="Advances theory: habitat filtering,\ncommunity assembly, and resilience\n\nApplied: improved biomonitoring\n(habitat-specific indicators)\nand management under change")
pE <- ggplot(imp_df, aes(x=1,y=1)) +
  geom_tile(aes(width=1.9, height=1.6), fill="#f7f7f7", color=NA) +
  geom_text(aes(label=label), size=3.5) +
  theme_void() +
  ggtitle("E. Theoretical & applied implications") +
  theme(plot.title=element_text(size=10, face="bold", hjust=0.02))

# Assemble panels using patchwork
final_plot <- (pA + pB) / (pC + pD) / pE + plot_layout(heights = c(1,1,0.8))

# Save as vector PDF and SVG
ggsave("RPBS_conceptual_figure.pdf", final_plot, width=10, height=7, 
       units="in", device = cairo_pdf)
svglite::svglite("RPBS_conceptual_figure.svg", width=9, height=12)
print(final_plot)
svglite::dev.off()
