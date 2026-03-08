library(terra)
library(grDevices)
library(dplyr)      
library(FSA)        
library(tidyr)      

# Loading data
input_dir <- "PATH/TO/LCC_TRANSITION_TIF_FILES"
habitats <- c("Broadleaved", "Coniferous", "Arable", "Grassland", "Heathland")

# Continuous color function
col_fun <- colorRamp(c("#734C00", "white", "#01665E"))  # brown → white → green

# Determine max y across all densities for uniform y-axis
max_y <- 0
for(i in 1:5){
 for(j in 1:5){
  f <- paste0("MK_tau_", i, j, ".tif")
  fp <- file.path(input_dir, f)
  if(file.exists(fp)){
   r <- rast(fp)
   vals <- values(r, na.rm=TRUE)
   if(length(vals) > 1){
    d <- density(vals, from=-1, to=1, na.rm=TRUE, n=2048) 
    max_y <- max(max_y, d$y)
   }
  }
 }
}
if(max_y == 0) max_y <- 1

# Set up plotting layout
par(mfrow=c(5,5), mar=c(2,2,2,1), oma=c(4,6,2,2)) 

# Plot each panel
for(i in 1:5){
 for(j in 1:5){
  
  filename <- paste0("MK_tau_", i, j, ".tif")
  filepath <- file.path(input_dir, filename)
  
  if(file.exists(filepath)){
   r <- rast(filepath)
   vals <- values(r, na.rm=TRUE)
   
   if(length(vals) > 1){
    
    mean_val <- mean(vals, na.rm=TRUE)
    d <- density(vals, from=-1, to=1, na.rm=TRUE, n=2048) 
    
    # Empty plot frame
    plot(NA, xlim=c(-1,1), ylim=c(0,max_y), xlab="", ylab="", bty="n", axes=FALSE)
    
    # Axis
    axis(1, cex.axis=1.5, col.axis="grey30", col="grey30")
    axis(2, cex.axis=1.5, col.axis="grey30", col="grey30")
    
    # Draw density in gradient color
    n_seg <- length(d$x) - 1
    for(k in 1:n_seg){
     x_seg <- c(d$x[k], d$x[k+1], d$x[k+1], d$x[k])
     y_seg <- c(0,0,d$y[k+1],d$y[k])
     x_mid <- mean(d$x[k:(k+1)])
     col_rgb <- col_fun((x_mid + 1)/2)
     polygon(x_seg, y_seg,
             col=rgb(col_rgb[1], col_rgb[2], col_rgb[3], maxColorValue=255),
             border=NA)
    }
    
    # Mean line
    abline(v=mean_val, lty=2, lwd=2)
    
   } else {
    # Empty axis if no data
    plot(NA, xlim=c(-1,1), ylim=c(0,max_y), xlab="", ylab="", bty="n", axes=FALSE)
    axis(1, cex.axis=1.8, col.axis="grey30", col="grey30")
    axis(2, cex.axis=1.8, col.axis="grey30", col="grey30")
   }
   
  } else {
   # Missing file → empty axis
   plot(NA, xlim=c(-1,1), ylim=c(0,max_y), xlab="", ylab="", bty="n", axes=FALSE)
   axis(1, cex.axis=1.8, col.axis="grey30", col="grey30")
   axis(2, cex.axis=1.8, col.axis="grey30", col="grey30")
  }
  
  # Column labels (top row)
  if(i==1){
   mtext(habitats[j], side=3, line=0.5, cex=1.3)
  }
  
  # Row labels
  if(j==1){
   mtext(habitats[i], side=2, line=3, cex=1.3, adj=1)
  }
 }
}

############ statistical differences between MK - LCC transitions
# Read rasters and build dataframe
df_list <- list()

for(i in 1:5){
 for(j in 1:5){
  filename <- paste0("MK_tau_", i, j, ".tif")
  filepath <- file.path(input_dir, filename)
  
  if(file.exists(filepath)){
   r <- rast(filepath)
   vals <- values(r, na.rm=TRUE)
   
   if(length(vals) > 0){
    temp_df <- data.frame(vals)  # values(r) returns column b1 by default
    # rename numeric column to MK
    names(temp_df)[1] <- "MK"
    
    temp_df$From <- factor(habitats[i], levels=habitats)
    temp_df$To <- factor(habitats[j], levels=habitats)
    temp_df$Transition <- factor(paste(habitats[i], "→", habitats[j]))
    
    df_list[[length(df_list)+1]] <- temp_df
   }
  }
 }
}

# Combine all transitions
df_list_nonempty <- df_list[sapply(df_list, nrow) > 0]
if(length(df_list_nonempty) > 0){
 df <- bind_rows(df_list_nonempty)
} else stop("No MK data found in any raster")

# Compute mean/median per transition
transition_summary <- df %>%
 group_by(From, To) %>%
 summarise(
  mean_MK = mean(MK, na.rm=TRUE),
  median_MK = median(MK, na.rm=TRUE),
  n = n(),
  .groups="drop"
 )

# Compute significance vs no-change (From = To)
sig_list <- list()

for(h in habitats){
 
 sub <- df %>% filter(From == h)
 
 if(length(unique(sub$To)) > 1){
  
  # Kruskal-Wallis
  kw <- kruskal.test(MK ~ To, data=sub)
  
  if(kw$p.value < 0.05){
   # Dunn post-hoc
   dunn <- dunnTest(MK ~ To, data=sub, method="bh")$res
   ref <- h  # no-change transition
   
   # keep only comparisons involving the reference
   sig <- dunn %>%
    filter(grepl(ref, Comparison)) %>%
    mutate(
     To = ifelse(grepl(paste0(ref," - "), Comparison),
                 gsub(paste0(ref," - "), "", Comparison),
                 gsub(paste0(" - ",ref), "", Comparison)),
     sig = ifelse(P.adj < 0.05, "*", "")
    ) %>%
    select(To, sig)
   
   sig$From <- h
   sig_list[[h]] <- sig
  }
 }
}

sig_df <- bind_rows(sig_list)

# Replace NA with blank for stable transitions
sig_df$sig[is.na(sig_df$sig)] <- ""

# Merge significance with transition summary
results_table <- transition_summary %>%
 left_join(sig_df, by=c("From","To")) %>%
 mutate(
  sig = ifelse(is.na(sig), "", sig),
  MK_mean = round(mean_MK,3),
  MK_label = paste0(MK_mean, sig)
 )

# Pivot to 5×5 transition matrix
transition_matrix <- results_table %>%
 select(From, To, MK_label) %>%
 pivot_wider(names_from = To, values_from = MK_label)

print(transition_matrix)

# Export tables
# Get current working directory
out_dir <- getwd()  

# Export long results table
write.csv(results_table,
          file = file.path(out_dir, "MK_transition_results.csv"),
          row.names = FALSE)

# Export 5x5 transition matrix
write.csv(transition_matrix,
          file = file.path(out_dir, "MK_transition_matrix.csv"),
          row.names = FALSE)

cat("Files exported to:", out_dir, "\n")


# Calculate proportion of LCC
pixel_area_km2 <- 30 * 30 / 1e6   # Landsat pixel area

# Build transition matrix
transition_mat <- matrix(
 0,
 nrow=5,
 ncol=5,
 dimnames=list(habitats, habitats)
)

for(i in 1:5){
 for(j in 1:5){
  
  filename <- paste0("MK_tau_", i, j, ".tif")
  filepath <- file.path(input_dir, filename)
  
  if(file.exists(filepath)){
   
   r <- rast(filepath)
   
   n_pixels <- global(!is.na(r), "sum", na.rm=TRUE)[1,1]
   
   area_km2 <- n_pixels * pixel_area_km2
   
   transition_mat[i,j] <- area_km2
  }
 }
}

# Calculate LCC
# persistence
diag_vals <- diag(transition_mat)

# row totals
row_totals <- rowSums(transition_mat)

# LCC area per class
lcc_area <- row_totals - diag_vals

# TOTAL landscape LCC
total_lcc <- sum(lcc_area)

# proportion relative to TOTAL LCC
lcc_prop <- lcc_area / total_lcc * 100

# Build table
transition_df <- as.data.frame(transition_mat)

transition_df$`Total LCC area (Km²) per land cover type` <- lcc_area
transition_df$`Proportion of LCC (%) per land cover type` <- lcc_prop

transition_df$From <- rownames(transition_mat)

transition_df <- transition_df %>%
 select(
  From,
  all_of(habitats),
  `Total LCC area (Km²) per land cover type`,
  `Proportion of LCC (%) per land cover type`
 )

# rounding for table
transition_df[,-1] <- round(transition_df[,-1],2)

print(transition_df)

# Export
write.csv(
 transition_df,
 file=file.path(input_dir,"Transition_matrix_1990_2023.csv"),
 row.names=FALSE
)

cat("Transition matrix exported\n")


##### Proportion per land cover type
# Function to round a vector so sum == 100
round_to_100 <- function(x, digits=2) {
 x_rounded <- round(x, digits)
 diff <- 100 - sum(x_rounded)
 
 if(diff != 0){
  # Adjust the largest value in the row
  max_idx <- which.max(x_rounded)
  x_rounded[max_idx] <- x_rounded[max_idx] + diff
 }
 
 return(x_rounded)
}

# Compute row-wise proportions
transition_pct <- transition_mat / rowSums(transition_mat) * 100

# Apply rounding adjustment row-wise
transition_pct_adj <- t(apply(transition_pct, 1, round_to_100, digits=2))

# Convert to data.frame
transition_pct_df <- as.data.frame(transition_pct_adj)
transition_pct_df$From <- rownames(transition_mat)

# Reorder columns
transition_pct_df <- transition_pct_df %>%
 select(From, all_of(habitats))

print(transition_pct_df)

# Optional: export
write.csv(
 transition_pct_df,
 file = file.path(input_dir, "Transition_matrix_pct_per_row_100.csv"),
 row.names = FALSE
)


## Proportion (%) of LCC in relation to study area
# total landscape area
total_area <- sum(transition_mat)

# persistence (no change)
persistence_area <- sum(diag(transition_mat))

# total land cover change
lcc_area <- total_area - persistence_area

# proportion of LCC relative to study area
change_percent <- (lcc_area / total_area) * 100

cat("Total study area:", round(total_area,2), "km²\n")
cat("Total LCC area:", round(lcc_area,2), "km²\n")
cat("Percent landscape changed:", round(change_percent,2), "%\n")

############ LCC map
# Binary change rasters list
change_rasters <- list()

for(i in 1:5){
 for(j in 1:5){
  
  filename <- paste0("MK_tau_", i, j, ".tif")
  filepath <- file.path(input_dir, filename)
  
  if(file.exists(filepath)){
   
   r <- rast(filepath)
   
   if(i == j){
    # No change
    change_rasters[[length(change_rasters)+1]] <- r * 0
   } else {
    # Change
    change_rasters[[length(change_rasters)+1]] <- r * 0 + 1
   }
  }
 }
}

# Stack rasters
change_stack <- rast(change_rasters)

# Collapse stack to single binary map
# If any transition pixel = 1 → change
binary_change <- app(change_stack, fun=max, na.rm=TRUE)

#Plot
dev.new(width=10, height=8)

plot(binary_change,
     col=c("grey80","red"),
     legend=FALSE,
     main="Land Cover Change",
     axes=FALSE)

legend("bottomleft",
       legend=c("No change","Change"),
       fill=c("grey80","red"),
       bty="n")

# Export tif
out_file <- file.path(input_dir, "Land_cover_change_map.tif")

# Write raster
writeRaster(binary_change,
            out_file,
            overwrite = TRUE,
            datatype = "INT1U",      
            gdal = c("COMPRESS=LZW"))

