library(dplyr) 
library(GCCfactor) 
library(kableExtra) 
library(ggplot2)
library(tidyr)
library(psych)

# Φόρτωση δεδομένων
data("UKhouse")
panel <- UKhouse

# Δημιουργία Y_list
Y_list <- panel2list( 
  panel,
  depvar_header = "dlPrice",
  i_header = "Region",
  j_header = "LPA_Type",
  t_header = "Date"
)

# Εκτίμηση μοντέλου
est_multi <- multilevel(Y_list, ic = "BIC3", standarise = TRUE, r_max = 5) 
est_summary <- summary(est_multi)
summary_df <- as.data.frame(est_summary)

# Υπολογισμός μέσου και τυπικής απόκλισης
mean_sd_by_region <- panel %>%
  group_by(Region) %>% 
  summarise(
    Mean = mean(dlPrice, na.rm = TRUE), 
    Std = sd(dlPrice, na.rm = TRUE), 
    .groups = "drop"
  )

group.index = c(10,2,6,3,8,4,5,1,7,9) 
all_regions <- unique(panel$Region)
selected_regions <- all_regions[group.index]

mean_sd_by_region <- mean_sd_by_region %>%
  slice(match(selected_regions, Region)) 

summary_N <- summary_df$N 
summaryN <- summary_N[group.index] 
summary_RIG <- summary_df$"RI global"
summaryRIG <- summary_RIG[group.index] 
summary_RIF <- summary_df$"RI local" 
summaryRIF <- summary_RIF[group.index]
summary_r <- summary_df$r
summaryr <- summary_r[group.index]

summary_df <- data.frame(
  Region = selected_regions,
  N = summaryN,
  Mean = round(mean_sd_by_region$Mean, 3),
  Std = round(mean_sd_by_region$Std, 3),
  r = round(summaryr,3), 
  RIG = round(summaryRIG, 3), 
  RIF = round(summaryRIF, 3) 
)

summary_row <- data.frame(
  Region = "Summary/Average",
  N = sum(summary_df$N),
  Mean = round(mean(mean_sd_by_region$Mean), 3),
  Std = round(mean(mean_sd_by_region$Std), 3),
  r = summary(est_multi)[11,2],
  RIG = round(mean(summaryRIG), 3),
  RIF = round(mean(summaryRIF), 3)
)

final_df <- bind_rows(summary_df, summary_row) 
print(final_df)
# Εκτίμηση μέσων από Gamma
means_vec <- numeric(est_multi$R)
for (j in 1:est_multi$R) {
  means_vec[j] <- mean(est_multi$Gamma[[j]][, 1], na.rm = TRUE)
}
means_matrix <- matrix(round(means_vec, 3), ncol = 1)
est_multiG <- round(est_multi$G,3) 
globalcommon <- est_multiG %*% t(means_matrix)

non_na_indices <- which(!sapply(est_multi$F, anyNA))
clean_lists <- est_multi$F[non_na_indices]
desired_order <- c(7, 2, 5, 6, 3, 4, 1)  
clean_lists <- clean_lists[desired_order]

describe(est_multi$G)
options(warn = -1)
describe(as.data.frame(do.call(cbind, clean_lists)))
options(warn = 0)
cormatF = round(cor(do.call(cbind, clean_lists), use="pairwise.complete.obs"),3) 
cormatG = round(cor(est_multi$G, use="pairwise.complete.obs"),3)
cormatGF = round(cor(cbind(est_multi$G,do.call(cbind, clean_lists)), use="pairwise.complete.obs"),3)
library(ggplot2)
library(zoo)
library(tidyr)
library(dplyr)
library(lemon)

# Δημιουργία ημερομηνιών
quarterly_dates <- seq(as.yearqtr("1996 Q1", format = "%Y Q%q"),
                       as.yearqtr("2021 Q2", format = "%Y Q%q"),
                       by = 0.25)
date_range <- as.Date(range(quarterly_dates))

# Εθνικά δεδομένα
df <- as.data.frame(est_multi$G)
df$time <- as.Date(quarterly_dates)

df_long <- df %>%
  pivot_longer(cols = -time, names_to = "series", values_to = "value")

# Περίοδοι ύφεσης
recession_periods <- data.frame(
  start = as.Date(c("2008-09-01", "2020-03-01")),
  end = as.Date(c("2009-06-01", "2021-03-01")),
  label = c("Financial Crisis", "COVID-19")
)

# --- Figure 2 ---
national_df <- data.frame(
  Date = as.Date(quarterly_dates),
  est_multi$G
)
colnames(national_df)[-1] <- paste0("Component_", 1:(ncol(national_df) - 1))
national_long <- pivot_longer(national_df, cols = -Date, names_to = "Component", values_to = "Value")

p1 <- ggplot(national_long, aes(x = Date, y = Value, color = Component)) +
  geom_rect(data = recession_periods, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(size = 0.8) +
  geom_text(data = recession_periods,
            aes(x = start + (end - start)/2,
                y = max(national_long$Value, na.rm = TRUE) * 0.85,
                label = label),
            inherit.aes = FALSE,
            size = 2.2, fontface = "bold", color = "black", angle = 90, vjust = 0.5) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", limits = date_range) +
  theme_minimal(base_size = 12) +
  labs(title = "Figure 2: Estimated national components", x = "Date", y = "Value")

print(p1)

# --- Regional setup ---
regional_dfs <- lapply(seq_along(clean_lists), function(i) {
  data.frame(
    Date = as.Date(quarterly_dates),
    Value = clean_lists[[i]],
    Region = paste0("Region_", i)
  )
})
regional_long <- do.call(rbind, regional_dfs)

plot_area <- function(regions, title_text) {
  ggplot(filter(regional_long, Region %in% regions),
         aes(x = Date, y = Value, color = Region)) +
    geom_rect(data = recession_periods, inherit.aes = FALSE,
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = "grey90", alpha = 0.6) +
    geom_line(size = 1) +
    geom_text(data = recession_periods,
              aes(x = start + (end - start)/2,
                  y = max(regional_long$Value, na.rm = TRUE) * 0.8,
                  label = label),
              inherit.aes = FALSE,
              size = 2.2, fontface = "bold", color = "black", angle = 90, vjust = 0.5) +
    scale_x_date(date_labels = "%Y", date_breaks = "2 years", limits = date_range) +
    scale_color_brewer(palette = "Set1") +
    theme_minimal(base_size = 12) +
    labs(title = title_text, x = "Date", y = NULL) +
    theme(legend.position = "bottom")
}

p_area1 <- plot_area(c("Region_1", "Region_2", "Region_3"),
                     "Figure 3b: Estimated regional components (Area 1)")
p_area2 <- plot_area(c("Region_4", "Region_5", "Region_6", "Region_7"),
                     "Figure 3c: Estimated regional components (Area 2)")

print(p_area1)
print(p_area2)

# Καθορισμός περιοχών
regional_long$Area <- case_when(
  regional_long$Region %in% c("Region_1", "Region_2", "Region_3") ~ "Area 1",
  regional_long$Region %in% c("Region_4", "Region_5", "Region_6", "Region_7") ~ "Area 2",
  TRUE ~ NA_character_
)

# Μέσος όρος ανά περιοχή
regional_avg <- regional_long %>%
  filter(!is.na(Area)) %>%
  group_by(Date, Area) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# --- Figure 3d ---
p_area12_avg <- ggplot(regional_avg, aes(x = Date, y = Value, color = Area)) +
  geom_rect(data = recession_periods, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(size = 1.2) +
  geom_text(data = recession_periods,
            aes(x = start + (end - start)/2,
                y = max(regional_avg$Value, na.rm = TRUE) * 0.85,
                label = label),
            inherit.aes = FALSE,
            size = 2.2, fontface = "bold", color = "black", angle = 90, vjust = 0.5) +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", limits = date_range) +
  scale_color_manual(values = c("Area 1" = "#1f77b4", "Area 2" = "#ff7f0e")) +
  theme_minimal(base_size = 12) +
  labs(title = "Figure 3d: Common regional components (Area 1 vs Area 2)",
       x = "Date", y = NULL) +
  theme(legend.position = "bottom")

print(p_area12_avg)

# --- Figure 4 ---
p_all_regions <- ggplot(regional_long, aes(x = Date, y = Value)) +
  geom_rect(data = recession_periods, inherit.aes = FALSE,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.6) +
  geom_line(color = "steelblue") +
  geom_text(data = recession_periods,
            aes(x = start + (end - start)/2, y = -3, label = label),
            inherit.aes = FALSE, size = 2.2, fontface = "bold", angle = 90, vjust = 0.5) +
  facet_wrap(~Region, ncol = 2, scales = "free_y") +
  scale_x_date(date_labels = "%Y", date_breaks = "2 years", limits = date_range) +
  theme_minimal(base_size = 11) +
  labs(title = "Figure 4: Estimated regional components for all regions", x = "Date", y = NULL)

print(p_all_regions)


library(reshape2)
library(ggplot2)

plot_group_mean_abs_corr <- function(Y_list, group_names = names(Y_list),
                                     title = "Group Mean Absolute Correlations") {
  n_groups <- length(Y_list)
  M_group  <- matrix(NA, n_groups, n_groups,
                     dimnames = list(group_names, group_names))
  for (i in seq_along(Y_list)) {
    for (j in seq_along(Y_list)) {
      if (i == j) {
        ci <- cor(Y_list[[i]], use = "pairwise.complete.obs")
        M_group[i,i] <- mean(abs(ci), na.rm = TRUE)
      } else if (i < j) {
        cij           <- cor(Y_list[[i]], Y_list[[j]], use = "pairwise.complete.obs")
        M_group[i,j]  <- mean(abs(cij), na.rm = TRUE)
        M_group[j,i]  <- M_group[i,j]
      }
    }
  }
  df_group <- reshape2::melt(M_group)
  colnames(df_group) <- c("Group1","Group2","MeanAbsCorrelation")
  
  p_group_heat <- ggplot(df_group,
                         aes(x = Group1, y = Group2, fill = MeanAbsCorrelation)) +
    geom_tile() +
    scale_y_discrete(limits = rev(levels(df_group$Group2))) +
    scale_fill_gradientn(colors = colorRampPalette(c("white","red"))(10),
                         limits = c(0, 1)) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid  = element_blank(),
      legend.title= element_blank()
    )
  
  print(p_group_heat)
  
}
print(plot_group_mean_abs_corr(Y_list, group_names = names(Y_list),
                               title = "Group Mean Absolute Correlations"))
Y_list_noG <- Y_list
for (i in seq_along(Y_list)) {
  G_common <- est_multi$G %*% t(est_multi$Gamma[[i]])
  Y_list_noG[[i]] <- Y_list[[i]] - G_common
}
Y_list_noGF <- Y_list
for (i in seq_along(Y_list)) {
  Y_list_noGF[[i]] <- est_multi$Resid[[i]]
}
print(plot_group_mean_abs_corr(Y_list_noG, names(Y_list), title = "(b) After Removing National Factors"))

# Μετά την αφαίρεση των εθνικών και τοπικών παραγόντων
print(plot_group_mean_abs_corr(Y_list_noGF, names(Y_list), title = "(c) After Removing National and Local Factors"))
