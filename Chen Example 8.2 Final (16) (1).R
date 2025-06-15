library(GrFA)      # (υποθέτουμε ότι περιέχει βοηθητικές συναρτήσεις, π.χ. est_num – εδώ χρησιμοποιούμε δική μας υλοποίηση IC μέσω BaiNg)
library(readxl) # read xlsx
library(dplyr) # use of pipe operator %>% and pull function
library(ggplot2)
library(reshape2)
library(corrplot)
# Ορισμός του working directory (τροποποίησε την διαδρομή κατάλληλα)
graphics.off(); rm(list = ls(all = TRUE)); cat("\014"); getwd()

groupnames <- c("foods", "beverages,alcohol and tobacco", "textiles", "clothing", "furniture", "office supplies",
                "gold and jewelry", "sports and entertainment goods", "cosmetics", "building materials and hardware",
                "medicine and healthcares", "home appliances", "daily necessities", "transportation and communication",
                "books, newspapers and publications", "fuel")

groupnames.short <- c("foods", "beverages", "textiles", "clothing", "furniture", "office",
                      "goldj", "sports", "cosmetics", "building",
                      "medicine", "homeapp", "dailynec", "transportation",
                      "books", "fuel")

seriesnames <- c("NationalIndex","BeiJ","TianJ","HeB","ShanX","NeiMG","LiaoN","JiL","HeLJ","ShangH","JiangS","ZheJ","AnH","FuJ","JiangX","ShanD","HeN","HuB","HuN","GuangD","GuangX","HaiN","ChongQ","SiC","GuiZ","YunN","XiZ","ShaanX","GanS","QingH","NingX","XinJ")

Nm = length(seriesnames)
Mgroups = length(groupnames)

# Ανάγνωση και προεπεξεργασία δεδομένων
raw <- read.csv("RPI3level_m.csv", header = FALSE, stringsAsFactors = FALSE)

Trow = nrow(raw)
Ntotal = ncol(raw)

# group.index = 1:Mgroups #as given in RPI3level_m.csv
group.index = c(2,15,10,4,9,13,1,16,5,7,12,11,6,8,3,14) #table 5 ordering as given in Chen (2023) paper
# group.index =  sample(Mgroups,size=Mgroups,replace = FALSE) #random permutation of index

X <- as.matrix(raw); X_standardized <- scale(X)

# Δημιουργία blocks (ομαδοποίηση κατά κατηγορία)

a1 = seq(1,Ntotal,Nm) ; a2 = seq(0,Ntotal,Nm)[-1]
category_list <- list()
for (i in 1:Mgroups) {category_list[[i]] <- X_standardized[,a1[i]:a2[i]]}

category_list <- category_list[group.index]
groupnames <- groupnames[group.index]
groupnames.short <- groupnames.short[group.index]

# Ορισμός των κριτηρίων
criteria <- c("PC1", "PC2", "PC3", "IC1", "IC2", "IC3", "AIC3", "BIC3", "ER", "GR")

cat("\014")
# Βρόχος για κάθε κατηγορία
for (cat_idx in 1:Mgroups) {
  
  # Εκτύπωση ονόματος κατηγορίας
  cat("\nCategory", cat_idx, ":\n")
  # Ορισμός των εύρων τιμών i (κmax) με βάση τη διάσταση της εκάστοτε κατηγορίας
  i_values <- 1:floor(ncol(category_list[[cat_idx]])/2)
  
  # Δημιουργία πίνακα δεδομένων για την κατηγορία
  results1 <- data.frame(matrix(ncol = length(criteria) + 1, nrow = length(i_values)))
  colnames(results1) <- c("kmax", paste0("est_", criteria))
  results1$kmax <- i_values
  
  # Υπολογισμός των τιμών για κάθε i και αποθήκευση στο data.frame
  for (i in seq_along(i_values)) {
    for (j in seq_along(criteria)) {
      
      # Κλήση της est_num για κάθε κριτήριο
      value <- est_num(category_list[[cat_idx]], kmax = i_values[i], type = criteria[j])
      
      # Αποθήκευση της τιμής, αποφυγή NA σφαλμάτων
      if (!is.null(value) && length(value) == 1) {
        results1[i, j + 1] <- value
      } else {
        results1[i, j + 1] <- NA
      }
    }
  }
  
  # Εκτύπωση του αποτελέσματος
  print(results1)
}

cat("\014")
out1 <- CP(category_list, rmax = 15, r0 = NULL, r = NULL, localfactor = T, type = "IC2")
out2 <- CP(category_list, rmax = 15, r0 = 8, r = rep(7,16), localfactor = T, type = "IC2") #Force 8 global factors
print(out1)
print(out2)

colnames(out2$Ghat)<-paste("gfactor", rep(1:out2$r0hat), sep="")

round( t(out2$Ghat)%*%out2$Ghat , 4 )

# SOME PLOTS
df1 = as.data.frame(out2$Ghat)
df1$Date <- seq(as.Date("2010-01-01"), by = "month", length.out = Trow)
df1m <- melt(df1, id.vars = "Date")

g1 <- ggplot(df1m, aes(x = Date, y = value))+
  geom_line() + facet_wrap(~variable,ncol=2,nrow=4)
g2 <- ggplot(df1m, aes(x = Date, y = value, col=variable))+
  geom_line()
g3 <- ggplot(df1, aes(x=Date, y=gfactor1)) + geom_line()

g1
g2
g3

# global factor correlation
round(cor(out2$Ghat),4)
# Example of cross-correlation between two local factors
round(cor(cbind(out2$Fhat[[1]] , out2$Fhat[[2]] )),4)

# TABLE 5
category.list.global.common <- list()
category.list.local.common <- list()
for (i in 1:Mgroups) {category.list.global.common[[i]] <-out2$Ghat%*%t(out2$loading_G[[i]])}
for (i in 1:Mgroups) {category.list.local.common[[i]] <-out2$Fhat[[i]]%*%t(out2$loading_F[[i]])}

multiple.func <- function(x) { mean( diag( var( x ) ) ) }
Table5row3 <- lapply(out2$residual,multiple.func)
Table5row1 <- lapply(category.list.global.common,multiple.func)
Table5row2 <- lapply(category.list.local.common,multiple.func)
Table5.f = rbind(Table5row1,Table5row2,Table5row3)
colnames(Table5.f) <- 1:Mgroups
rownames(Table5.f) <- c("global","local","idiosyncratic")
Table5.f
groupnames

groupnames.short
colnames(Table5.f) <- groupnames.short
Table5.f

# OTHER PLOTS
# Barplot of variance explained by global factors per group series in a group
c.group = 2
strtitle = paste("Group ",c.group,": ",groupnames[c.group],"\n"," Variance explained by global factors",sep="") 
data1.group <- data.frame(
  name=seriesnames,
  value=diag(var(category.list.global.common[[c.group]])))

# Create a custom theme
custom_theme <- theme(
  panel.background = element_rect(fill = "lightgray"),
  panel.grid.major = element_line(color = "white"),
  panel.grid.minor = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  axis.title = element_text(size = 12, face = "bold"),
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
)

barplot.p <- ggplot(data=data1.group, aes(x=name, y=value)) +
  geom_bar(stat="identity", width=0.75, fill="steelblue") +
  custom_theme +
  labs(title=strtitle, 
       x=" ", y = " ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  theme(plot.title = element_text(size = 10))

barplot.p 
out3 <- CP(category_list, rmax = 15, r0 = 8, r = rep(1,16), localfactor = TRUE, type = "ER") 

# Υπολογισμός των κοινών παραγόντων για κάθε κατηγορία
category.list.global.common <- list()
category.list.local.common <- list()

for (i in 1:Mgroups) {
  category.list.global.common[[i]] <- out3$Ghat %*% t(out3$loading_G[[i]])
  category.list.local.common[[i]] <- out3$Fhat[[i]] %*% t(out3$loading_F[[i]])
}

# Συνάρτηση για υπολογισμό διακύμανσης
multiple.func <- function(x) { mean(diag(var(x))) }

# Υπολογισμός των τιμών
Table5row1 <- sapply(category.list.global.common, multiple.func)  # Global
Table5row2 <- sapply(category.list.local.common, multiple.func)   # Local
Table5row3 <- Table5row1 + Table5row2  # Total

# Δημιουργία πίνακα
Table5.f <- rbind(Global = Table5row1, Local = Table5row2, Total = Table5row3)

# Ονομασίες στηλών
colnames(Table5.f) <- 1:Mgroups

# Στρογγυλοποίηση στις 2 δεκαδικές
Table5.f <- round(Table5.f, 2)
print(Table5.f)


# **Heatmaps για τα residuals**





# Τοποθέτηση του άξονα x
# Σχεδίαση του heatmap





# Ορισμός των τιμών για τους άξονες του δεύτερου heatmap

# Ορισμός των τιμών για τους άξονες


# Τοποθέτηση του άξονα x



# --- ΕΝΣΩΜΑΤΩΣΗ HEATMAPS (όπως στον δεύτερο κώδικα) ---

# Υποθέτουμε ότι το X_standardized έχει οριστεί ήδη (X_standardized <- scale(X))
# και ότι έχετε υπολογίσει τα global common factors (globalCommon)
# σύμφωνα με την ανάλυση (π.χ., με βάση το out3 ή με άλλο τρόπο).
# Αν δεν υπάρχει ήδη, βεβαιωθείτε ότι ορίζετε globalCommon όπως στο δεύτερο script:
#   globalCommon <- Fhat %*% t(Lambdahat)
# (προσαρμόστε τον υπολογισμό στα δεδομένα σας)
globalCommon_list <- list()
for(i in 1:Mgroups){
  globalCommon_list[[i]] <- out2$Ghat %*% t(out2$loading_G[[i]])
}
globalCommon <- do.call(cbind, globalCommon_list)

Y_global <- X_standardized - globalCommon
# Υπολογισμός των correlation matrices με απόλυτες τιμές
plotHeatmapFull <- function(corrMat, main_title = "", colorbar = TRUE) {
  corrAbs <- abs(corrMat)
  dimMat <- ncol(corrAbs)
  
  color_palette <- colorRampPalette(c("white", "blue", "cyan", "yellow", "red"))(64)
  
  if (colorbar) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4, 1.5))
  } else {
    layout(matrix(1))
  }
  
  par(mar = c(5, 5, 4, 2))
  
  # Heatmap πλήρους μεγέθους (512x512)
  image(
    x = 1:dimMat,
    y = 1:dimMat,
    z = t(corrAbs[dimMat:1, ]),  # αντιστροφή y
    col = color_palette,
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_title
  )
  
  # Θέσεις και ετικέτες για άξονες
  ticks <- seq(50, 500, by = 100)
  
  # Άξονας x (αριστερά προς δεξιά)
  axis(1, at = ticks, labels = ticks)
  
  # Άξονας y (από πάνω προς τα κάτω → rev)
  axis(2, at = ticks, labels = rev(ticks), las = 2)
  
  # Χρωματική μπάρα
  if (colorbar) {
    par(mar = c(5, 2, 4, 4))
    color_levels <- seq(0, 1, length.out = 64)
    image(
      1,
      color_levels,
      matrix(color_levels, nrow = 1),
      col = color_palette,
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    axis(4, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), las = 2)
  }
}

# Υπολογισμός συσχετίσεων
# Υπολογισμός συσχετίσεων
corr_data <- cor(X_standardized)
corr_global <- cor(Y_global)
residual_local <- do.call(cbind, out3$residual)
corr_globloc <- cor(residual_local)

# Χρωματική παλέτα
color_palette <- colorRampPalette(c("white", "blue", "cyan", "yellow", "red"))(64)

# Δημιουργία layout 3 γραμμών × 2 στηλών (heatmap + colorbar)
layout(matrix(c(1,2, 3,4, 5,6), nrow=3, byrow=TRUE), widths = c(4,1.2), heights = c(1,1,1))

# Συνάρτηση για heatmap με χρωματική μπάρα
plotHeatmapWithColorbar <- function(corrMat, main_title = "") {
  corrAbs <- abs(corrMat)
  dimMat <- ncol(corrAbs)
  
  # --- Plot Heatmap ---
  par(mar = c(4, 4, 4, 2))  # bottom, left, top, right
  image(
    x = 1:dimMat,
    y = 1:dimMat,
    z = t(corrAbs[dimMat:1, ]),
    col = color_palette,
    zlim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main_title
  )
  ticks <- seq(50, dimMat, by = 100)
  axis(1, at = ticks, labels = ticks)
  axis(2, at = ticks, labels = rev(ticks), las = 2)
  
  # --- Plot Colorbar ---
  par(mar = c(4, 2, 4, 4))
  color_levels <- seq(0, 1, length.out = 64)
  image(
    1,
    color_levels,
    matrix(color_levels, nrow = 1),
    col = color_palette,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(4, at = seq(0, 1, by = 0.5), labels = seq(0, 1, by = 0.5), las = 2)
}

# Εκτύπωση των 3 heatmaps με colorbars στο ίδιο παράθυρο
plotHeatmapWithColorbar(corr_data, "(α) Συσχέτιση Αρχικών Δεδομένων")
plotHeatmapWithColorbar(corr_global, "(β) Μετά την αφαίρεση καθολικών παραγόντων")
plotHeatmapWithColorbar(corr_globloc, "(γ) Μετά την αφαίρεση καθολικών και τοπικών παραγόντων")

