library(GrFA)
library(ggplot2)
library(reshape2)

# Καθαρισμός περιβάλλοντος
graphics.off() ; rm(list = ls(all = TRUE)) ; cat("\014")

# Διάβασμα των δεδομένων
data_developed <- read.csv("dev1.csv", header = F)
data_emerging <- read.csv("eme1.csv", header = F)

# Προσθήκη χρονικού δείκτη
data_developed$Time <- seq(as.Date("1996-07-01"), by = "quarter", length.out = nrow(data_developed))
data_emerging$Time <- seq(as.Date("1996-07-01"), by = "quarter", length.out = nrow(data_emerging))

# Προσθήκη στήλης κατηγορίας
data_developed$Category <- "Ανεπτυγμένες"
data_emerging$Category <- "Αναπτυσσόμενες"

# Μετατροπή σε long format
data_developed_long <- reshape2::melt(data_developed, id.vars = c("Time", "Category"))
data_emerging_long <- reshape2::melt(data_emerging, id.vars = c("Time", "Category"))

# Συγχώνευση των δεδομένων
data_combined <- rbind(data_developed_long, data_emerging_long)

# Δημιουργία γραφήματος χωρίς legend μεταβλητών
ggplot(data_combined, aes(x = Time, y = value, group = variable, color = Category)) +
  geom_line() +  # Χρησιμοποιούμε μόνο το χρώμα για να ξεχωρίζουν οι ομάδες
  labs(
    title = "Χρονοσειρές για Ανεπτυγμένες και Αναπτυσσόμενες Αγορές",
    x = "Χρόνος: Τριμηνιαία δεδομένα",
    y = "Τιμή *",
    color = "Κατηγορία"
  ) +
  theme_minimal() +
  theme(legend.position = "right") +  # Διατήρηση legend μόνο για τις κατηγορίες
  scale_x_date(
    date_labels = "Q%3/%y",
    breaks = seq(as.Date("1996-07-01"), by = "36 months", length.out = ceiling(nrow(data_developed) / 12))
  ) +
  guides(color = guide_legend(title = "Κατηγορία"), linetype = "none")  # Αφαίρεση legend μεταβλητών
