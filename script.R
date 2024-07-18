################################################################################

read_out_file <- "read_out.xlsx"     # adjust names if filename deviates
design_file <- "design.xlsx"         # adjust names if filename deviates
show_time_plots <- F    # for confirming linear range detection, makes script sig. slower
result_dir <- "./results/"
image_type <- ".jpg"

sample_volume_in_well_in_µl <- 40
assay_volume_in_µl <- 200
dilution_factor <- sample_volume_in_well_in_µl/assay_volume_in_µl
light_path_length_in_cm <- 0.5
extinction_coeff_in_M_per_cm <- 18600

boxplots <- T # if False, dotplots will be created

# ! TAKE CARE

# ! make sure the input files are not open in another program
# ! make sure the time column in the readout file contains the word "min" somewhere
# ! the sample_table     has to be marked with a "#1"
# ! the treatment_table  has to be marked with a "#2"


# you can run the code! :D

################################################################################

# Load libraries
library(readxl)
library(ggplot2)
library(segmented)
library(dplyr)

# set up the environment--------------------------------------------------------

# Set working directory to the current directory the file is in
wor_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wor_dir)

# Create and set result_dir
dir.create(result_dir, showWarnings = FALSE)

# Read input files
read_out <- read_xlsx(read_out_file)
design <- read_xlsx(design_file)

NA_list <- c("Empty", "empty", "emtpy",
             "NA", "na",
             "Na", "EMPTY", "")

# subset and unify the files for further processing-----------------------------

# data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
header <- which(is.na(read_out[,ncol(read_out)-5])) #check for empty cells at the far right
data <- as.data.frame(read_out[-header,])

#getting the times out of the excel and make it rownames
bool <- apply(data, c(1,2), function(x){grepl("min", x)}, simplify = T)
row_name_col <- which(colSums(bool) != 0) #the column that contains the rownames

if (is.na(data[1, row_name_col])){    # set min as the header for the column if
  data[1, row_name_col] <- "min"      # it is not already set
}

rownames(data) <- data[, row_name_col] 
data <- data[,(row_name_col+1):ncol(data)]

#setting the colnames
colnames(data) <- data[1,]
data <- data[2:nrow(data),] #get rid of the column headers of the data file

data <- as.data.frame(lapply(data, function(x) as.numeric(as.character(x))))
#this somehow (idk how) also gets rid of the "min" in the rownames

# design ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subset_table <- function(design, hash){
  #this function is used to extract the annotation tables from the big excel.
  #It works by searching the table using the hash and then subsetting the table
  #based on assumed size of 8x12
  
  tab_start <- which(design==hash, arr.ind = T) + c(1, 0)
  tab_end <- tab_start + c(7, 12)
  
  tab <- as.data.frame(design[tab_start[1]:tab_end[1],
                       tab_start[2]:tab_end[2]])
  
  if ("A" %in% tab[,1]){
    rownames(tab) <- tab[,1]
    tab <- tab[,-1]
  } else {
    rownames(tab) <- LETTERS[1:nrow(tab)]
  }
  
  colnames(tab) <- c(1:ncol(tab))
  
  return(tab)
}

sample_table <- subset_table(design, "#1")
treatment_table <- subset_table(design, "#2")

wells <- colnames(data)

# Create a function to map well names to respective meta data
get_meta <- function(well_name, table) {
  row <- substr(well_name, 1, 1)             # subsets the first letter
  col <- as.numeric(substr(well_name, 2, 3)) # subsets the number
  return(table[row, col])                    # returns the element in the cell
}

samples <- sapply(wells, function(well) get_meta(well, sample_table))
treatments <- sapply(wells, function(well) get_meta(well, treatment_table))

well_df <- data.frame(well = wells,
                      sample = samples,
                      treatment = treatments,
                      slope = rep(NA, length(wells)),
                      SEAP = rep(NA, length(wells)))

#exclude wells without sample
well_df <- well_df[!well_df$sample %in% NA_list,]
well_df$sample <- as.factor(well_df$sample)
well_df$treatment <- as.factor(well_df$treatment)

#get the linear range using the segmented_package-------------------------------

for (sample in colnames(data)){
           
  df <- data.frame(x = 1:nrow(data),
                   y = data[[sample]])
  
  if(all(abs(df$y) > 0.1)){ #then this well is not empty
  
    # Fit a linear model
    lm_fit <- lm(y ~ x, data = df)
    
    # Apply the segmented package to detect breakpoints
    seg_fit <- segmented(lm_fit, seg.Z = ~ x, npsi = 1)
    
    slope <- max(coef(seg_fit)[2:3])
    
    well_df[sample, "slope"] <- slope
    
    if (show_time_plots){
      breakpoints <- seg_fit$psi[, "Est."]
      
      p<- ggplot(df, aes(x = x, y = y)) +
        geom_point() +
        geom_line(aes(y = fitted(lm_fit)), color = "blue", linetype = "dashed") +
        geom_line(aes(y = fitted(seg_fit)), color = "red") +
        geom_vline(xintercept = breakpoints, linetype = "dotted", color = "black") +
        theme_minimal() +
        labs(title = paste(sample, "Segmented Linear Regression"), x = "X", y = "Y")
      print(p)
    }
  }
}

absorption <- extinction_coeff_in_M_per_cm*light_path_length_in_cm

well_df$SEAP <- (well_df$slope/absorption)*(1/dilution_factor)*1e6

# Calculate mean and standard error for each sample/treatment combination
final <- well_df %>%
  group_by(sample, treatment) %>%
  summarize(mean_SEAP = mean(SEAP),
            SE_SEAP = sd(SEAP)/sqrt(n())) #standard error of the mean

# Significance tests -----------------------------------------------------------

unq_samples <- unique(samples)
unq_samples <- unq_samples[!unq_samples %in% NA_list]

#getting the number of treatments
n_treatments <- length(unique((treatments[!treatments %in% NA_list])))

#if there are two treatments:
if(n_treatments <= 1){
  
  #no statistical tests can be done
  
} else if (n_treatments == 2) {

  test_used <- "Wilcoxon rank sum exact test"
  
  perform_statistical_test <- function(sample_name) {
    # Filter the data for the specified sample
    sample_data <- well_df %>% filter(sample == sample_name)
    
    # Perform the Wilcoxon rank-sum test
    test_result <- wilcox.test(SEAP ~ treatment, data = sample_data)
    
    return(test_result$p.value)
  }
} else if (n_treatments > 2) {
  
  test_used <- "Kruskal-Wallis rank sum test"
  
  perform_statistical_test <- function(sample_name) {
    # Filter the data for the specified sample
    sample_data <- well_df %>% filter(sample == sample_name)
    
    # Perform the Kruskal-Wallis test
    test_result <- kruskal.test(SEAP ~ treatment, data = sample_data)
    
    return(test_result$p.value)
  }
}
  
test_results <- lapply(unq_samples, perform_statistical_test)
names(test_results) <- unq_samples

# Example: Perform the test for "Sample 1"
#test_result <- perform_statistical_test("Sample 1")
#print(test_result)

#exporting the resulting data---------------------------------------------------

write.csv(final, paste0(result_dir, "SEAP_readout.csv"))
write.csv(test_results, paste0(result_dir, "p_values.csv"))

#plot and safe the plots--------------------------------------------------------

custom_save <- function(plot, path, name, size=1){
  ggsave(paste0(path, name), plot,
         width = 7000*size,
         height = 4000*size,
         units = "px",
         dpi = 500)
}

#plot each ind. sample

for (sample in unq_samples){
  df <- well_df[which(well_df$sample == sample),]
  if (boxplots){
    p <- ggplot(df, aes(x = treatment, y = SEAP)) +
          geom_boxplot() +
          labs(x = "condition", y = "SEAP_value", title = sample, 
               subtitle = paste0(test_used, ": p = ", round(test_results[[sample]], digits = 3))) +
          theme_classic()
  } else {
    p <- ggplot(df, aes(x = treatment, y = SEAP)) +
      geom_boxplot() +
      labs(x = "condition", y = "SEAP_value", title = sample, 
           subtitle = paste0(test_used, ": p = ", round(test_results[[sample]], digits = 3))) +
      theme_classic()
  }
  print(p)
  custom_save(p, result_dir, paste0(sample, image_type), size = .5)
}

# Plotting
p <- ggplot(final, aes(x=sample, y=mean_SEAP, fill=treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=mean_SEAP-SE_SEAP, ymax=mean_SEAP+SE_SEAP), width=.2, position=position_dodge(.9)) +
  labs(title="SEAP_values by Treatment and Sample", x="Treatment", y="Mean SEAP",
       subtitle = "error bars represent standard error of the mean (SEM)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
custom_save(p, result_dir, paste0("Barplot_by_sample", image_type), size = 1)

p <- ggplot(final, aes(x=treatment, y=mean_SEAP, fill=sample)) +
      geom_bar(stat="identity", position=position_dodge()) +
      geom_errorbar(aes(ymin=mean_SEAP-SE_SEAP, ymax=mean_SEAP+SE_SEAP), width=.2, position=position_dodge(.9)) +
      labs(title="SEAP_values by Treatment and Sample", x="Treatment", y="Mean SEAP",
           subtitle = "error bars represent standard error of the mean (SEM)") +
      theme_minimal()+ 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
custom_save(p, result_dir, paste0("Barplot_by_treatment", image_type), size = 1)
