# importing required packages 
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)     
library(stringr)

multiplesheets <- function(fname) { 
  
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  
  # assigning names to data frames 
  names(data_frame) <- sheets 
  
  # print data frame 
  print(data_frame) 
} 

# specifying the path name 
path <- "/Users/bwremjoe/Dropbox/03May24_version/Simulation_Mastersheet.xlsx"
all_data <- multiplesheets(path)

data_to_plot <- all_data[["Sim330to321"]][,1:11] %>% pivot_longer(cols = !TimeSteps) %>% 
  mutate(producer=str_count(name,"1"))

num_cg <- 6
cc <- scales::seq_gradient_pal("#000000", "#00CC00", "Lab")(seq(0,1,length.out=4))

all_data[["Sim330to321"]][,1:11] %>% pivot_longer(cols = !TimeSteps) %>% 
  mutate(producer=str_count(name,"1"))
  
plot <- all_data[["Sim330to321"]][,1:11] %>% pivot_longer(cols = !TimeSteps) %>% 
  mutate(producer=str_count(name,"1")) %>% 
  ggplot(aes(x=TimeSteps,y=value,grp=name,col=paste0(producer,'-producer'))) +
  geom_line() +
  theme_classic() +
  xlab("Time") +
  ylab("Abundance") +
  scale_y_continuous(breaks=c(0,20000)) +
  scale_colour_manual(values=cc,name="") +
  theme(axis.title.y = element_text(vjust = -6.0))

ggsave("330to321.png", plot=plot,dpi=70)
        