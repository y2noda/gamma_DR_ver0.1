library(tidyverse)
library(ggplot2)
library(gt)


## アウトプット
# export_data <- data.frame(results.beta_t)
# # write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/beta_hat_m00.csv")
# # 
# export_data <- data.frame(results.gamma)
# write_csv2(export_data, file = "~/Projects/gamma_DR_ver0.1/results/gamma_hat_m00.csv")

## input
df <- read_csv2("~/Projects/gamma_DR_ver0.1/results/beta_hat_m03.csv") 
df %>% summary()


## data加工
make_df_fn <- function(df, est_type, mp_num, file_path){
  
  beta <- read_csv2(file_path) 
  
  n <- nrow(beta)
  
  type <- rep(est_type, n)
  mp <- rep(mp_num, n)
  
  new_df <- data.frame(type, mp, beta)
  
  names(new_df) <- c("Type","MP","Beta")
  
  results <- rbind(df, new_df) 
  
  return(results)
}

est_types <- c("naive", "hat")
mp_nums <- c("00","005","01","015","02","03")

df <- data.frame()
for (est_type in est_type) {
  for(mp_num in mp_nums){
    file_path <- paste("~/Projects/gamma_DR_ver0.1/results/", "beta_",est_type ,"_m" ,mp_num ,".csv" ,sep="")
    df <- make_df_fn(df, est_type,mp_num,file_path)
  }
}
# write_csv2(df, file = "~/Projects/gamma_DR_ver0.1/results/all_beta.csv")

df %>% dim


## plot

# path <- "~/gamma_DR_ver0.1/results/all_beta.csv"
path <- "~/Projects/gamma_DR_ver0.1/results/all_beta.csv"
dat <- read_csv2(path)


p2 <- ggplot(dat, aes(x = MP, y = Beta, fill = factor(Type))) +
  geom_boxplot() +
  scale_y_continuous(name = "Beta") +
  scale_x_discrete(labels = abbreviate, name = "Miss probability")

p2

dat <- dat %>% 
  group_by(Type, MP) %>% 
  summarise(mean = mean(Beta),
            sd = sd(Beta)) %>% 
  ungroup()

beta_table <- dat %>%
  gt(rowname_col = "Type", groupname_col = "MP")

beta_table

# gtsave(beta_table, "~/Projects/gamma_DR_ver0.1/fig/beta_table.tex")



### gamma

gamma_list <- NULL

mp_nums <- c("00","005","01","015","02","03")

for (mp_num in mp_nums) {
  file_path <- paste("~/Projects/gamma_DR_ver0.1/results/", "gamma_hat_m" ,mp_num ,".csv" ,sep="")
  gmean <- read_csv2(file_path)$results.gamma %>% as.vector() %>% mean()
  
  gamma_list <- append(gamma_list,gmean) 
}

# 
# path <- "~/Projects/gamma_DR_ver0.1/results/gamma_hat_m03.csv"
# path <- "~/gamma_DR_ver0.1/results/gamma_hat_m03.csv"
# 
# mp_list <- c("0.05","0.1","0.15","0.2","0.3")

gamma_df <- data.frame("miss prob"=mp_nums, "gamma.mean"=gamma_list) 

# write_csv2(gamma_df, file = "~/Projects/gamma_DR_ver0.1/results/all_gamma_mean.csv")

path <- "~/Projects/gamma_DR_ver0.1/results/all_gamma_mean.csv"
# path <- "~/gamma_DR_ver0.1/results/all_gamma_mean.csv"
dat <- read_csv2(path)

gamma_table <- dat %>% gt()

gamma_table

gtsave(gamma_table, "~/Projects/gamma_DR_ver0.1/fig/gamma_table.tex")