# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal());library("qqplotr")

# Helper Fn ---------------------------------------------------------------
sim_as_tb <- function(sim, sim_type = "hldr"){
  data_name <- names(sim) |> 
    stringr::str_split_i(pattern = "_", 2)
  out <- list()
  if(sim_type == "hldr"){
    for(i in 1:length(data_name)){
      out[[i]] <- lapply(sim[[i]], function(x){
        meds <- apply(x, 2, median)
        min_dim <- which.min(meds)
        x[, min_dim]
      }) |> 
        do.call(cbind, args = _) |> 
        as_tibble() |> 
        mutate("dataset" = rep(data_name[[i]], times = n()))
    }
  } else if(sim_type == "qda"){
    for(i in 1:length(data_name)){
      out[[i]] <- sim[[i]] |> 
        as_tibble() |> 
        mutate("dataset" = rep(data_name[[i]], times = n()))
    }
  }
  out |> 
    do.call(rbind, args =_)
}
# HLDR Boxplots -----------------------------------------------------------
sims_hldr <- miscset::lload("./saved_sims", pattern = "_hldr.RData") |> 
  sim_as_tb() |> 
  filter(dataset != "bands") |> 
  mutate(
  "dataset" = fct_relevel(as_factor(dataset), 
                          "autism", "bc", "div","ion",  "spect", "peng",
                          "wheat", "ecoli", "beans") |> 
    fct_recode("Autism" = "autism",  "Breast Cancer" = "bc", 
               "Divorce" = "div", "Ionosphere" = "ion",  
               "SPECT Heart" = "spect", "Penguins" = "peng",
               "Wheat Seeds" = "wheat", "Ecoli" = "ecoli", 
               "Dry Beans" = "beans") 
  )

sims_long <-  sims_hldr |> 
  select(-QDA) |> 
  rename("MLE" = none, 
         "MRY" = "Molstad") |> 
  pivot_longer(-dataset, names_to = "method", values_to = "cer") |> 
  mutate("method" = fct_relevel(as_factor(method), 
                                "MLE", "Haff", "Wang", "Bodnar", "MRY"))

meds_hldr <- sims_hldr |> 
  summarise("med" = median(QDA), .by = dataset)
  

outliers <- sims_long |>
  group_by(method, dataset) |>
  filter(cer < quantile(cer, .25) - 1.5*IQR(cer) | 
           cer > quantile(cer, .75) + 1.5*IQR(cer)) |>
  ungroup() |>
  mutate("x_dum" = as.numeric(method)
  )

sims_long |> 
  ggplot(aes(method, cer, fill = method))+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, linewidth = .2)+
  geom_point(aes(x = x_dum, cer, shape = method), data = outliers, 
             alpha = .5)+
  geom_hline(aes(yintercept = med, color = "1"), data = meds_hldr)+
  scale_color_manual(values = "red", labels = "QDF")+
  facet_wrap(~dataset, scales = "free")+
  scale_fill_manual(values = c("purple", "#6495ED", "#21918c", "#5ec962", "#fde725"))+
  labs(y = expression(bar("CER")), 
       color = "", fill = "Precision", shape = "Precision")+
  theme(strip.background =element_rect(fill="lightgrey"), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "bottom")+
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(override.aes = list(alpha = 1)))
# 6x5


# QDA Boxplots ------------------------------------------------------------

## Didnt Use in Paper
sims_qda <- miscset::lload("./saved_sims", pattern = "_qda.RData") |> 
  sim_as_tb(sim_type = "qda")

sims_qda |> 
  summarise(across(where(is.numeric), median), .by = dataset) |> 
  pivot_longer(-dataset) |> 
  ggplot(aes(as_factor(name), value, group = 1))+
  geom_line()+
  facet_wrap(~dataset, scales = "free")

meds_qda <- sims_qda |> 
  summarise("med" = median(QDA), .by = dataset)
sims_qda |> 
  select(-QDA) |> 
  pivot_longer(-dataset) |> 
  ggplot(aes(name, value, fill = name))+
  geom_boxplot()+
  geom_hline(aes(yintercept = med, color = "1"), data = meds_qda)+
  scale_color_manual(values = "red", labels = "QDA")+
  facet_wrap(~dataset, scales = "free")+
  scale_fill_manual(values = c("purple", "#3b528b", "#21918c", "#5ec962", "#fde725"))+
  labs(color = "", fill = "SY Shrinkage")+
  theme(strip.background =element_rect(fill="lightgrey"))
