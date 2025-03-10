# Packages ----------------------------------------------------------------
library("tidyverse");theme_set(theme_minimal())

# Helper Fn ---------------------------------------------------------------
sim_as_tb <- function(sim, sim_type = "data_sdr"){
  out <- list()
  if(sim_type == "sim"){
    name <- names(sim)
    for (i in 1:length(name)) {
      out[[i]] <- sim[[i]] |> 
        mutate("method" = as_factor(method)) |> 
        group_by(ni, method) |> 
        group_split() |> 
        lapply(function(x){
          meds <- x |> 
            select(-method, -ni) |> 
            apply(2, median)    
          x |> 
            select(c(which.min(meds), method, ni)) |> 
            setNames(c("error", "method", "ni"))
        }) |> 
        do.call(rbind, args = _) |> 
        mutate("config" = rep(name[[i]], n()), 
               "ni" = rep(c("p + 1", "2p", "6p"), each = n()/3)
        )
    }
  } else if(sim_type == "data_sdr"){
    name <- names(sim) |> 
      stringr::str_split_i(pattern = "_", 2)
    for(i in 1:length(name)){
      out[[i]] <- lapply(sim[[i]], function(x){
        meds <- apply(x, 2, median)
        min_dim <- which.min(meds)
        x[, min_dim]
      }) |> 
        do.call(cbind, args = _) |> 
        as_tibble() |> 
        mutate("dataset" = rep(name[[i]], times = n()))
    }
  } else if(sim_type == "qda"){
    name <- names(sim) |> 
      stringr::str_split_i(pattern = "_", 2)
    for(i in 1:length(name)){
      out[[i]] <- sim[[i]] |> 
        as_tibble() |> 
        mutate("dataset" = rep(name[[i]], times = n()))
    }
  } else if(sim_type == "comp"){
    name <- names(sim) |> 
      stringr::str_split_i(pattern = "_", 2)
    for(i in 1:length(name)){
      out[[i]] <- sim[[i]] |> 
        as_tibble() |> 
        mutate("dataset" = rep(name[[i]], times = n()))
    }
  }
  out |> 
    do.call(rbind, args =_)
}


# MC Sim SSDR Boxplots -------------------------------------------------------
sims <- miscset::lload("./saved_sims", pattern = "config") |> 
  sim_as_tb(sim_type = "sim") |> 
  mutate("config" = case_when(config == "config1" ~ "Configuration 1", 
                              config == "config2" ~ "Configuration 2",
                              config == "config3" ~ "Configuration 3",
                              config == "config4" ~ "Configuartion 4") |> 
           as_factor(), 
         "ni" = as_factor(ni))

qdf_meds <- sims |> 
  filter(method == "qdf") |> 
  group_by(ni, config) |> 
  summarise("med" = median(error)) |> 
  ungroup()

outliers <- sims |>
  filter(method != "qdf") |> 
  group_by(method, ni, config) |>
  filter(error < quantile(error, .25) - 1.5*IQR(error) | 
           error > quantile(error, .75) + 1.5*IQR(error)) |>
  ungroup() |>
  mutate("x_dum" = as.numeric(method)
  )

sims |> 
  filter(method != "qdf") |> 
  ggplot(aes(method, error, fill = method))+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, linewidth = .2)+
  geom_point(aes(x = x_dum, error, shape = method), data = outliers,
             alpha = .25)+
  geom_hline(aes(yintercept = med, color = "1"), data = qdf_meds)+
  scale_color_manual(values = "red", labels = "QDA")+
  scale_fill_manual(values = c("purple", "#6495ED", "#21918c", "#5ec962", "#fde725"), 
                    labels = c(bquote("SSDR"["S"^{-1}]), expression("SSDR"["Haff"]), 
                               expression("SSDR"["Wang"]), expression("SSDR"["Bod"]), 
                               expression("SSDR"["MRY"])))+
  scale_shape_manual(values = c(1, 2, 3, 4, 5), 
                     labels = c(bquote("SSDR"["S"^{-1}]), expression("SSDR"["Haff"]), 
                                expression("SSDR"["Wang"]), expression("SSDR"["Bod"]), 
                                expression("SSDR"["MRY"])))+
  facet_grid(cols = vars(ni), rows = vars(config), scales = "fixed",
             labeller = label_bquote(cols = n[i] == .(as.character(ni)))
  )+
  theme(strip.background = element_rect(fill="lightgrey"), 
        legend.position = "bottom", 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank())+
  labs(y = expression("CER"), 
       color = "", fill = "", shape = "")+
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(override.aes = list(alpha = 1)))


# Real Data Application Competitors --------------------------------------------------------------------
sims_comp <- miscset::lload("./saved_sims", pattern = "comp") |> 
  sim_as_tb(sim_type = "comp")

sims_sdr <- miscset::lload("./saved_sims", pattern = "_sdrs.RData") |> 
  sim_as_tb() |> 
  select("Molstad", "QDA")

tb <- sims_comp |> 
  bind_cols(
    sims_sdr
  ) |> 
  relocate(QDA, .before = everything()) |> 
  relocate(Molstad, .before = dataset) |> 
  mutate(
    "dataset" = fct_relevel(as_factor(dataset), 
                            "autism", "bc", "div","ion",  "spect", "peng",
                            "wheat", "ecoli", "beans") |> 
      fct_recode("Autism" = "autism",  "Breast Cancer" = "bc", 
                 "Divorce" = "div", "Ionosphere" = "ion",  
                 "SPECT Heart" = "spect", "Penguins" = "peng",
                 "Wheat Seeds" = "wheat", "Ecoli" = "ecoli", 
                 "Dry Beans" = "beans") 
  ) |> 
  rename("MRY" = "Molstad", 
         "MSDA" = "msda", 
         "QDAP" = "qdap", 
         "SIRL" = "lassoSir", 
         "ENDS" = "ends")

qda_med <- tb |> 
  summarise("med" = median(QDA), .by = dataset)


sim_long <- tb |> 
  select(-QDA) |> 
  pivot_longer(-dataset, names_to = "method", values_to = "cer") |> 
  mutate("method" = fct_relevel(as_factor(method), 
                                "MSDA", "SIRL", "ENDS","QDAP", "MRY")
  ) |> 
  drop_na()

outliers <- sim_long |> 
  group_by(method, dataset) |>
  filter(cer < quantile(cer, .25) - 1.5*IQR(cer) | 
           cer > quantile(cer, .75) + 1.5*IQR(cer)) |>
  ungroup() |>
  mutate("x_dum" = as.numeric(method)
  )

sim_long |> 
  ggplot(aes(method, cer, fill = method))+
  stat_boxplot(geom='errorbar')+
  geom_boxplot(outlier.shape = NA, linewidth = .2)+
  geom_point(aes(x = x_dum, cer, shape = method), data = outliers, 
             alpha = .5)+
  geom_hline(aes(yintercept = med, color = "1"), data = qda_med)+
  facet_wrap(~dataset, scales = "free_y", nrow = 3)+
  scale_shape_manual(values = 1:5, labels = c("MSDA", "SIRL", "ENDS","QDAP", expression("SSDR"["MRY"])))+
  scale_color_manual(values = "red", labels = "QDA")+
  scale_fill_manual(values = c("#F28E2B", "#B07AA1", "#4E79A7", "#76B7B2", "#fde725"), 
                    labels = c("MSDA", "SIRL", "ENDS","QDAP", expression("SSDR"["MRY"])))+
  theme(strip.background =element_rect(fill="lightgrey"), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        legend.position = "bottom")+
  ggh4x::facetted_pos_scales(y = list(
    dataset ==  "SPECT Heart" ~ scale_y_continuous(limits = c(0.15, 0.6)),
    dataset %in% c("Ionosphere", "Autism", "Ecoli") ~
      scale_y_continuous(limits = c(0, 0.25)), 
    dataset %in% c("Dry Beans", "Penguins", "Divorce") ~
      scale_y_continuous(limits = c(0, 0.25)), 
    dataset %in% c("Breast Cancer", "Wheat Seeds") ~
      scale_y_continuous(limits = c(0, 0.25))
  ))+
  labs(color = "", fill = "", y = expression(bar("CER")), shape = "")+
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(override.aes = list(alpha = 1)))

# Real Data Application SSDR Boxplots (Not Used) -----------------------------------------------------------
sims_sdr <- miscset::lload("./saved_sims", pattern = "_sdrs.RData") |> 
  sim_as_tb() |> 
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

sims_long <-  sims_sdr |> 
  select(-QDA) |> 
  rename("MLE" = none, 
         "MRY" = "Molstad") |> 
  pivot_longer(-dataset, names_to = "method", values_to = "cer") |> 
  mutate("method" = fct_relevel(as_factor(method), 
                                "MLE", "Haff", "Wang", "Bodnar", "MRY"))

meds_sdr <- sims_sdr |> 
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
  geom_hline(aes(yintercept = med, color = "1"), data = meds_sdr)+
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


