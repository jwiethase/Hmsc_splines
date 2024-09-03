setwd(file.path(dirname(rstudioapi::getSourceEditorContext()$path),"../"))
rm(list = ls())
library(gridExtra)
library(stringr)
library(Hmsc)
source('source/misc_functions.R')

model_choices <- c("splines", "poly")

for(mod in model_choices){
   # Load data ===================
   load('spatCampaign_2023_HPC_TFfitted_chains4_samples1000_thin10_10.RData')
   m_PA_full <- thin_list[["10"]][['PA']][[mod]]
   m_COP_full <- thin_list[["10"]][['COP']][[mod]]
   
   ## Overall richness & evenness plots ===================
   gradient_plot_list <- list()
   chosen_covars <- c("Salinity", "Temp", "TN_ug_L", "TP_ug_L")
   
   for(i in 1:length(chosen_covars)){
      print(chosen_covars[i])
      if(grepl("spline", mod)){
            Gradient_i = constructGradient_JHW(m_PA_full, 
                                               focalVariables = c(chosen_covars[i], str_subset(names(m_PA_full$XData), paste0(chosen_covars[i], "_spline"))),
                                               ngrid = 50,
                                               type = "marginal")
      } else {
            Gradient_i = constructGradient_JHW(m_PA_full, 
                                               focalVariables = chosen_covars[i],
                                               ngrid = 50,
                                               type = "marginal")
      }
      
      pred_PA_i = predict(m_PA_full, Gradient = Gradient_i, expected = TRUE)
      pred_COP_i = predict(m_COP_full, Gradient = Gradient_i, expected = TRUE)
      pred_combined_i <- Map('*', pred_PA_i, lapply(pred_COP_i, function(x){exp(x)}))
      
      gradient_plot_list[[paste0(chosen_covars[i], "_S")]] <- plotGradient_JHW(m_PA_full, 
                                                                               Gradient_i, 
                                                                               pred = pred_PA_i, 
                                                                               measure = "S")
      gradient_plot_list[[paste0(chosen_covars[i], "_E")]] <- plotGradient_JHW(m_COP_full, 
                                                                               Gradient_i, 
                                                                               pred = pred_combined_i, 
                                                                               measure = "E", 
                                                                               back_transform_Y = 'log')
      gradient_plot_list[[paste0(chosen_covars[i], "_Y")]] <- plotGradient_JHW(m_PA_full, 
                                                                               Gradient_i, 
                                                                               pred = pred_PA_i, 
                                                                               measure = "Y", 
                                                                               index = 1)
   }
   
   png(filename = paste0('figures/gradient_plot_SE_', mod, '.png'), width = 37.5, height = length(gradient_plot_list)*5, units = "cm", res = 300)
   do.call(grid.arrange, c(gradient_plot_list, ncol = 3))
   dev.off()
}



