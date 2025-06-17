# Fit A-T parameters
# by Camille Sicangco
# Created 31 Jan 2025

ACi_params = read.csv("data/in/raw/PPC-TGlob_V1.0.csv") %>% 
  filter(Species == "Eucalyptus parramattensis" & 
           Date %in% c("4/10/2016", "5/10/2016", "6/10/2016", "7/10/2016"),
         !(Chamber == 12 & Tleaf > 41)#,
         #Chamber != 5
         ) %>% 
  mutate(TsK = Tleaf + 273.15)
ACi_params.l = split(ACi_params, ACi_params$Chamber)

ACi_params %>% ggplot(aes(x = Tleaf, y = Jmax)) +
  geom_point() +
  theme_classic() +
  facet_wrap(vars(Chamber))

## Fit Vcmax temperature response
fitpeaked(filter(ACi_params.l[[5]], Tleaf < 41), return = "Peak") # Fails for 5
VvT_df = data.frame(do.call(rbind, Map(fitpeaked, ACi_params.l[-5], return = "Peak")))
VvT_df = VvT_df %>% dplyr::mutate(Chamber = row.names(VvT_df),
                                  .before = 1)
rownames(VvT_df) <- NULL

## Fit Jmax temperature response 
JvT_df = data.frame(do.call(rbind, Map(fitpeakedJ, ACi_params.l, return = "Peak")))
JvT_df = JvT_df %>% dplyr::mutate(Chamber = row.names(JvT_df),
                                  .before = 1)
rownames(JvT_df) <- NULL

Tresp_fits = merge(VvT_df, JvT_df, by = c("Chamber"), all = TRUE) %>% 
  mutate(T_treatment = 
           ifelse(Chamber %in% c(1,3,5,7,9,11), "ambient", "elevated"), 
         .after = Chamber)
Tresp_fits %>% 
  #group_by(T_treatment) %>% 
  summarise(across(c(Vcmax25, EaV, delsV, ToptV, Jmax25, EaJ, delsJ, ToptJ), 
                   ~ mean(.x, na.rm = TRUE)))
