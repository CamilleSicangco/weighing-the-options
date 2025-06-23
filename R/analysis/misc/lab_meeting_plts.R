Tair_vec = seq(10,70, by = 1)

kmax = calc_kmax(kmax_25 = 0.5, Tair = Tair_vec)
data.frame(Tair = Tair_vec, kmax) %>% 
  ggplot(aes(x = Tair, y = kmax)) +
  geom_line() +
  theme_classic() +
  ylab(expression("k"[max] * " (mmol m"^-2*" s"^-1*" MPa"^-1 *")")) +
  xlab(expression("T"[air]*" (\u00B0C)")) +
  theme(text = element_text(size = 14))

A_out = Photosyn_custom(Tleaf = Tair_vec,
                Tcrit = 43.4, T50 = 44.6,
                #Vcmax=34,EaV=62307,EdVC=2e5,delsC=639,
                #Jmax = 60,EaJ=33115,EdVJ=2e5,delsJ=635, Rd0 = 0.92,
                new_JT = FALSE)
Anet = A_out[,2]
Rd = A_out[,8]
Agross = Anet + Rd
Aj = A_out[,6]
plot(Tair_vec, Aj)

linetype = c(Anet = "solid", Agross = "dashed", Rd = "dotted")
data.frame(Tair = rep(Tair_vec, 3), 
           A = c(Anet, Agross, Rd),
           var = rep(c("Anet", "Agross", "Rd"), each = length(Tair_vec))) %>%
  ggplot(aes(x = Tair, y = A, linetype = var)) +
  geom_line() +
  theme_classic() +
  ylab(expression("A or R"[d]*" ("*mu*"mol m"^-2*" s"^-1*")")) +
  xlab(expression("T"[air]*" (\u00B0C)")) +
  theme(text = element_text(size = 14)) +
  scale_linetype_manual(values = linetype)

Jold = TJmax(Tair_vec,EaJ=33115,EdVJ=2e5,delsJ=635)
Jnew = TJmax_updated(Tair_vec,EaJ=33115,EdVJ=2e5,delsJ=635, Tcrit = 43.4, T50 = 47.6)
data.frame(Tair = rep(Tair_vec, 2), 
           Jmax = c(Jold, Jnew),
           type = rep(c("old", "new"), each = length(Tair_vec))) %>%
  ggplot(aes(x = Tair, y = Jmax, linetype = type)) +
  geom_line() +
  theme_classic() +
  ylab("Jmax/Jmax25") +
  xlab(expression("T"[leaf]*" (\u00B0C)")) +
  theme(text = element_text(size = 14))

T50_49.6 = df %>% 
  filter(ID %in% c("HC_varkmax", "CG_net_newJT", "TC")) %>% 
  filter(P != min(df$P)) %>% 
  pivot_wider(names_from = "ID", values_from = "cost_gain") %>% 
  mutate(Profit = CG_net_newJT - (HC_varkmax + TC)) %>% 
  select(P, Profit)
T50_44.5 = df2 %>% 
  filter(ID %in% c("HC_varkmax", "CG_net_newJT", "TC")) %>% 
  filter(P != min(df2$P)) %>% 
  pivot_wider(names_from = "ID", values_from = "cost_gain") %>% 
  mutate(Profit = CG_net_newJT - (HC_varkmax + TC)) %>% 
  select(P, Profit)

bind_rows(T50_49.6, T50_44.5, .id = "T50") %>% 
  ggplot(aes(x = P, y = Profit, linetype = T50)) + 
  geom_line(linewidth = 1) +
  scale_linetype_manual(values = c("solid", "dashed"),
                        labels = c("49.6",
                                   "44.5"
                        )) +
  theme_classic() +
  xlab("Leaf water potential (-MPa)")
