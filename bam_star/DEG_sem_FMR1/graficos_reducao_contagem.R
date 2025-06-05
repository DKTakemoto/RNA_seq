library(ggplot2)
library(dplyr)

# Counts da imagem
before <- c(53390, 52302, 47978, 50003, 240898, 213306, 270456, 223489, 4254, 3535, 4466, 4231)
after  <- c(53318, 52243, 47908, 49934, 232478, 205230, 262181, 214862, 4153, 3471, 4383, 4166)

# Nomes reais das amostras
samples <- c(
  "C21", "C22", "C23", "C24",  # C2 grupo
  "C31", "C32", "C33", "C34",  # C3 grupo
  "N1", "N2", "N3", "N4"       # N grupo
)

# Detectar o grupo pelo prefixo
groups <- ifelse(grepl("^C2", samples), "C2",
                 ifelse(grepl("^C3", samples), "C3", "N"))

# Criar o dataframe
df <- data.frame(
  sample = rep(samples, each = 2),
  condition = rep(c("Antes", "Depois"), times = length(samples)),
  value = c(rbind(before, after)),
  group = rep(groups, each = 2)
)

# Gráfico pareado
ggplot(df, aes(x = condition, y = value, group = sample, color = group)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1, alpha = 0.7) +
  geom_text(data = df %>% filter(condition == "Depois"),
            aes(label = group),
            hjust = -0.1,
            size = 3.5,
            show.legend = FALSE) +
  theme_minimal() +
  labs(title = "Expressão do gene FMR1 antes e depois da remoção",
       x = "Condição",
       y = "Contagem") +
  theme(legend.position = "top")
