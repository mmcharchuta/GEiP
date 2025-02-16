# Genetyka Ewolucyjna i Populacyjna, ćwiczenia 2
# stwórz wykresy na podstawie plików .bcf:
  # 1) quality
  # 2) mapping quality
  # 3) depth of coverage
# 
stats <- read.table("C:/Users/mikcha1/Downloads/Stats_QualMQDP_ruf_09.txt",
                    col.names = c("QUAL", "MQ", "DP")) # nadaj nazwy kolumnom

# dla C_pyg_26:
# stats <- read.table("C:/Users/mikcha1/Downloads/GEiP/Stats_QualMQDP.txt",
# col.names = c("QUAL", "MQ", "DP"))

# wersja podstawowa: base R graphics:
hist(stats$QUAL, breaks = 30)
hist(stats$MQ, breaks = 30)
hist(stats$DP, breaks = 30)

# przedstaw 3 zmienne na jednym wykresie:
  # stwórz nowy zestaw danych w formacie długim (long format): pierwsza kolumna zawiera wartość, a druga nazwę zmiennej

library(tidyr) # wczytaj paczkę do manipulacji danymi. Jesli jest taka potrzeba, zainstaluj ją komendą install.packages("tidyr")

pivot_longer(stats, # pod nazwą stats_longer,
             cols = c("QUAL", "MQ", "DP"), # wybierz zmienne do długiego formatu (u nas- wszystkie) 
             names_to = "zmienna",  # nazwa kolumny z nazwą zmiennej
             values_to = "wartość") # nazwa kolumny z wartością zmiennej

# zapisz output jako stats_longer
stats_longer <- pivot_longer(stats, # pod nazwą stats_longer,
                             cols = c("QUAL", "MQ", "DP"), # wybierz zmienne do długiego formatu (u nas- wszystkie) 
                             names_to = "zmienna",  # nazwa kolumny z nazwą zmiennej
                             values_to = "wartość") # nazwa kolumny z wartością zmiennej

# użyjemy biblioteki ggplot2 do generowania wykresu. Jesli jest taka potrzeba, zainstaluj ją komendą install.packages("ggplot2")
library(ggplot2) 

ggplot(stats_longer, # nazwa danych
        aes(wartość)) + # na wykresie naniesiemy tylko zmienną wartość
  geom_histogram() + # rysowanie histogramu
  facet_wrap(~zmienna, # na osobnych panelach przedstaw osobno wartości z kolumny "zmienna"
             nrow = 3) # użyj 3 osobnych wierszy
             
# zauważ, że dla wszystkich paneli wartość na osi X jest taka sama, co utrudnia odczytanie informacji. Zmienimy to dodając parametr scales = "free_x"
ggplot(stats_longer,
       aes(wartość)) + 
  geom_histogram() + 
  facet_wrap(~zmienna, 
             nrow = 3, 
             scales = "free_x") # oś X może przyjmować różne wartości dla każdego panelu

# upiększ wykres
ggplot(stats_longer,
       aes(wartość)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~zmienna, 
             nrow = 3, 
             scales = "free") +
  theme_bw() + # zmiana parametrów graficznych
  ggtitle("BCF file, C_ruf_09") # dodanie tytułu

# zapisz wykres jako my_plot
my_plot <- ggplot(stats_longer,
                 aes(wartość)) + 
  geom_histogram() + 
  facet_wrap(~zmienna, 
             nrow = 3, 
             scales = "free_x") +
  theme_bw() + # zmiana parametrów graficznych
  ggtitle("BCF file, C_pyg_26")

stats_longer %>%
  dplyr::filter(zmienna == "DP", wartość < 30) %>%
ggplot(
       aes(wartość)) + 
  geom_histogram(binwidth = 1) + 
  facet_wrap(~zmienna, 
             nrow = 3, 
             scales = "free_x") +
  geom_vline(xintercept = 20) 
  theme_bw() + # zmiana parametrów graficznych
  ggtitle("BCF file, individual A")

ggsave("Stats_QualMQDP_C_ruf_09.png", 
       width = 5, # szerokość (domyślnie w calach, możesz zmienić na cm dodając opcję units = "cm")
       height = 4) # wysokość

