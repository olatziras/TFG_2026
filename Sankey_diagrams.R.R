library(tidyverse)
library(readxl)
library(networkD3)
library(htmlwidgets)

# ══════════════════════════════════════════════════════════════
#  CONFIGURACIÓN 
# ══════════════════════════════════════════════════════════════

COLORES <- c(
  "K. pneumoniae"  = "#440154",
  "S. typhimurium" = "#90d743",
  "E. coli"        = "#31688e",
  "K. variicola"   = "#fde725",
  "K. oxytoca"     = "#35b779"
)

FUENTE        <- "Arial"   # Arial, Georgia, Helvetica, Verdana, etc.
TAMANO_TEXTO  <- 14        # tamaño en píxeles
ANCHO_NODO    <- 25        # grosor de las barras (px)
SEPARACION    <- 15        # espacio vertical entre nodos (px)
OPACIDAD_LINK <- 0.45      # transparencia de los flujos (0 = invisible, 1 = sólido)
ANCHO_GRAFICO <- 950       # ancho del gráfico (px)
ALTO_GRAFICO  <- 650       # alto del gráfico (px)

# ══════════════════════════════════════════════════════════════
#  CARGA DE DATOS
# ══════════════════════════════════════════════════════════════

ruta_archivo <- file.choose()

if (str_detect(ruta_archivo, "\\.xlsx$")) {
  hojas  <- excel_sheets(ruta_archivo)
  hoja   <- if ("Hoja2" %in% hojas) "Hoja2" else hojas[1]
  df_raw <- read_excel(ruta_archivo, sheet = hoja)
} else {
  df_raw <- read_csv(ruta_archivo)
}

# ══════════════════════════════════════════════════════════════
#  CLASIFICACIÓN (mismas reglas que el chord diagram)
# ══════════════════════════════════════════════════════════════

clasificar_species <- function(raw) {
  conteo <- str_detect(raw, "S_typhimurium") +
            str_detect(raw, "E_coli")        +
            str_detect(raw, "K_pneumoniae")  +
            str_detect(raw, "K_variicola")   +
            str_detect(raw, "K_oxytoca")

  case_when(
    conteo > 1                    ~ NA_character_,  # muro: ambigüedad
    str_detect(raw, "/")          ~ NA_character_,  # muro: barra
    str_detect(raw, "pGUR|pAF")   ~ "TUNEL",        # túnel: se salta
    str_detect(raw, "S_typhimurium") ~ "S. typhimurium",
    str_detect(raw, "E_coli")        ~ "E. coli",
    str_detect(raw, "K_pneumoniae")  ~ "K. pneumoniae",
    str_detect(raw, "AQC")  ~ "K. pneumoniae",
    str_detect(raw, "K_variicola")   ~ "K. variicola",
    str_detect(raw, "K_oxytoca")     ~ "K. oxytoca",
    TRUE ~ NA_character_            # muro: desconocido
  )
}

# ══════════════════════════════════════════════════════════════
#  EXTRACCIÓN DE CADENAS DESDE EL @
# ══════════════════════════════════════════════════════════════

spacer_cols <- names(df_raw)[str_detect(names(df_raw), "Spacer")]

df_largo <- df_raw %>%
  pivot_longer(cols = all_of(spacer_cols),
               names_to  = "Spacer_Pos",
               values_to = "Species_Raw") %>%
  mutate(
    Spacer_Num = as.numeric(str_extract(Spacer_Pos, "\\d+")),
    Raw_Safe   = replace_na(as.character(Species_Raw), ""),
    Raw_Sin_At = str_remove(Raw_Safe, "^@"),
    Species    = clasificar_species(Raw_Sin_At)
  )

arrays_con_at <- df_largo %>%
  filter(str_starts(Raw_Safe, "@")) %>%
  pull(Array_ID) %>%
  unique()

extraer_cadena <- function(array_data) {
  pos_at <- which(str_starts(replace_na(array_data$Raw_Safe, ""), "@"))
  if (length(pos_at) == 0) return(NULL)
  pos_at        <- pos_at[1]
  spacer_at_num <- array_data$Spacer_Num[pos_at]

  anteriores <- array_data %>%
    filter(Spacer_Num < spacer_at_num,
           !Raw_Safe %in% c("-", "", NA)) %>%
    arrange(desc(Spacer_Num)) %>%
    filter(is.na(Species) | Species != "TUNEL")

  cadena <- "K. pneumoniae"
  for (i in seq_len(nrow(anteriores))) {
    sp <- anteriores$Species[i]
    if (is.na(sp)) break
    cadena <- c(cadena, sp)
  }

  if (length(cadena) < 2) return(NULL)

  tibble(
    Origen        = cadena[-length(cadena)],
    Destino       = cadena[-1],
    Nivel_Origen  = seq_along(cadena[-1]) - 1L,
    Nivel_Destino = seq_along(cadena[-1]),
    Array_ID      = array_data$Array_ID[1]
  )
}

df_cadenas <- df_largo %>%
  filter(Array_ID %in% arrays_con_at) %>%
  group_by(Array_ID) %>%
  group_split() %>%
  map_dfr(extraer_cadena)

if (nrow(df_cadenas) == 0) stop("No se encontraron flujos válidos.")

# ══════════════════════════════════════════════════════════════
#  CONSTRUCCIÓN DE NODOS Y LINKS
# ══════════════════════════════════════════════════════════════

conteo_nodos <- bind_rows(
  df_cadenas %>% distinct(Array_ID, Especie = Origen,  Nivel = Nivel_Origen),
  df_cadenas %>% distinct(Array_ID, Especie = Destino, Nivel = Nivel_Destino)
) %>%
  distinct() %>%
  count(Especie, Nivel, name = "n")

nodos_df <- conteo_nodos %>%
  arrange(Nivel, Especie) %>%
  mutate(
    Label = paste0(Especie, " (n=", n, ")"),             # Texto visible (con espacios)
    GrupoColor = str_replace(Especie, " ", "_"),         
    id    = row_number() - 1L
  )

flujos <- df_cadenas %>%
  count(Origen, Nivel_Origen, Destino, Nivel_Destino, name = "Frecuencia")

links_df <- flujos %>%
  left_join(nodos_df, by = c("Origen"  = "Especie", "Nivel_Origen"  = "Nivel")) %>%
  rename(source = id) %>%
  left_join(nodos_df, by = c("Destino" = "Especie", "Nivel_Destino" = "Nivel")) %>%
  rename(target = id) %>%
  select(source, target, value = Frecuencia)

# ══════════════════════════════════════════════════════════════
#  COLORES
# ══════════════════════════════════════════════════════════════

nombres_sin_espacio <- str_replace(names(COLORES), " ", "_")

color_js <- paste0(
  'd3.scaleOrdinal()',
  '.domain([', paste0('"', nombres_sin_espacio, '"', collapse = ", "), '])',
  '.range([',  paste0('"', COLORES, '"', collapse = ", "), '])'
)

# ══════════════════════════════════════════════════════════════
#  GENERAR SANKEY
# ══════════════════════════════════════════════════════════════

sankey <- sankeyNetwork(
  Links        = as.data.frame(links_df),
  Nodes        = as.data.frame(nodos_df),
  Source       = "source",
  Target       = "target",
  Value        = "value",
  NodeID       = "Label",
  NodeGroup    = "GrupoColor",    
  colourScale  = color_js,
  fontSize     = TAMANO_TEXTO,
  fontFamily   = FUENTE,
  nodeWidth    = ANCHO_NODO,
  nodePadding  = SEPARACION,
  LinkGroup    = NULL,
  sinksRight   = FALSE,
  width        = ANCHO_GRAFICO,
  height       = ALTO_GRAFICO
)

# Aplicar opacidad a los links vía JS inline
sankey <- htmlwidgets::onRender(sankey, paste0('
  function(el) {
    d3.select(el).selectAll(".link")
      .style("stroke-opacity", ', OPACIDAD_LINK, ');
  }
'))
# ══════════════════════════════════════════════════════════════
#  GUARDAR
# ══════════════════════════════════════════════════════════════

nombre_archivo <- tools::file_path_sans_ext(basename(ruta_archivo))
nombre_salida  <- paste0("Sankey_", nombre_archivo, ".html")

saveWidget(sankey, file = nombre_salida, selfcontained = TRUE)
message("Guardado: ", nombre_salida)

sankey
