library(tidyverse)
library(circlize)
library(viridis)
library(readxl) # Librería específica para Excel

# 1. Seleccionar el archivo (Soporta .xlsx o .csv)
ruta_archivo <- file.choose()

if (str_detect(ruta_archivo, "\\.xlsx$")) {
  # Si eliges el Excel, leemos la Hoja2
  df_raw <- read_excel(ruta_archivo, sheet = "Hoja2")
} else {
  # Si eliges el CSV
  df_raw <- read_csv(ruta_archivo)
}

# 2. Procesamiento estricto (Gestión de muros y túneles)
df_procesado <- df_raw %>%
  # Pasamos a formato largo
  pivot_longer(cols = contains("Spacer"), 
               names_to = "Spacer_Pos", 
               values_to = "Species_Raw") %>%
  # Limpieza de nombres y clasificación
  mutate(
    # Extraemos el número para el orden cronológico
    Spacer_Num = as.numeric(str_extract(Spacer_Pos, "\\d+")),
    
    # Creamos un string seguro por si hay NAs originales
    Raw_Safe = replace_na(Species_Raw, ""),
    
    # Contamos a cuántas de las 5 especies clave mapea este spacer
    Conteo_Especies = str_detect(Raw_Safe, "S_typhimurium") +
      str_detect(Raw_Safe, "E_coli") +
      str_detect(Raw_Safe, "K_pneumoniae") +
      str_detect(Raw_Safe, "K_variicola") +
      str_detect(Raw_Safe, "K_oxytoca"),
    
    # Normalización a las 5 especies deseadas
    Species = case_when(
      # PRIMERO: Los checkeos de ambigüedad (Muros)
      Conteo_Especies > 1 ~ NA_character_,                
      str_detect(Raw_Safe, "/") ~ NA_character_,          
      
      # SEGUNDO: Checkeo del plásmido propio (Túnel)
      str_detect(Raw_Safe, "pGUR") ~ "IGNORAR_Y_CONECTAR", 
      str_detect(Raw_Safe, "pAF") ~ "IGNORAR_Y_CONECTAR", 
      
      # TERCERO: Las bacterias
      str_detect(Raw_Safe, "S_typhimurium") ~ "S. typhimurium",
      str_detect(Raw_Safe, "E_coli") ~ "E. coli",
      str_detect(Raw_Safe, "K_pneumoniae") ~ "K. pneumoniae",
      str_detect(Raw_Safe, "K_variicola") ~ "K. variicola",
      str_detect(Raw_Safe, "K_oxytoca") ~ "K. oxytoca",
      
      # CUARTO: Todo lo demás
      TRUE ~ NA_character_ # Otros plásmidos y desconocidos se vuelven NA (Muro)
    )
  ) %>%
  # ¡FILTRADO CLAVE ANTES DEL LEAD!:
  # Eliminamos pGUR de la tabla. Al no estar en la secuencia, la bacteria
  # anterior se emparejará directamente con la que venga después.
  # Mantenemos los NA intactos para que sigan bloqueando flujos falsos.
  filter(is.na(Species) | Species != "IGNORAR_Y_CONECTAR")

# 3. Calcular flujos (Trayectorias)
flujos_resumen <- df_procesado %>%
  group_by(Array_ID) %>%
  arrange(Array_ID, desc(Spacer_Num)) %>% #aqui vamos del ultimo spacer al primero (más cercano al leader)
  # Calculamos el destino (pGUR ya no estorba, pero los NA sí frenan el flujo)
  mutate(Destino = lead(Species)) %>%
  ungroup() %>%
  # FILTRO FINAL: Solo nos quedamos con saltos limpios entre las 5 especies.
  filter(!is.na(Species), !is.na(Destino), Species != Destino) %>%
  count(Origen = Species, Destino, name = "Frecuencia")

# 4. Colores y Orden
especies_finales <- c("S. typhimurium", "E. coli", "K. pneumoniae", "K. variicola", "K. oxytoca")
colores_especies <- viridis(5, begin = 0, end = 1)
names(colores_especies) <- especies_finales

# 5. Extraer el nombre del archivo sin la ruta ni la extensión
nombre_archivo <- tools::file_path_sans_ext(basename(ruta_archivo))

# Crear el nombre del archivo de salida
nombre_salida <- paste0("ChordDiagram_5_Especies_", nombre_archivo, "_FINAL.png")

# 6. Generar el Gráfico usando el nombre dinámico
png(nombre_salida, width = 12, height = 12, units = "in", res = 300)

circos.clear()
# Ajustamos los márgenes (canvas.xlim/ylim) para que los nombres no se corten
circos.par(gap.after = 8, canvas.xlim = c(-1.3, 1.3), canvas.ylim = c(-1.3, 1.3))

chordDiagram(
  x = flujos_resumen,
  grid.col = colores_especies,
  order = especies_finales,
  directional = 1,
  direction.type = "diffHeight", 
  diffHeight = -0.04,            
  transparency = 0.35,
  annotationTrack = "grid", 
  preAllocateTracks = list(track.height = 0.15) 
)

# PISTA 1: Nombres de las especies en cursiva
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(x = CELL_META$xcenter, 
              y = CELL_META$ylim + mm_y(6), 
              labels = CELL_META$sector.index, 
              facing = "bending.inside", 
              niceFacing = TRUE, 
              cex = 1.4, 
              font = 3) # Cursiva
}, bg.border = NA)

# PISTA 2: Ejes numéricos
circos.track(track.index = 2, panel.fun = function(x, y) {
  circos.axis(h = "top", labels.cex = 0.7, major.tick.length = mm_y(1.5))
}, bg.border = NA)

dev.off()
circos.clear()