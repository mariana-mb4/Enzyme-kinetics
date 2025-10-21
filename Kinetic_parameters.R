# Instalar librerías si es que no se tienen ya instaladas
# drc es para modelos no lineales de dosis-respuesta
# renz es solamente para cinetica enzimatica

# install.packages("drc")
# install.packages("renz")
# install.packages("tidyverse")
# install.packages("dplyr")
# install.packages("ggplot2")

library(tidyverse)
# Dataframe
datos_1 <- data.frame(Concentracion = c(0.00003, 0.00006, 0.0001, 0.0002, 0.0004, 0.0006, 0.0008, 0.001, 0.0012),
                    v_inicial1 = c(0.0491, 0.0842, 0.1327, 0.2183, 0.2802, 0.3088, 0.3142, 0.3379, 0.3402),
                    v_inicial2 = c(0.0441, 0.095, 0.1485, 0.2385, 0.2949, 0.3167, 0.3345, 0.3389, 0.3423),
                    v_inicial3 = c(0.0447, 0.0936, 0.1565, 0.2284, 0.282, 0.3008, 0.3252, 0.3319, 0.3155))
print(datos_1)

#Transformar el dataframe para poder graficar múltiples variables en el eje x
datos_2 <- datos_1 |>
  pivot_longer(
    cols = starts_with("v_inicial"), # Selecciona las columnas a "pivotar"
    names_to = "experimento",        # Crea la nueva columna para los nombres
    values_to = "v0"        # Crea la nueva columna para los valores
  )

# Media y desviacion estandar
datos_resumen <- datos_2 |> group_by(Concentracion) |> 
  summarise(media_v0 = mean(v0, na.rm = TRUE),
            sd_v0 = sd(v0, na.rm = TRUE))

########## DRC
library(drc)
# Estimar los valores de Km y Vmax
# Vmax como la v_inicial más alta
vmax_predicha <- max(datos_2$v0, na.rm = TRUE)
# Km = Vmax/2
km_predicha <- vmax_predicha/2
print(c(vmax_predicha, km_predicha))

# Ajustar el modelo de Michaelis-Menten con drm() función para dosis-response
modelo_ajustado <- drm(v0 ~ Concentracion, data = datos_2, fct = MM.2(), 
                       start = c(Vmax = vmax_predicha, Km = km_predicha))

summary(modelo_ajustado)

# Obtener los parámetros cinéticos coef()
parametros <- coef(modelo_ajustado)
Vmax <- parametros["d:(Intercept)"]
Km <- parametros["e:(Intercept)"]

# Calcular el valor de R^2 forma manual
y_experimental <- datos_2$v0
y_predicho <- predict(modelo_ajustado)
ss_total <- sum((y_experimental - mean(y_experimental))^2)
ss_residual <- sum((y_experimental - y_predicho)^2)
r_cuadrado <- 1 - (ss_residual / ss_total)

# Constantes
cat("Velocidad Máxima (Vmax):", round(Vmax, 8), "\n")
cat("Constante de Michaelis (Km):", round(Km, 8), "\n")
cat("Valor de R^2:", round(r_cuadrado, 6), "\n")

########## RENZ
library(renz)
#Gráfico con RENZ
dir.MM(datos_2[ , c(1,3)], unit_S = "mM", unit_v = "nM/min")

# Visualizar gráfico

library(dplyr)
library(ggplot2)

# Dataframe con la curva de ajuste
rango_concentracion <- seq(min(datos_2$Concentracion), max(datos_2$Concentracion), length = 100)
v0_pred <- (Vmax * rango_concentracion) / (Km + rango_concentracion)
predicciones <- data.frame(S = rango_concentracion, 
                           v0_pred = v0_pred)
concentraciones_breaks <- unique(datos_2$Concentracion)

#Generar grafico con ggplot
ggplot(datos_2, aes(Concentracion, v0)) +
         geom_line(data = predicciones, aes(x = S, y = v0_pred), size = 0.8, colour = "darkgray") +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult=1),
               geom="errorbar", color= "green3", width=0.0001, linewidth = 0.7) +
  stat_summary(fun.y=mean, geom="point", size = 3, color="navyblue") +
  scale_x_continuous(limits = c(0, 0.0015),
                     breaks = concentraciones_breaks,
                     labels = scales::number_format(accuracy = 0.0001)) +
  labs(title = "Ajuste no lineal del modelo de Michaelis-Menten",
              subtitle = paste0("Vmax =", round(Vmax, 6), "        Km =", round(Km, 6)),
              x = "Concentración de Sustrato (mM)", y = "Velocidad Inicial") +
         theme_light() + theme(plot.title = element_text(hjust = 0.5),
                                 plot.subtitle = element_text(hjust = 0.5),
                            axis.text.x = element_text(angle = 45, hjust = 1))

# Grafico sencillo
plot(modelo_ajustado, log = '', type = "all", 
     main = "Ajuste no lineal del modelo de Michaelis-Menten", 
     xlab = "Concentración de Sustrato (mM)", ylab = "Velocidad Inicial", pch = 19)

points(datos_2$Concentracion, y = y_predicho,
       col = "green3", pch = 19)

legend("bottomright", 
       legend = c("Puntos observados", "Puntos predichos"),
       pch = c(19, 19), 
       col = c("black", "green3"),
       bty = "n")
