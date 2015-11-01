/* Preprocessor variables and directives */
/* Para verificar lel correcto funcionamiento del termostato */
/* Guarda datos de velocidades para sacar histograma */
#undef CONTROL_TEMP
/* when not defined, better to do: #undef VAR */
#define THERM 1 /*0=NVE, no thermostat 1=langevin thermostat */
/* Defino variable para indicar si se graban las trayectorias de las partículas */
#undef GRABA_TRAYECTORIA
