/* Preprocessor variables and directives */
/* Para verificar lel correcto funcionamiento del termostato */
/* Guarda datos de velocidades para sacar histograma */
#undef CONTROL_TEMP
/* when not defined, better to do: #undef VAR */
#define THERM 1 /*0=NVE, no thermostat 1=langevin thermostat */
#define DEBUG 0
/* Defino variable para indicar si se graban las trayectorias de las partículas */
#define GRABA_TRAYECTORIA
#define CORR_PAR
#define PABLO 1 /* 0=no se consideran agregados hechos por pablo */
#undef LUIS /*para el bloque de correlaciones porque no me compila*/
#undef MODOS_VIB
#undef MODOS_VIB_EQUIVALENTES
