/* Preprocessor variables and directives */
/* Para verificar lel correcto funcionamiento del termostato */
/* Guarda datos de velocidades para sacar histograma */
#undef CONTROL_TEMP
/* when not defined, better to do: #undef VAR */
#define THERM 0 /*0=NVE, no thermostat 1=langevin thermostat */
#define DEBUG 0
/* Defino variable para indicar si se graban las trayectorias de las partículas */
#define GRABA_TRAYECTORIA
#undef CORR_PAR
#define PABLO 0 /* 0=no se consideran agregados hechos por pablo */
#undef LUIS /*para el bloque de correlaciones porque no me compila*/
