/*=================================================*/
/* Preprocessor variables and directives           */
/*=================================================*/
/* Para verificar lel correcto funcionamiento del termostato */
/* Guarda datos de velocidades para sacar histograma */
#undef CONTROL_TEMP

/* when not defined, better to do: #undef VAR */
#define THERM 1 /*0=NVE, no thermostat 1=langevin thermostat */

/* Variable para indicar si se graban las trayectorias de las partículas */
#define GRABA_TRAYECTORIA

/* Variable para calcular la función de correlación entre pares g(r) */
#define CORR_PAR

/*==========================================*/
/* Variables para hacer cálculo vibracional */
#define MODOS_VIB
