#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CARTESIAN
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LimO3
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            3

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ENRG0                          0
#define  DNST0                          1
#define  GAMMA                          2

/* [Beg] user-defined constants (do not change this line) */

#define  INITIAL_SMOOTHING              YES

/* [End] user-defined constants (do not change this line) */
