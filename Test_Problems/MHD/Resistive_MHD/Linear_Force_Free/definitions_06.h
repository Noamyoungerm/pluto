#define  PHYSICS                        MHD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        0
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            1

/* -- physics dependent declarations -- */

#define  EOS                            ISOTHERMAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   CONSTRAINED_TRANSPORT
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    EXPLICIT
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  DENSITY                        0

/* [Beg] user-defined constants (do not change this line) */

#define  INCLUDE_JDIR                   NO
#define  LIMITER                        MC_LIM
#define  CT_EMF_AVERAGE                 ARITHMETIC
#define  CHECK_DIVB_CONDITION           TRUE

/* [End] user-defined constants (do not change this line) */
