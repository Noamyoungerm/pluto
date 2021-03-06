#define  PHYSICS                        MHD
#define  DIMENSIONS                     2
#define  GEOMETRY                       CYLINDRICAL
#define  BODY_FORCE                     NO
#define  COOLING                        SNEq
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  CHARACTERISTIC_TRACING
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            6

/* -- physics dependent declarations -- */

#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  DIVB_CONTROL                   NO
#define  BACKGROUND_FIELD               NO
#define  AMBIPOLAR_DIFFUSION            NO
#define  RESISTIVITY                    NO
#define  HALL_MHD                       NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  ETA                            0
#define  JET_VEL                        1
#define  SIGMA_Z                        2
#define  SIGMA_PHI                      3
#define  PERT_AMPLITUDE                 4
#define  PERT_PERIOD                    5

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  UNIT_DENSITY                   (CONST_amu*200.0)
#define  UNIT_LENGTH                    (2.5e15)
#define  UNIT_VELOCITY                  (1.e5)

/* [End] user-defined constants (do not change this line) */
