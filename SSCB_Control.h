#ifndef SSCB_CONTROL_H
#define SSCB_CONTROL_H

#include <stdbool.h>

enum {
    SSCB_V_OV_LIMIT = 0,
    SSCB_V_OV_DELAY = 1,
    SSCB_V_UV_LIMIT = 2,
    SSCB_V_UV_DELAY = 3,
    SSCB_V_NOMINAL = 4,
    SSCB_V_SETTINGS_COUNT = 5,

    SSCB_TCC_FWD_PICKUP = 0,
    SSCB_TCC_FWD_DELAY = 1,
    SSCB_TCC_FWD_INST_PICKUP = 2,
    SSCB_TCC_FWD_INST_DELAY = 3,
    SSCB_TCC_REV_PICKUP = 4,
    SSCB_TCC_REV_DELAY = 5,
    SSCB_TCC_REV_INST_PICKUP = 6,
    SSCB_TCC_REV_INST_DELAY = 7,
    SSCB_TCC_SETTINGS_COUNT = 8
};

typedef struct {
    double VP2N2;
    double Ip;
    double In;
    double V_Settings[SSCB_V_SETTINGS_COUNT];
    double TCC_Settings[SSCB_TCC_SETTINGS_COUNT];
} SSCB_ControlInputs;

typedef struct {
    double ov_timer_s;
    double uv_timer_s;
    double forward_timer_s;
    double forward_inst_timer_s;
    double reverse_timer_s;
    double reverse_inst_timer_s;
    bool ov_trip;
    bool uv_trip;
    bool forward_trip;
    bool forward_inst_trip;
    bool reverse_trip;
    bool reverse_inst_trip;
} SSCB_ControlState;

void SSCB_Control_reset(SSCB_ControlState *state);
int SSCB_Control_step(SSCB_ControlState *state,
                      const SSCB_ControlInputs *inputs,
                      double dt_s);

#endif
