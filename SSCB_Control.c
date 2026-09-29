#include <math.h>
#include <stddef.h>
#include <string.h>

#include "SSCB_Control.h"

/*
 * Portable C99 approximation of the SSCB_Control PSCAD component.
 *
 * Assumptions (the XML contains PSCAD blocks, not their runtime semantics):
 * - logic_mult Type=0 is AND and Type=1 is OR, as specified by the user.
 * - Settings are ordered as follows:
 *     V_Settings:   [OV limit (kV), OV delay (s), UV limit (kV),
 *                    UV delay (s), nominal voltage (kV)]
 *     TCC_Settings: [forward pickup (kA), forward pickup delay (s),
 *                    forward instantaneous pickup (kA), instantaneous delay (s),
 *                    reverse pickup (kA), reverse pickup delay (s),
 *                    reverse instantaneous pickup (kA), instantaneous delay (s)]
 * - VP2N2 is in kV; Ip and In are directional current magnitudes in kA.
 *   Current inputs are treated as magnitudes for their respective directions.
 * - Protection elements are high-side trips: the output is a closed/permissive
 *   signal (1) only while no trip is latched. Each pickup must persist for its
 *   configured delay; nonpositive pickup settings disable that element.
 * - The XML's detailed TCC energy/integral and reset behavior is represented
 *   here by pickup delays and instantaneous elements. Trip latches clear only
 *   when SSCB_Control_reset() is called.
 *
 * This is a behavioral starting point, not a bit-exact PSCAD code export.
 */

void SSCB_Control_reset(SSCB_ControlState *state)
{
    if (state != NULL) {
        memset(state, 0, sizeof(*state));
    }
}

static bool delay_elapsed(bool active, double delay_s, double dt_s,
                          double *timer_s)
{
    if (!active) {
        *timer_s = 0.0;
        return false;
    }
    if (delay_s <= 0.0) {
        *timer_s = 0.0;
        return true;
    }
    *timer_s += dt_s;
    return *timer_s >= delay_s;
}

static bool setting_is_finite(const double *values, size_t count)
{
    size_t i;
    for (i = 0; i < count; ++i) {
        if (!isfinite(values[i])) {
            return false;
        }
    }
    return true;
}

/*
 * Advance one sample and return the active-high PE_SW permissive output.
 * Invalid pointers, non-finite inputs/settings, or a negative/non-finite
 * timestep fail closed (return 0) without advancing the state.
 */
int SSCB_Control_step(SSCB_ControlState *state,
                      const SSCB_ControlInputs *inputs,
                      double dt_s)
{
    const double *v = inputs != NULL ? inputs->V_Settings : NULL;
    const double *tcc = inputs != NULL ? inputs->TCC_Settings : NULL;

    if (state == NULL || inputs == NULL || !isfinite(dt_s) || dt_s < 0.0 ||
        !isfinite(inputs->VP2N2) || !isfinite(inputs->Ip) ||
        !isfinite(inputs->In) ||
        !setting_is_finite(v, SSCB_V_SETTINGS_COUNT) ||
        !setting_is_finite(tcc, SSCB_TCC_SETTINGS_COUNT)) {
        return 0;
    }

    if (v[SSCB_V_OV_LIMIT] > 0.0) {
        state->ov_trip |= delay_elapsed(
            inputs->VP2N2 > v[SSCB_V_OV_LIMIT], v[SSCB_V_OV_DELAY], dt_s,
            &state->ov_timer_s);
    } else {
        state->ov_timer_s = 0.0;
    }

    if (v[SSCB_V_UV_LIMIT] > 0.0) {
        state->uv_trip |= delay_elapsed(
            inputs->VP2N2 < v[SSCB_V_UV_LIMIT], v[SSCB_V_UV_DELAY], dt_s,
            &state->uv_timer_s);
    } else {
        state->uv_timer_s = 0.0;
    }

    if (tcc[SSCB_TCC_FWD_PICKUP] > 0.0) {
        state->forward_trip |= delay_elapsed(
            fabs(inputs->Ip) > tcc[SSCB_TCC_FWD_PICKUP],
            tcc[SSCB_TCC_FWD_DELAY], dt_s, &state->forward_timer_s);
    } else {
        state->forward_timer_s = 0.0;
    }
    if (tcc[SSCB_TCC_FWD_INST_PICKUP] > 0.0) {
        state->forward_inst_trip |= delay_elapsed(
            fabs(inputs->Ip) > tcc[SSCB_TCC_FWD_INST_PICKUP],
            tcc[SSCB_TCC_FWD_INST_DELAY], dt_s,
            &state->forward_inst_timer_s);
    } else {
        state->forward_inst_timer_s = 0.0;
    }

    if (tcc[SSCB_TCC_REV_PICKUP] > 0.0) {
        state->reverse_trip |= delay_elapsed(
            fabs(inputs->In) > tcc[SSCB_TCC_REV_PICKUP],
            tcc[SSCB_TCC_REV_DELAY], dt_s, &state->reverse_timer_s);
    } else {
        state->reverse_timer_s = 0.0;
    }
    if (tcc[SSCB_TCC_REV_INST_PICKUP] > 0.0) {
        state->reverse_inst_trip |= delay_elapsed(
            fabs(inputs->In) > tcc[SSCB_TCC_REV_INST_PICKUP],
            tcc[SSCB_TCC_REV_INST_DELAY], dt_s,
            &state->reverse_inst_timer_s);
    } else {
        state->reverse_inst_timer_s = 0.0;
    }

    return !(state->ov_trip || state->uv_trip || state->forward_trip ||
             state->forward_inst_trip || state->reverse_trip ||
             state->reverse_inst_trip);
}
