VERSION:
3.001

#include "SSCB_Control.h"

STATIC:

	double dt;
	double ovTimer, uvTimer;
	double fwdTimer, fwdInstTimer;
	double revTimer, revInstTimer;
	int ovTrip, uvTrip;
	int fwdTrip, fwdInstTrip;
	int revTrip, revInstTrip;

RAM_FUNCTIONS:
	static int delayElapsed(int active, double delay, double step, double *timer)
	{
		if (!active)
		{
			*timer = 0.0;
			return 0;
		}

		if (delay <= 0.0)
		{
			*timer = 0.0;
			return 1;
		}

		*timer += step;
		return (*timer >= delay);
	}

RAM:
	dt = getTimeStep();

	ovTimer = 0.0;
	uvTimer = 0.0;
	fwdTimer = 0.0;
	fwdInstTimer = 0.0;
	revTimer = 0.0;
	revInstTimer = 0.0;

	ovTrip = 0;
	uvTrip = 0;
	fwdTrip = 0;
	fwdInstTrip = 0;
	revTrip = 0;
	revInstTrip = 0;

CODE_FUNCTIONS:
#include <builtin_gcc.h>

CODE:
	/*
	 * Assumptions for the approximation:
	 * V_Settings = [OV limit (kV), OV delay (s), UV limit (kV),
	 *               UV delay (s), nominal voltage (kV)].
	 * TCC_Settings = [forward pickup (kA), forward delay (s),
	 *                 forward instantaneous pickup (kA), instantaneous delay (s),
	 *                 reverse pickup (kA), reverse delay (s),
	 *                 reverse instantaneous pickup (kA), instantaneous delay (s)].
	 * Ip and In are directional current magnitudes. Nonpositive pickup limits
	 * disable that element. Trips latch until the simulation is reinitialized.
	 * This approximates the PSCAD logic; it does not reproduce the detailed
	 * TCC energy/integral behavior.
	 */

	dt = getTimeStep();

	if (V_Settings[0] > 0.0)
	{
		ovTrip |= delayElapsed(VP2N2 > V_Settings[0],
			V_Settings[1], dt, &ovTimer);
	}
	else
	{
		ovTimer = 0.0;
	}

	if (V_Settings[2] > 0.0)
	{
		uvTrip |= delayElapsed(VP2N2 < V_Settings[2],
			V_Settings[3], dt, &uvTimer);
	}
	else
	{
		uvTimer = 0.0;
	}

	if (TCC_Settings[0] > 0.0)
	{
		fwdTrip |= delayElapsed(fabs(Ip) > TCC_Settings[0],
			TCC_Settings[1], dt, &fwdTimer);
	}
	else
	{
		fwdTimer = 0.0;
	}

	if (TCC_Settings[2] > 0.0)
	{
		fwdInstTrip |= delayElapsed(fabs(Ip) > TCC_Settings[2],
			TCC_Settings[3], dt, &fwdInstTimer);
	}
	else
	{
		fwdInstTimer = 0.0;
	}

	if (TCC_Settings[4] > 0.0)
	{
		revTrip |= delayElapsed(fabs(In) > TCC_Settings[4],
			TCC_Settings[5], dt, &revTimer);
	}
	else
	{
		revTimer = 0.0;
	}

	if (TCC_Settings[6] > 0.0)
	{
		revInstTrip |= delayElapsed(fabs(In) > TCC_Settings[6],
			TCC_Settings[7], dt, &revInstTimer);
	}
	else
	{
		revInstTimer = 0.0;
	}

	/* PE_SW is high while the breaker is permissive and no trip has latched. */
	PE_SW = !(ovTrip || uvTrip || fwdTrip || fwdInstTrip ||
		revTrip || revInstTrip);
