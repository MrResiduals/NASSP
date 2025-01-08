/***************************************************************************
  This file is part of Project Apollo - NASSP

  MCC Calculations

  Project Apollo is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Project Apollo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with Project Apollo; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  See http://nassp.sourceforge.net/license/ for more details.

  **************************************************************************/

#include "MCC_Calculations.h"
#include "rtcc.h"

MCC_Calculations::MCC_Calculations(RTCC *r) : RTCCModule(r)
{

}

bool MCC_Calculations::CreateEphemeris(EphemerisData sv, double EphemerisLeftLimitGMT, double EphemerisRightLimitGMT, EphemerisDataTable2 &ephem)
{
	EMSMISSInputTable in;
	PLAWDTOutput weights;

	in.AnchorVector = sv;
	in.EphemerisLeftLimitGMT = EphemerisLeftLimitGMT;
	in.EphemerisRightLimitGMT = EphemerisRightLimitGMT;

	in.EphemerisBuildIndicator = true;
	if (sv.RBI == BODY_EARTH)
	{
		in.ECIEphemerisIndicator = true;
		in.ECIEphemTableIndicator = &ephem;
	}
	else
	{
		in.MCIEphemerisIndicator = true;
		in.MCIEphemTableIndicator = &ephem;
	}
	in.useInputWeights = true;

	//TBD
	weights.KFactor = 0.0;
	weights.CC.set(0);
	weights.CSMWeight = 1.0;

	in.DensityMultiplier = 0.0;
	in.WeightsTable = &weights;
	in.VehicleCode = RTCC_MPT_CSM; //Not used

	pRTCC->EMSMISS(&in);

	ephem.Header.TUP = 1; //Only has to be non-zero

	return (in.NIAuxOutputTable.ErrorCode != 0);
}
double MCC_Calculations::EnvironmentChange(EphemerisDataTable2 &ephem, double gmt_estimate, int option, bool present, bool terminator)
{
	ManeuverTimesTable MANTIMES;
	LunarStayTimesTable LUNRSTAY;
	RTCC::EMMENVInputTable in;
	RTCC::EMMENVOutputTable out;

	in.GMT = gmt_estimate;
	in.option = option;
	in.present = present;
	in.terminator = terminator;

	pRTCC->EMMENV(ephem, MANTIMES, &LUNRSTAY, in, out);

	return out.T_Change;
}

double MCC_Calculations::Sunrise(EphemerisDataTable2 &ephem, double gmt_estimate)
{
	return EnvironmentChange(ephem, gmt_estimate, 0, true, false);
}

double MCC_Calculations::TerminatorRise(EphemerisDataTable2 &ephem, double gmt_estimate)
{
	return EnvironmentChange(ephem, gmt_estimate, 0, true, true);
}

bool MCC_Calculations::LongitudeCrossing(EphemerisDataTable2 &ephem, double lng, double gmt_estimate, double &gmt_cross)
{
	ManeuverTimesTable MANTIMES;
	EphemerisDataTable2 ephem_true;
	EphemerisData2 sv;
	double dErr;
	int in, out, iErr;

	gmt_cross = 0.0;

	//Convert to ECT or MCT table
	ephem_true = ephem;

	in = ephem.Header.CSI;
	out = in + 1;

	iErr = pRTCC->ELVCNV(ephem.table, in, out, ephem_true.table);
	if (iErr) return true;

	ephem_true.Header.CSI = out;

	//Longitude crossing
	dErr = pRTCC->RLMTLC(ephem_true, MANTIMES, lng, gmt_estimate, gmt_cross, sv);

	return (dErr < 0.0);
}