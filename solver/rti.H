if(rtiActive)
{
	const volScalarField& Tgas = T;
	volScalarField dTg = Tgas - T0;
	const volScalarField Ugas = mag(U);
	const scalar deltaT =  mesh.time().deltaTValue();
	
	//RTI1
	volScalarField dTe1 = sqrt(Ugas*dimU)/RTI1_*(dTg/dimT-(1+C1_/(sqrt(Ugas*dimU)+SMALL))*dTeOld1)*deltaT+dTeOld1;
	Tlink1 = T0 + dTe1*dimT;
	dTeOld1 = dTe1;
	forAll (tAct1, cellI) {
		if ( Tlink1[cellI] >= Tact1_ ) {
			tAct1[cellI] = runTime.value();
		}
		if ( pos(tActOld1[cellI]-tAct1[cellI]) == 0) {
			tAct1[cellI] = tActOld1[cellI];
		}
	}
	tActOld1 = tAct1;
	//RTI2
	volScalarField dTe2 = sqrt(Ugas*dimU)/RTI2_*(dTg/dimT-(1+C2_/(sqrt(Ugas*dimU)+SMALL))*dTeOld2)*deltaT+dTeOld2;
	Tlink2 = T0 + dTe2*dimT;
	dTeOld2 = dTe2;
	forAll (tAct2, cellI) {
		if ( Tlink2[cellI] >= Tact2_ ) {
			tAct2[cellI] = runTime.value();
		}
		if ( pos(tActOld2[cellI]-tAct2[cellI]) == 0) {
			tAct2[cellI] = tActOld2[cellI];
		}
	}
	tActOld2 = tAct2;
	//RTI3
	volScalarField dTe3 = sqrt(Ugas*dimU)/RTI3_*(dTg/dimT-(1+C3_/(sqrt(Ugas*dimU)+SMALL))*dTeOld3)*deltaT+dTeOld3;
	Tlink3 = T0 + dTe3*dimT;
	dTeOld3 = dTe3;
	forAll (tAct3, cellI) {
		if ( Tlink3[cellI] >= Tact3_ ) {
			tAct3[cellI] = runTime.value();
		}
		if ( pos(tActOld3[cellI]-tAct3[cellI]) == 0) {
			tAct3[cellI] = tActOld3[cellI];
		}
	}
	tActOld3 = tAct3;
	//RTI4
	volScalarField dTe4 = sqrt(Ugas*dimU)/RTI4_*(dTg/dimT-(1+C4_/(sqrt(Ugas*dimU)+SMALL))*dTeOld4)*deltaT+dTeOld4;
	Tlink4 = T0 + dTe4*dimT;
	dTeOld4 = dTe4;
	forAll (tAct4, cellI) {
		if ( Tlink4[cellI] >= Tact4_ ) {
			tAct4[cellI] = runTime.value();
		}
		if ( pos(tActOld4[cellI]-tAct4[cellI]) == 0) {
			tAct4[cellI] = tActOld4[cellI];
		}
	}
	tActOld4 = tAct4;
	//RTI5
	volScalarField dTe5 = sqrt(Ugas*dimU)/RTI5_*(dTg/dimT-(1+C5_/(sqrt(Ugas*dimU)+SMALL))*dTeOld5)*deltaT+dTeOld5;
	Tlink5 = T0 + dTe5*dimT;
	dTeOld5 = dTe5;
	forAll (tAct5, cellI) {
		if ( Tlink5[cellI] >= Tact5_ ) {
			tAct5[cellI] = runTime.value();
		}
		if ( pos(tActOld5[cellI]-tAct5[cellI]) == 0) {
			tAct5[cellI] = tActOld5[cellI];
		}
	}
	tActOld5 = tAct5;
}
