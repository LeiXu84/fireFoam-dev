    volScalarField T0
    (
        IOobject
        (
            "T0",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("T0", dimTemperature, scalar(298.15))
    );

    volScalarField Tlink1
    (
        IOobject
        (
            "Tlink1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tlink1", dimTemperature, scalar(298.15))
    );

    volScalarField tAct1
    (
        IOobject
        (
            "tAct1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tAct1", dimTime, scalar(1e3))
    );

    volScalarField Tlink2
    (
        IOobject
        (
            "Tlink2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tlink2", dimTemperature, scalar(298.15))
    );

    volScalarField tAct2
    (
        IOobject
        (
            "tAct2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tAct2", dimTime, scalar(1e3))
    );

    volScalarField Tlink3
    (
        IOobject
        (
            "Tlink3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tlink3", dimTemperature, scalar(298.15))
    );

    volScalarField tAct3
    (
        IOobject
        (
            "tAct3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tAct3", dimTime, scalar(1e3))
    );

    volScalarField Tlink4
    (
        IOobject
        (
            "Tlink4",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tlink4", dimTemperature, scalar(298.15))
    );

    volScalarField tAct4
    (
        IOobject
        (
            "tAct4",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tAct4", dimTime, scalar(1e3))
    );

    volScalarField Tlink5
    (
        IOobject
        (
            "Tlink5",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Tlink5", dimTemperature, scalar(298.15))
    );

    volScalarField tAct5
    (
        IOobject
        (
            "tAct5",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tAct5", dimTime, scalar(1e3))
    );


// RTI dictionary
IOdictionary rtiPropertiesDict
(
    IOobject
    (
        "rtiProperties",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT, //MUST_READ,
        IOobject::NO_WRITE
    )
);
Switch rtiActive
(
    rtiPropertiesDict.lookupOrDefault<Switch>("active",false)
);
scalar RTI1_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("RTI1",1.0)
);
scalar C1_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("C1",1.0)
);
scalar Tact1_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("Tact1",298.0)
);
scalar RTI2_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("RTI2",1.0)
);
scalar C2_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("C2",1.0)
);
scalar Tact2_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("Tact2",298.0)
);
scalar RTI3_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("RTI3",1.0)
);
scalar C3_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("C3",1.0)
);
scalar Tact3_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("Tact3",298.0)
);
scalar RTI4_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("RTI4",1.0)
);
scalar C4_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("C4",1.0)
);
scalar Tact4_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("Tact4",298.0)
);
scalar RTI5_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("RTI5",1.0)
);
scalar C5_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("C5",1.0)
);
scalar Tact5_
(
    rtiPropertiesDict.lookupOrDefault<scalar>("Tact5",298.0)
);

    volScalarField dTeOld1
    (
        IOobject
        (
            "dTeOld1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dTeOld1", dimless, scalar(0))
    );
    volScalarField dTeOld2
    (
        IOobject
        (
            "dTeOld2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dTeOld2", dimless, scalar(0))
    );
    volScalarField dTeOld3
    (
        IOobject
        (
            "dTeOld3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dTeOld3", dimless, scalar(0))
    );
    volScalarField dTeOld4
    (
        IOobject
        (
            "dTeOld4",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dTeOld4", dimless, scalar(0))
    );
    volScalarField dTeOld5
    (
        IOobject
        (
            "dTeOld5",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("dTeOld5", dimless, scalar(0))
    );
    volScalarField tActOld1
    (
        IOobject
        (
            "tActOld1",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tActOld1", dimTime, scalar(1e3))
    );
    volScalarField tActOld2
    (
        IOobject
        (
            "tActOld2",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tActOld2", dimTime, scalar(1e3))
    );
    volScalarField tActOld3
    (
        IOobject
        (
            "tActOld3",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tActOld3", dimTime, scalar(1e3))
    );
    volScalarField tActOld4
    (
        IOobject
        (
            "tActOld4",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tActOld4", dimTime, scalar(1e3))
    );
    volScalarField tActOld5
    (
        IOobject
        (
            "tActOld5",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("tActOld5", dimTime, scalar(1e3))
    );
dimensionedScalar dimT("dimT",dimTemperature,1);
dimensionedScalar dimU("dimU",dimTime/dimLength,1);
