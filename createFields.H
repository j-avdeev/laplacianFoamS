    Info<< "Reading field T\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );
    
    
//     Info<< "Creating field F\n" << endl;
//     volScalarField F
//     (
// 	IOobject
// 	(
// 	    "F",
// 	    runTime.timeName(),
// 	    mesh,
// 	    IOobject::MUST_READ,
// 	    IOobject::AUTO_WRITE
// 	),
// 	mesh
//     );
    
    Info<< "Creating field S\n" << endl;
    volScalarField S
    (
	IOobject
	(
	    "S",
	    runTime.timeName(),
	    mesh,
	    IOobject::MUST_READ,
	    IOobject::AUTO_WRITE
	),
	mesh
    );