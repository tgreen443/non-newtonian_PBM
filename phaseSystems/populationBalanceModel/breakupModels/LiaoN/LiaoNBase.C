/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "LiaoNBase.H"
#include "fvcGrad.H"
#include "dispersedPhaseInterface.H"
#include "dispersedDragModel.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::LiaoNBase::LiaoNBase
(
    const populationBalanceModel& popBal,
    const dictionary& dict
)
:
    populationBalance_(popBal),
    kolmogorovLengthScale_
    (
        IOobject
        (
            "kolmogorovLengthScale",
            populationBalance_.time().timeName(),
            populationBalance_.mesh()
        ),
        populationBalance_.mesh(),
        dimensionedScalar
        (
            "kolmogorovLengthScale",
            dimLength,
            Zero
        )
    ),
    shearStrainRate_
    (
        IOobject
        (
            "shearStrainRate",
            populationBalance_.time().timeName(),
            populationBalance_.mesh()
        ),
        populationBalance_.mesh(),
        dimensionedScalar
        (
            "shearStrainRate",
            dimVelocity/dimLength,
            Zero
        )
    ),
    eddyStrainRate_
    (
        IOobject
        (
            "eddyStrainRate",
            populationBalance_.time().timeName(),
            populationBalance_.mesh()
        ),
        populationBalance_.mesh(),
        dimensionedScalar
        (
            "eddyStrainRate",
            dimVelocity/dimLength,
            Zero
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::LiaoNBase::precompute()
{
    const dimensionedScalar minEpsilon_
    (
        "minEpsilon",
        dimensionSet(0,2,-3,0,0,0,0),
        1e-9
    );
    kolmogorovLengthScale_ =
        pow025
        (
            pow3(populationBalance_.continuousPhase().thermo().nu())
           /(populationBalance_.continuousTurbulence().epsilon(),minEpailon_)
        );

    shearStrainRate_ =
        sqrt(2.0)
       *mag(symm(fvc::grad(populationBalance_.continuousPhase().U())));

    eddyStrainRate_ =
           sqrt
           (
               populationBalance_.continuousPhase().rho()
              *populationBalance_.continuousTurbulence().epsilon()
              /populationBalance_.continuousPhase().thermo().mu()
           );

    if (uTerminal_.empty())
    {
        const fvMesh& mesh = populationBalance_.mesh();
        const uniformDimensionedVectorField& g =
            mesh.lookupObject<uniformDimensionedVectorField>("g");

        const dimensionedScalar nuc
        (
            "nuc",
            dimViscosity,
            gAverage(populationBalance_.continuousPhase().thermo().nu()())
        );

        const dimensionedScalar rhoc
        (
            "rhoc",
            dimDensity,
            gAverage(populationBalance_.continuousPhase().rho()())
        );

        const dimensionedScalar rhod
        (
            "rhod",
            dimDensity,
            gAverage(populationBalance_.sizeGroups()[1].phase().rho()())
        );

        const dimensionedScalar sigma
        (
            "sigma",
            dimForce/dimLength,
            gAverage
            (
                populationBalance_.sigmaWithContinuousPhase
                (
                    populationBalance_.sizeGroups()[1].phase()
                )()
            )
        );

        for(int m = 0; m < populationBalance_.sizeGroups().size(); m++)
        {
            const sizeGroup& f = populationBalance_.sizeGroups()[m];

           // Construct dispersedPhaseInterface for this size group
            const dispersedPhaseInterface interface
            (
                f.phase(), // dispersed phase
                populationBalance_.continuousPhase() // continuous phase
            );

            // Lookup drag model for the interface
            if (!mesh.foundObject<dragModel>(IOobject::groupName(dragModel::typeName, interface.name())))
            {
                FatalErrorInFunction
                    << "No drag model found for interface " << interface.name()
                    << exit(FatalError);
            }

            const dragModels::dispersedDragModel& drag =
                mesh.lookupObject<dragModels::dispersedDragModel>
                (
                    IOobject::groupName(dragModel::typeName, interface.name())
                );
            dimensionedScalar uTerminal("uTerminal", dimVelocity, 0.2);
            dimensionedScalar Cd("Cd", dimless, 0.44);
            dimensionedScalar CdEllipse("CdEllipse", dimless, 1);

            dimensionedScalar Re(uTerminal*f.dSph()/nuc);
            const dimensionedScalar Eo
            (
                mag(g)*mag(rhoc - rhod)*sqr(f.dSph())/sigma
            );

            dimensionedScalar F("F", dimForce/dimArea, 1);
            dimensionedScalar dF("dF", dimForce/dimArea/dimVelocity, 1);
            const dimensionedScalar uTerminalX("uTerminalX", dimVelocity, 1e-5);
            dimensionedScalar ReX("ReX", dimless, Re.value());
            dimensionedScalar CdX("CdX", dimless, Cd.value());
            dimensionedScalar dCd("dCd", Cd.dimensions()/dimVelocity, Zero);

            int n = 0;

            while(mag(F.value()) >= 1.0e-05 && n++ <= 20)
            {
                Re = uTerminal*f.dSph()/nuc;

                Cd =dimensionedScalar
                (
                    "Cd",
                    dimless,
                    gAverage(drag.Cd(Re)())
                );
                    

                F =
                    4.0/3.0*(rhoc - rhod)*mag(g)*f.dSph()
                  - rhoc*Cd*sqr(uTerminal);

                ReX = (uTerminal + uTerminalX)*f.dSph()/nuc;

                CdX = dimensionedScalar
                (
                    "CdX",
                    dimless,
                    gAverage(drag.Cd(ReX)())
                );

                dCd = (CdX - Cd) / uTerminalX;

                dF = -(2*rhoc*uTerminal*Cd + rhoc*sqr(uTerminal)*dCd);

                uTerminal -= F/dF;
            }

            uTerminal_.append(new dimensionedScalar("uTerminal", uTerminal));

            Cd_.append(new dimensionedScalar("Cd", Cd));
        }
    }
}


// ************************************************************************* //
