/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018-2022 OpenFOAM Foundation
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

#include "PB.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace diameterModels
    {
        namespace coalescenceModels
        {
            defineTypeNameAndDebug(PB, 0);
            addToRunTimeSelectionTable(
                coalescenceModel,
                PB,
                dictionary);
        }
    }
}

using Foam::constant::mathematical::pi;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::coalescenceModels::PB::
    PB(
        const populationBalanceModel &popBal,
        const dictionary &dict)
    : coalescenceModel(popBal, dict),
      LiaoNBase(popBal, dict),
      a_(dimensionedScalar::lookupOrDefault("a", dict, dimless, 1)),
      PMax_(dimensionedScalar::lookupOrDefault("PMax", dict, dimless, 0.8)),
      C1_(dimensionedScalar::lookupOrDefault("C1", dict, dimless, 0.356)),
      h0_(
          dimensionedScalar::lookupOrDefault(
              "h0",
              dict,
              dimLength,
              1e-4)),
      hf_(
          dimensionedScalar::lookupOrDefault(
              "hf",
              dict,
              dimLength,
              1e-8)),
      turbulence_(dict.lookup("turbulence")),
      buoyancy_(dict.lookup("buoyancy")),
      laminarShear_(dict.lookup("laminarShear")),
      eddyCapture_(dict.lookup("eddyCapture")),
      wakeEntrainment_(dict.lookup("wakeEntrainment")),
      CPack_(
          IOobject(
              "CPack",
              popBal_.time().timeName(),
              popBal_.mesh()),
          popBal_.mesh(),
          dimensionedScalar(
              "CPack",
              dimless,
              Zero)),
      CPackMax_(
          dimensionedScalar::lookupOrDefault("CPackMax", dict, dimless, 1e5)),
      dCrit_(
          IOobject(
              "dCrit",
              popBal_.time().timeName(),
              popBal_.mesh()),
          popBal_.mesh(),
          dimensionedScalar(
              "dCrit",
              dimLength,
              Zero)),
      uRelBuoy_(
          IOobject(
              "uRelBuoy",
              popBal_.time().timeName(),
              popBal_.mesh()),
          popBal_.mesh(),
          dimensionedScalar(
              "uRelBuoy",
              dimVelocity,
              Zero)),
      correctionFunctionType_(
          dict.lookupOrDefault<word>("correctionFunction", "none"))
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diameterModels::coalescenceModels::PB::precompute()
{
    LiaoNBase::precompute();

    CPack_ = min(PMax_ / max(PMax_ - popBal_.alphas(), SMALL), CPackMax_);

    const uniformDimensionedVectorField &g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    dCrit_ =
        4 * sqrt(
                popBal_.sigmaWithContinuousPhase(popBal_.sizeGroups()[1].phase())() / (mag(g) * (popBal_.continuousPhase().rho() - popBal_.sizeGroups()[1].phase().rho()())));
}

void Foam::diameterModels::coalescenceModels::PB::
    addToCoalescenceRate(
        volScalarField &coalescenceRate,
        const label i,
        const label j)
{
    const phaseModel &continuousPhase = popBal_.continuousPhase();
    const sizeGroup &fi = popBal_.sizeGroups()[i];
    const sizeGroup &fj = popBal_.sizeGroups()[j];
    const uniformDimensionedVectorField &g =
        popBal_.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const dimensionedScalar rij(1 / (1 / fi.dSph() + 1 / fj.dSph()));
    dimensionedScalar Aij(pi * 0.25 * sqr(fi.dSph() + fj.dSph()));
    const dimensionedScalar minEpsilon_("minEpsilon", dimensionSet(0, 2, -3, 0, 0, 0, 0), 1e-20);

    const volScalarField rd_(

        pow(30, 0.75) * pow(popBal_.continuousPhase().thermo().nu(), 3 / 4) * pow(max(popBal_.continuousTurbulence().epsilon(), minEpsilon_), -0.25));
    const volScalarField rl_(

        pow(2 / 3, 1.5) * pow(popBal_.continuousTurbulence().k(), 1.5) / max(popBal_.continuousTurbulence().epsilon(), minEpsilon_));

    const volScalarField vsqri(
        2 * pow(popBal_.continuousPhase().thermo().nu(), 2) *
        pow(pow(fi.dSph(), 2) / (pow(fi.dSph(), 2) + pow(rd_, 2)), 2.0 / 3.0) *
        pow(pow(fi.dSph(), 2) / (pow(fi.dSph(), 2) + pow(rl_, 2)), 1 / 3));
    const volScalarField vsqrj(
        2 * pow(popBal_.continuousPhase().thermo().nu(), 2) *
        pow(pow(fj.dSph(), 2) / (pow(fj.dSph(), 2) + pow(rd_, 2)), 2.0 / 3.0) *
        pow(pow(fj.dSph(), 2) / (pow(fj.dSph(), 2) + pow(rl_, 2)), 1 / 3));
    const volScalarField hinit_(
        0.014 * rij * continuousPhase.rho() * sqrt(vsqri + vsqrj) / popBal_.sigmaWithContinuousPhase(fi.phase()));
    const volScalarField collisionEfficiency(
        exp(
            -sqrt(
                pow3(rij) * continuousPhase.rho() / (16 * popBal_.sigmaWithContinuousPhase(fi.phase()))) *
            log(h0_ / hf_) * cbrt(popBal_.continuousTurbulence().epsilon()) / pow(rij, 2 / 3)));

    if (turbulence_)
    {
        coalescenceRate +=
            (C1_ * pi * sqr(fi.dSph() + fj.dSph()) * cbrt(vsqri + vsqrj)) * collisionEfficiency;
    }

    if (buoyancy_)
    {
        const dimensionedScalar Sij(pi / 4 * sqr(fi.dSph() + fj.dSph()));

        coalescenceRate +=
            (Sij * mag(
                       sqrt(
                           2.14 * popBal_.sigmaWithContinuousPhase(fi.phase()) / (continuousPhase.rho() * fi.dSph()) + 0.505 * mag(g) * fi.dSph()) -
                       sqrt(
                           2.14 * popBal_.sigmaWithContinuousPhase(fi.phase()) / (continuousPhase.rho() * fj.dSph()) + 0.505 * mag(g) * fj.dSph()))) *
            collisionEfficiency;
    }

    if (laminarShear_)
    {
        coalescenceRate +=
            pow3(fi.dSph() + fj.dSph()) / 6 * shearStrainRate_ * collisionEfficiency;
    }

    if (eddyCapture_)
    {
        const volScalarField uRelEddy(
            1 * 0.5 / pi * (fi.dSph() + fj.dSph()) * eddyStrainRate_);

        coalescenceRate +=
            pos0(kolmogorovLengthScale_ - (fi.dSph() + fj.dSph())) * CPack_ * 0.5 * pi * 0.25 * sqr(fi.dSph() + fj.dSph()) * uRelEddy * collisionEfficiency;
    }

    if (wakeEntrainment_)
    {
        const dimensionedScalar uRelWakeI(
            1 * uTerminal_[i] * cbrt(Cd_[i]));

        const dimensionedScalar uRelWakeJ(
            1 * uTerminal_[j] * cbrt(Cd_[j]));

        coalescenceRate +=
            CPack_ * 0.125 * pi * (sqr(fi.dSph()) * uRelWakeI * pos0(fi.dSph() - 0.5 * dCrit_) * (pow6(fi.dSph() - 0.5 * dCrit_) / (pow6(fi.dSph() - 0.5 * dCrit_) + pow6(0.5 * dCrit_))) + sqr(fj.dSph()) * uRelWakeJ * pos0(fj.dSph() - 0.5 * dCrit_) * (pow6(fj.dSph() - 0.5 * dCrit_) / (pow6(fj.dSph() - 0.5 * dCrit_) + pow6(0.5 * dCrit_))));
    }

    coalescenceRate *= a_;
}

// ************************************************************************* //
