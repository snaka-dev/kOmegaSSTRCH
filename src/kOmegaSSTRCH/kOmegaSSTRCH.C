/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

#include "kOmegaSSTRCH.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> kOmegaSSTRCH<BasicTurbulenceModel>::F4() const
{
    return 1 / (1 + cRC_*Ri_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
kOmegaSSTRCH<BasicTurbulenceModel>::kOmegaSSTRCH
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kOmegaSST<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
	, type // added for test
    ),

    cRC_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "cRC",
            this->coeffDict_,
            1.4
        )
    ),
    S_
    (
        IOobject
        (
            "S",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ, //READ_IF_PRESENT
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("S",dimless/dimTime,0.)
    ),
    W_
    (
        IOobject
        (
            "W",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ, //READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("W",dimless/dimTime,0.)
    ),
    Ri_
    (
        IOobject
        (
            "Ri", 
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ, //READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,  //scalar(0.)
        dimensionedScalar("Ri",dimless,scalar(0.0)) 
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool kOmegaSSTRCH<BasicTurbulenceModel>::read()
{
    if (kOmegaSST<BasicTurbulenceModel>::read())
    {
        cRC_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
void kOmegaSSTRCH<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    fv::options& fvOptions(fv::options::New(this->mesh_));
    
    volScalarField& k_ = this->k_;
    volScalarField& omega_ = this->omega_;
    const dimensionedScalar& alphaOmega2_ = this->alphaOmega2_;
    
    BasicTurbulenceModel::correct();

    volScalarField::Internal divU(fvc::div(fvc::absolute(this->phi(), U)));

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField S2(2*magSqr(symm(tgradU())));
    volScalarField::Internal GbyNu0((tgradU() && dev(twoSymm(tgradU()))));
    volScalarField::Internal G(this->GName(), nut*GbyNu0);

    W_ = sqrt(2*magSqr(skew(tgradU())));
    S_ = sqrt(S2);
    Ri_ = W_/(S_+dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL))
         *(W_/(S_+dimensionedScalar ("VSMALL", dimless/dimTime, VSMALL))-1);
	 
    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField F1(this->F1(CDkOmega));
    volScalarField F23(this->F23());
    volScalarField F4(this->F4());
    
    {
        volScalarField::Internal gamma(this->gamma(F1));
        volScalarField::Internal beta(this->beta(F1));

        // Turbulent frequency equation
        tmp<fvScalarMatrix> omegaEqn
        (
            fvm::ddt(alpha, rho, omega_)
          + fvm::div(alphaRhoPhi, omega_)
          - fvm::laplacian(alpha*rho*this->DomegaEff(F1), omega_)
         ==
            alpha()*rho()*gamma*this->GbyNu(GbyNu0, F23(), S2())
          - fvm::SuSp((2.0/3.0)*alpha()*rho()*gamma*divU, omega_)
          //- fvm::Sp(alpha()*rho()*beta*omega_(), omega_)
          - fvm::Sp(F4*alpha()*rho()*beta*omega_(), omega_)    // F4 is included.
          - fvm::SuSp
            (
                alpha()*rho()*(F1() - scalar(1))*CDkOmega()/omega_(),
                omega_
            )
          + this->Qsas(S2(), gamma, beta)
          + this->omegaSource()
          + fvOptions(alpha, rho, omega_)
        );

        omegaEqn.ref().relax();
        fvOptions.constrain(omegaEqn.ref());
        omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
        solve(omegaEqn);
        fvOptions.correct(omega_);
        bound(omega_, this->omegaMin_);
    }

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*this->DkEff(F1), k_)
     ==
        alpha()*rho()*this->Pk(G)
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*this->epsilonByk(F1, tgradU()), k_)
      + this->kSource()
      + fvOptions(alpha, rho, k_)
    );

    tgradU.clear();

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    this->correctNut(S2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
