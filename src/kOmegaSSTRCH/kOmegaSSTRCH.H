/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is NOT a part of OpenFOAM.

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


Class
    Foam::RASModels::kOmegaSSTRCH

Group
    grpRASTurbulence

Description
    Menter SST Two-Equation Model with Hellsten's Simplified Rotation/Curvature
    Correction (SST-RC-Hellsten) 
        https://turbmodels.larc.nasa.gov/sst.html

    This work was presented on the OpenCAE Symposium 2017@Nagoya, Japan.
    Nakayama, K. and Nakagawa, S. (2017)
    Implementation of k-omega SST model to rotation and curbature on OpenFOAM
        http://www.opencae.or.jp/activity/symposium/opencae_symposium2017/

    Codes are NOT tested throughout.

    References:
    \verbatim
    
        Hellsten, A. (1998).
        Some Improvements in Menter's k-omega SST Turbulence Model.
        AIAA Paper 98-2554.
        
        Mani, M., Ladd, J. A., and Bower, W. W. (2004).
        Rotation and Curvature Correction Assessment for One- and Two-Equation
        Turbulence Models.
        Journal of Aircraft, 41(2), 268-273.
    \endverbatim

    The model coefficients are
    \verbatim
        kOmegaSSTRCHCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1.0;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1.0;
            c1          10.0;
            F3          no;

            // Default RCH coefficients
            cRC         1.4;
        }
    \endverbatim

SourceFiles
    kOmegaSSTRCH.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSSTRCH_H
#define kOmegaSSTRCH_H

#include "kOmegaSST.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class kOmegaSSTRCH Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class kOmegaSSTRCH
:
    public kOmegaSST<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        kOmegaSSTRCH(const kOmegaSSTRCH&);
        void operator=(const kOmegaSSTRCH&);


protected:

    // Protected data

        // Model constants

            //- c_RC = 1.4
            dimensionedScalar cRC_;

        // Fields

            //- SQRT(2SijSij)
            volScalarField S_;

            //- SQRT(2WijWij)
            volScalarField W_;

            //- Richardson number
            volScalarField Ri_;


    // Protected Member Functions

        //- F4 = 1 / (1 + cRC*Ri)
        tmp<volScalarField> F4() const;

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("kOmegaSSTRCH");


    // Constructors

        //- Construct from components
        kOmegaSSTRCH
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~kOmegaSSTRCH()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "kOmegaSSTRCH.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
