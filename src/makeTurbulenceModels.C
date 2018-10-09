#include "IncompressibleTurbulenceModel.H"
#include "incompressible/transportModel/transportModel.H"
#include "addToRunTimeSelectionTable.H"

#include "laminarModel.H"
// above line can be leplaced with #include "fvm.H" and #include "fvc.H"
#include "RASModel.H"

#include "fvOptions.H"

// from makeTurbulenceModel.H 

#define makeTurbulenceModelTypes(Alpha, Rho, baseModel, BaseModel, Transport)  \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef BaseModel<Transport> Transport##BaseModel;                     \
        typedef RASModel<Transport##BaseModel> RAS##Transport##BaseModel;      \
    }

#define makeBaseTurbulenceModel(Alpha, Rho, baseModel, BaseModel, Transport)   \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        typedef TurbulenceModel                                                \
        <                                                                      \
            Alpha,                                                             \
            Rho,                                                               \
            baseModel,                                                         \
            Transport                                                          \
        > Transport##baseModel;                                                \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            Transport##baseModel,                                              \
            dictionary                                                         \
        );                                                                     \
                                                                               \
                                                                               \
        defineNamedTemplateTypeNameAndDebug(RAS##Transport##BaseModel, 0);     \
                                                                               \
        defineTemplateRunTimeSelectionTable                                    \
        (RAS##Transport##BaseModel, dictionary);                               \
                                                                               \
                                                                               \
    }

#define makeTemplatedTurbulenceModel(BaseModel, SType, Type)                   \
    defineNamedTemplateTypeNameAndDebug                                        \
        (Foam::SType##Models::Type<Foam::BaseModel>, 0);                       \
                                                                               \
    namespace Foam                                                             \
    {                                                                          \
        namespace SType##Models                                                \
        {                                                                      \
            typedef Type<BaseModel> Type##SType##BaseModel;                    \
                                                                               \
            addToRunTimeSelectionTable                                         \
            (                                                                  \
                SType##BaseModel,                                              \
                Type##SType##BaseModel,                                        \
                dictionary                                                     \
            );                                                                 \
        }                                                                      \
    }


// makeTurbulenceModelTypes and makeRASModel(Type) 
//     from turbulentTransportModels.H
makeTurbulenceModelTypes
(
    geometricOneField,
    geometricOneField,
    incompressibleTurbulenceModel,
    IncompressibleTurbulenceModel,
    transportModel
);


#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (transportModelIncompressibleTurbulenceModel, RAS, Type)


// makeBaseTurbulenceModel from turbulentTransportModels.C

makeBaseTurbulenceModel
(
    geometricOneField,
    geometricOneField,
    incompressibleTurbulenceModel,
    IncompressibleTurbulenceModel,
    transportModel
);

// -------------------------------------------------------------------------- //
// RAS models
// -------------------------------------------------------------------------- //

#include "kOmegaSSTRCH.H"
makeRASModel(kOmegaSSTRCH);

