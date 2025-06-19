import { BoundaryConditionsFormDataProvider } from "./providers/BoundaryConditionsFormDataProvider";
import QENEBContextProvider from "./providers/by_application/espresso/QENEBContextProvider";
import QEPWXContextProvider from "./providers/by_application/espresso/QEPWXContextProvider";
import NWChemTotalEnergyContextProvider from "./providers/by_application/nwchem/NWChemTotalEnergyContextProvider";
import VASPContextProvider from "./providers/by_application/vasp/VASPContextProvider";
import VASPNEBContextProvider from "./providers/by_application/vasp/VASPNEBContextProvider";
import { CollinearMagnetizationContextProvider } from "./providers/CollinearMagnetizationContextProvider";
import { HubbardContextProviderLegacy } from "./providers/HubbardContextProviderLegacy";
import { HubbardJContextProvider } from "./providers/HubbardJContextProvider";
import { HubbardUContextProvider } from "./providers/HubbardUContextProvider";
import { HubbardVContextProvider } from "./providers/HubbardVContextProvider";
import { IonDynamicsContextProvider } from "./providers/IonDynamicsContextProvider";
import { MLSettingsContextProvider } from "./providers/MLSettingsContextProvider";
import { MLTrainTestSplitContextProvider } from "./providers/MLTrainTestSplitContextProvider";
import { NEBFormDataProvider } from "./providers/NEBFormDataProvider";
import { NonCollinearMagnetizationContextProvider } from "./providers/NonCollinearMagnetizationContextProvider";
import { PlanewaveCutoffsContextProvider } from "./providers/PlanewaveCutoffsContextProvider";
import { PointsGridFormDataProvider } from "./providers/PointsGridFormDataProvider";
import {
    ExplicitPointsPath2PIBAFormDataProvider,
    ExplicitPointsPathFormDataProvider,
    PointsPathFormDataProvider,
} from "./providers/PointsPathFormDataProvider";

export default {
    BoundaryConditionsFormDataProvider,
    MLSettingsContextProvider,
    MLTrainTestSplitContextProvider,
    NEBFormDataProvider,
    PlanewaveCutoffsContextProvider,
    PointsGridFormDataProvider,
    PointsPathFormDataProvider,
    ExplicitPointsPathFormDataProvider,
    ExplicitPointsPath2PIBAFormDataProvider,
    HubbardJContextProvider,
    HubbardUContextProvider,
    HubbardVContextProvider,
    HubbardContextProviderLegacy,
    IonDynamicsContextProvider,
    CollinearMagnetizationContextProvider,
    NonCollinearMagnetizationContextProvider,
    VASPContextProvider,
    VASPNEBContextProvider,
    QEPWXContextProvider,
    QENEBContextProvider,
    NWChemTotalEnergyContextProvider,
};
