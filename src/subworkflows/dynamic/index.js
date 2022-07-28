import { getSurfaceEnergySubworkflowUnits } from "./surfaceEnergy";
import { getQpointIrrep } from "./espresso/getQpointIrrep";

const dynamicSubworkflowsByApp = {
    espresso: { getQpointIrrep },
}

export {
    getSurfaceEnergySubworkflowUnits,
    dynamicSubworkflowsByApp,
};
