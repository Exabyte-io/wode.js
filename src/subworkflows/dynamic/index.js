import { getQpointIrrep } from "./espresso/getQpointIrrep";
import { getSurfaceEnergySubworkflowUnits } from "./surfaceEnergy";

const dynamicSubworkflowsByApp = {
    espresso: { getQpointIrrep },
};

export { getSurfaceEnergySubworkflowUnits, dynamicSubworkflowsByApp };
