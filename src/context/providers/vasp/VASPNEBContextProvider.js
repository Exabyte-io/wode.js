import { jobContextMixin } from "@mat3ra/code/dist/js/context/JobContextMixin";
import { materialsContextMixin } from "@mat3ra/code/dist/js/context/MaterialsContextMixin";
// import type { MethodDataContextMixinType } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import { methodDataContextMixin } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
// import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import { workflowContextMixin } from "@mat3ra/code/dist/js/context/WorkflowContextMixin";
import { Made } from "@mat3ra/made";

// import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
// import type { JobContextMixinType } from "../../mixins/JobContextMixin";
import {
    // type MaterialContextMixinType,
    materialContextMixin,
} from "../../mixins/MaterialContextMixin";
// import type { MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
import {
    // type MaterialsSetContextMixinType,
    materialsSetContextMixin,
} from "../../mixins/MaterialsSetContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";
// import type { WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
// import { ExecutableContextProvider } from "../../providers";
// import type { Material } from "../espresso/QENEBContextProvider";
import VASPContextProvider from "./VASPContextProvider";

// type Base = typeof ExecutableContextProvider &
//     Constructor<MaterialContextMixinType> &
//     Constructor<MaterialsContextMixinType> &
//     Constructor<MaterialsSetContextMixinType> &
//     Constructor<MethodDataContextMixinType> &
//     Constructor<WorkflowContextMixinType> &
//     Constructor<JobContextMixinType>;

export default class VASPNEBContextProvider extends ExecutableContextProvider {
    _materials = [];

    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.initMaterialContextMixin();
        this.initMaterialsContextMixin();
        this.initMaterialsSetContextMixin();
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initJobContextMixin();
    }

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const VASPContexts = sortedMaterials.map((material) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new VASPContextProvider(config).getData();
        });

        return {
            FIRST_IMAGE: VASPContexts[0].POSCAR_WITH_CONSTRAINTS,
            LAST_IMAGE: VASPContexts[VASPContexts.length - 1].POSCAR_WITH_CONSTRAINTS,
            INTERMEDIATE_IMAGES: VASPContexts.slice(1, VASPContexts.length - 1).map(
                (data) => data.POSCAR_WITH_CONSTRAINTS,
            ),
        };
    }
}

materialContextMixin(VASPNEBContextProvider.prototype);
materialsContextMixin(VASPNEBContextProvider.prototype);
materialsSetContextMixin(VASPNEBContextProvider.prototype);
methodDataContextMixin(VASPNEBContextProvider.prototype);
workflowContextMixin(VASPNEBContextProvider.prototype);
jobContextMixin(VASPNEBContextProvider.prototype);
