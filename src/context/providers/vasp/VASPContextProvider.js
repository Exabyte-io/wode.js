// import type { MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
// import type { WorkflowContextMixinType } from "../../mixins/WorkflowContextMixin";
// import type { Material } from "../espresso/QENEBContextProvider";
import { jobContextMixin } from "@mat3ra/code/dist/js/context/JobContextMixin";
import { materialsContextMixin } from "@mat3ra/code/dist/js/context/MaterialsContextMixin";
// import type { MethodDataContextMixinType } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
import { methodDataContextMixin } from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
// import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
import { workflowContextMixin } from "@mat3ra/code/dist/js/context/WorkflowContextMixin";
// import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { Made } from "@mat3ra/made";

// import type { JobContextMixinType } from "../../mixins/JobContextMixin";
import {
    // type MaterialContextMixinType,
    materialContextMixin,
} from "../../mixins/MaterialContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";

// export type Base = typeof ExecutableContextProvider &
//     Constructor<MaterialContextMixinType> &
//     Constructor<MaterialsContextMixinType> &
//     Constructor<MethodDataContextMixinType> &
//     Constructor<WorkflowContextMixinType> &
//     Constructor<JobContextMixinType>;

export default class VASPContextProvider extends ExecutableContextProvider {
    static Material = Made.Material;

    _material = undefined;

    _materials = [];

    constructor(config) {
        super(config);
        this.initJobContextMixin();
        this.initMaterialsContextMixin();
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initMaterialContextMixin();
    }

    // eslint-disable-next-line class-methods-use-this
    buildVASPContext(material) {
        return {
            // TODO: figure out whether we need two separate POSCARS, maybe one is enough
            POSCAR: material.getAsPOSCAR(true, true),
            POSCAR_WITH_CONSTRAINTS: material.getAsPOSCAR(true),
        };
    }

    getDataPerMaterial() {
        if (!this.materials || this.materials.length <= 1) return {};
        return { perMaterial: this.materials.map((material) => this.buildVASPContext(material)) };
    }

    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        // consider adjusting so that below values are read from PlanewaveDataManager
        // ECUTWFC;
        // ECUTRHO;

        return {
            ...this.buildVASPContext(this.material),
            ...this.getDataPerMaterial(),
        };
    }
}

materialContextMixin(VASPContextProvider.prototype);
materialsContextMixin(VASPContextProvider.prototype);
methodDataContextMixin(VASPContextProvider.prototype);
workflowContextMixin(VASPContextProvider.prototype);
jobContextMixin(VASPContextProvider.prototype);
