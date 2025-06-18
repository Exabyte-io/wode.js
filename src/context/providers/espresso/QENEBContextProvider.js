import { jobContextMixin } from "@mat3ra/code/dist/js/context/JobContextMixin";
import {
    // type MethodDataContextMixinType,
    methodDataContextMixin,
} from "@mat3ra/code/dist/js/context/MethodDataContextMixin";
// import type { ContextProviderConfig } from "@mat3ra/code/dist/js/context/provider";
// import type { InMemoryEntity } from "@mat3ra/code/dist/js/entity";
// import type { OrderedInMemoryEntityInSet } from "@mat3ra/code/dist/js/entity/set/ordered/OrderedInMemoryEntityInSetMixin";
// import type { Constructor } from "@mat3ra/code/dist/js/utils/types";
import { Made } from "@mat3ra/made";
// import type { MaterialMixin } from "@mat3ra/made/dist/js/materialMixin";
import lodash from "lodash";

// import type { JobContextMixinType } from "../../mixins/JobContextMixin";
// import type { MaterialContextMixinType } from "../../mixins/MaterialContextMixin";
import { materialContextMixin } from "../../mixins/MaterialContextMixin";
// import type { MaterialsContextMixinType } from "../../mixins/MaterialsContextMixin";
import { materialsContextMixin } from "../../mixins/MaterialsContextMixin";
import {
    // type MaterialsSetContextMixinType,
    materialsSetContextMixin,
} from "../../mixins/MaterialsSetContextMixin";
import {
    // type WorkflowContextMixinType,
    workflowContextMixin,
} from "../../mixins/WorkflowContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";
import QEPWXContextProvider from "./QEPWXContextProvider";

// export type Base = typeof ExecutableContextProvider &
//     Constructor<MaterialContextMixinType> &
//     Constructor<MaterialsContextMixinType> &
//     Constructor<MethodDataContextMixinType> &
//     Constructor<WorkflowContextMixinType> &
//     Constructor<JobContextMixinType> &
//     Constructor<MaterialsSetContextMixinType>;

// export type Material = MaterialMixin & InMemoryEntity & OrderedInMemoryEntityInSet;

export default class QENEBContextProvider extends ExecutableContextProvider {
    static Material = Made.Material;

    _material = undefined;

    _materials = [];

    _materialsSet = undefined;

    constructor(config) {
        super(config);
        this.initJobContextMixin();
        this.initMaterialsContextMixin();
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initMaterialContextMixin();
        this.initMaterialsSetContextMixin();
    }

    getData() {
        const sortedMaterials = this.sortMaterialsByIndexInSet(this.materials);
        const PWXContexts = sortedMaterials.map((material) => {
            const context = { ...this.config.context, material };
            const config = { ...this.config, context };
            return new QEPWXContextProvider(config).getData();
        });

        return {
            ...lodash.omit(PWXContexts[0], [
                "ATOMIC_POSITIONS",
                "ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS",
            ]),
            FIRST_IMAGE: PWXContexts[0].ATOMIC_POSITIONS,
            LAST_IMAGE: PWXContexts[PWXContexts.length - 1].ATOMIC_POSITIONS,
            INTERMEDIATE_IMAGES: PWXContexts.slice(1, PWXContexts.length - 1).map(
                (data) => data.ATOMIC_POSITIONS,
            ),
        };
    }
}

materialContextMixin(QENEBContextProvider.prototype);
materialsContextMixin(QENEBContextProvider.prototype);
methodDataContextMixin(QENEBContextProvider.prototype);
workflowContextMixin(QENEBContextProvider.prototype);
jobContextMixin(QENEBContextProvider.prototype);
materialsSetContextMixin(QENEBContextProvider.prototype);
