import { jobContextMixin } from "../../../mixins/JobContextMixin";
import { materialContextMixin } from "../../../mixins/MaterialContextMixin";
import { materialsContextMixin } from "../../../mixins/MaterialsContextMixin";
import { materialsSetContextMixin } from "../../../mixins/MaterialsSetContextMixin";
import { methodDataContextMixin } from "../../../mixins/MethodDataContextMixin";
import { workflowContextMixin } from "../../../mixins/WorkflowContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";
import VASPContextProvider from "./VASPContextProvider";

export default class VASPNEBContextProvider extends ExecutableContextProvider {
    _materials = [];

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
