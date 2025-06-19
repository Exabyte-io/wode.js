import lodash from "lodash";

import { jobContextMixin } from "../../../mixins/JobContextMixin";
import { materialContextMixin } from "../../../mixins/MaterialContextMixin";
import { materialsContextMixin } from "../../../mixins/MaterialsContextMixin";
import { materialsSetContextMixin } from "../../../mixins/MaterialsSetContextMixin";
import { methodDataContextMixin } from "../../../mixins/MethodDataContextMixin";
import { workflowContextMixin } from "../../../mixins/WorkflowContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";
import QEPWXContextProvider from "./QEPWXContextProvider";

export default class QENEBContextProvider extends ExecutableContextProvider {
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
