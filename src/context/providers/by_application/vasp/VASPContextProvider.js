import { jobContextMixin } from "../../../mixins/JobContextMixin";
import { materialContextMixin } from "../../../mixins/MaterialContextMixin";
import { materialsContextMixin } from "../../../mixins/MaterialsContextMixin";
import { methodDataContextMixin } from "../../../mixins/MethodDataContextMixin";
import { workflowContextMixin } from "../../../mixins/WorkflowContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";

export default class VASPContextProvider extends ExecutableContextProvider {
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
