// eslint-disable-next-line no-unused-vars
import { createSubworkflowByName } from "../subworkflows";

export const RelaxationLogicMixin = (superclass) =>
    class extends superclass {
        // TODO: figure out how to avoid circular dependency on import in the platform webapp and re-enable or remove
        // get _allRelaxationSubworkflows() {
        //     /*
        //     * NOTE: dynamic import is used here to avoid import errors on application start.
        //     *       When `Workflow` is used in `AccountSelector` it attempts to import relaxation subworkflows.
        //     *       The latter is constructed through using DAOProvider. Because DAOProvider is not yet defined at the point
        //     *       when `AccountSelector` is imported the whole logic fails and prevent application start.
        //      */
        //     console.log("_allRelaxationSubworkflows", this.constructor._allRelaxationSubworkflows);
        //     return this.constructor._allRelaxationSubworkflows ? this.constructor._allRelaxationSubworkflows : {
        //         espresso: createSubworkflowByName({ appName: "espresso", swfName: "variable_cell_relaxation" }),
        //         vasp: createSubworkflowByName({ appName: "vasp", swfName: "variable_cell_relaxation" }),
        //     };
        // }

        get relaxationSubworkflow() {
            // deciding on the application based on the first subworkflow
            const firstSubworkflow = this.subworkflows[0];
            return this._allRelaxationSubworkflows[firstSubworkflow.application.name];
        }

        isRelaxationSubworkflow(subworkflow) {
            return Object.values(this._allRelaxationSubworkflows)
                .map((sw) => sw.systemName)
                .includes(subworkflow.systemName);
        }

        get hasRelaxation() {
            return Boolean(
                this.subworkflows.find((subworkflow) => {
                    return this.isRelaxationSubworkflow(subworkflow);
                }),
            );
        }

        toggleRelaxation() {
            if (this.hasRelaxation) {
                const relaxSubworkflow = this.subworkflows.find((sw) =>
                    this.isRelaxationSubworkflow(sw),
                );
                this.removeSubworkflow(relaxSubworkflow.id);
            } else {
                const vcRelax = this.relaxationSubworkflow;
                this.addSubworkflow(vcRelax, true);
            }
        }
    };
