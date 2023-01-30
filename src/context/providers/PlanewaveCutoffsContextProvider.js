import { Application } from "@exabyte-io/ade.js";
import { ApplicationContextMixin, ContextProvider } from "@exabyte-io/code.js/dist/context";
import { mix } from "mixwith";

const cutoffConfig = {
    vasp: {}, // assuming default cutoffs for VASP
    espresso: {
        // assuming the default GBRV set of pseudopotentials is used
        wavefunction: 40,
        density: 200,
    },
};

export class PlanewaveCutoffsContextProvider extends mix(ContextProvider).with(
    ApplicationContextMixin,
) {
    static Application = Application;

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            wavefunction: {},
            density: {},
        };
    }

    get defaultData() {
        return {
            wavefunction: this.defaultECUTWFC,
            density: this.defaultECUTRHO,
        };
    }

    get _cutoffConfigPerApplication() {
        return cutoffConfig[this.application.name];
    }

    get defaultECUTWFC() {
        return this._cutoffConfigPerApplication.wavefunction || null;
    }

    get defaultECUTRHO() {
        return this._cutoffConfigPerApplication.density || null;
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: " ",
            description:
                "Planewave cutoff parameters for electronic wavefunctions and density. Units are specific to simulation engine.",
            type: "object",
            properties: {
                wavefunction: {
                    type: "number",
                    default: this.defaultECUTWFC,
                },
                density: {
                    type: "number",
                    default: this.defaultECUTRHO,
                },
            },
        };
    }
}
