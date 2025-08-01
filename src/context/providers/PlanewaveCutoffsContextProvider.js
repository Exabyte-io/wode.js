import ContextProvider from "@exabyte-io/ade.js/dist/js/context/ContextProvider";

import { applicationContextMixin } from "../mixins/ApplicationContextMixin";

const cutoffConfig = {
    vasp: {
        // assuming default cutoffs for VASP as specified in POTCAR
        cutoffUnit: "eV",
    },
    espresso: {
        // assuming the default GBRV set of pseudopotentials is used
        wavefunction: 40,
        density: 200,
        cutoffUnit: "Ry",
    },
};

export class PlanewaveCutoffsContextProvider extends ContextProvider {
    constructor(config) {
        super(config);
        this.initApplicationContextMixin();
    }

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
        const descriptionText = `Planewave cutoff parameters for electronic \
            wavefunctions and density. Cutoffs are expressed in \
            ${this._cutoffConfigPerApplication.cutoffUnit} for \
            ${this.application.name}.`;
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: " ",
            description: descriptionText,
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

applicationContextMixin(PlanewaveCutoffsContextProvider.prototype);
