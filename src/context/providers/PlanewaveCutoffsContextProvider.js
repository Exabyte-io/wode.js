import ContextProvider from "@exabyte-io/ade.js/dist/js/context/ContextProvider";

import { applicationContextMixin } from "../mixins/ApplicationContextMixin";
import { methodDataContextMixin } from "../mixins/MethodDataContextMixin";

const cutoffConfig = {
    vasp: {}, // assuming default cutoffs for VASP
    espresso: {
        // assuming the default GBRV set of pseudopotentials is used
        wavefunction: 40,
        density: 200,
    },
};

export class PlanewaveCutoffsContextProvider extends ContextProvider {
    constructor(config) {
        super(config);
        this.initMethodDataContextMixin();
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
        const [ecutwfc] = this.highestCutoffsFromPseudos;

        if (ecutwfc > 0) {
            return ecutwfc;
        }

        return this._cutoffConfigPerApplication.wavefunction || null;
    }

    get defaultECUTRHO() {
        const [, ecutrho] = this.highestCutoffsFromPseudos; // destructure and select second item

        if (ecutrho > 0) {
            return ecutrho;
        }

        return this._cutoffConfigPerApplication.density || null;
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
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

methodDataContextMixin(PlanewaveCutoffsContextProvider.prototype);
applicationContextMixin(PlanewaveCutoffsContextProvider.prototype);
