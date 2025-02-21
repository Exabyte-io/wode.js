import { Application } from "@exabyte-io/ade.js";
import {
    ApplicationContextMixin,
    ContextProvider,
    MethodDataContextMixin,
} from "@mat3ra/code/dist/js/context";
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
    MethodDataContextMixin,
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

    get _cutoffsFromPseudos() {
        const pseudos = this.methodData?.pseudo || [];
        let ecutwfc = 0;
        let ecutrho = 0;

        pseudos.forEach((data) => {
            // set the highest cutoff of all elements
            if (data?.cutoffs?.wavefunction?.standard > ecutwfc) {
                ecutwfc = data.cutoffs.wavefunction.standard;
            }

            if (data?.cutoffs?.density?.standard > ecutrho) {
                ecutrho = data.cutoffs.density.standard;
            } else {
                // if rho cutoff is not present, set it based on wfc cutoff
                // if it is ultrasoft pseudopotential set rho cutoff 8 times
                // that of wfc cutoff, otherwise 4 times that of wfc cutoff
                const rhoMultiplier = this.methodData?.pseudo?.type === "us" ? 8 : 4;
                ecutrho = Math.max(ecutrho, ecutwfc * rhoMultiplier);
            }
        });

        return [ecutwfc, ecutrho];
    }

    get defaultECUTWFC() {
        if (["espresso", "qe"].includes(this.application.shortName)) {
            const [ecutwfc] = this._cutoffsFromPseudos;

            if (ecutwfc > 0) {
                return ecutwfc;
            }
        }

        return this._cutoffConfigPerApplication.wavefunction || null;
    }

    get defaultECUTRHO() {
        if (["espresso", "qe"].includes(this.application.shortName)) {
            const [, ecutrho] = this._cutoffsFromPseudos;

            if (ecutrho > 0) {
                return ecutrho;
            }
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
