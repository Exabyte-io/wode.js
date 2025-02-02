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

    getCutoffsFromPseudos = () => {
        const pseudos = this.methodData.pseudo || [];
        let ecutwfc = 0;
        let ecutrho = 0;

        pseudos.forEach((data) => {
            // set the highest cutoff of all elements
            if (data?.cutoffs?.wfc?.standard > ecutwfc) {
                ecutwfc = data.cutoffs.wfc.standard;
            }

            if (data?.cutoffs?.rho?.standard > ecutrho) {
                ecutrho = data.cutoffs.rho.standard;
            } else if (this.methodData.pseudo?.type === "us" && ecutwfc * 8 > ecutrho) {
                // if rho cutoff is not present, set it based on wfc cutoff
                // if it is ultrasoft pseudopotential set rho cutoff 8 times that of wfc cutoff
                ecutrho = ecutwfc * 8;
            } else if (ecutwfc * 4 > ecutrho) {
                // if it is not ultrasoft pseudopotential set rho cutoff 4 times that of wfc cutoff
                ecutrho = ecutwfc * 4;
            }
        });

        return [ecutwfc, ecutrho];
    };

    get defaultECUTWFC() {
        if (["espresso", "qe"].includes(this.application.shortName)) {
            const [ecutwfc] = this.getCutoffsFromPseudos();

            if (ecutwfc > 0) {
                return ecutwfc;
            }
        }

        return this._cutoffConfigPerApplication.wavefunction || null;
    }

    get defaultECUTRHO() {
        if (["espresso", "qe"].includes(this.application.shortName)) {
            const [, ecutrho] = this.getCutoffsFromPseudos();

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
