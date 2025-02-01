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
        let ecutwfc = 0;
        (this.methodData.pseudo || []).forEach((element) => {
            // set the highest of all elements
            if (
                element.cutoffs &&
                element.cutoffs.wfc &&
                element.cutoffs.wfc.standard &&
                element.cutoffs.wfc.standard > ecutwfc
            ) {
                ecutwfc = element.cutoffs.wfc.standard;
            }
        });

        let ecutrho = 0;
        (this.methodData.pseudo || []).forEach((element) => {
            if (
                element.cutoffs &&
                element.cutoffs.rho &&
                element.cutoffs.rho.standard &&
                element.cutoffs.rho.standard > ecutrho
            ) {
                ecutrho = element.cutoffs.rho.standard;
                // if rho cutoff is not present, set it based on wfc cutoff
            } else if (this.methodData.pseudo.type === "us") {
                // if it is ultrasoft pseudopotential set rho cutoff 8 times that of wfc cutoff
                ecutrho = ecutwfc * 8;
            } else {
                // if it is not ultrasoft pseudopotential set rho cutoff 4 times that of wfc cutoff
                ecutrho = ecutwfc * 4;
            }
        });

        return [ecutwfc, ecutrho];
    };

    get defaultECUTWFC() {
        const [ecutwfc] = this.getCutoffsFromPseudos();

        if (ecutwfc > 0) {
            return ecutwfc;
        }

        return this._cutoffConfigPerApplication.wavefunction || null;
    }

    get defaultECUTRHO() {
        const [, ecutrho] = this.getCutoffsFromPseudos();

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
