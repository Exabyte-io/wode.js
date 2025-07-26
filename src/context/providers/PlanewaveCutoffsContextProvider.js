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
        this.initApplicationContextMixin();
        this.initMethodDataContextMixin();
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

    /**
     * @summary Find cutoff values corresponding to wavefunction or charge
     * density with accuracy level in the pseudo method data
     * @param data {Object}: the data object from pseudo method data
     * @param cutoffEntity {String}: "wavefunction", or "density"
     * @param accuracyLevel {String}: "standard", "high", or "low"
     * @return {number}: if cutoff value present returns value else return 0
     */
    _cutoffFromPseudoData = (data, cutoffEntity, accuracyLevel = "standard") => {
        const cutoff = data?.cutoffs?.[cutoffEntity] || [];
        return cutoff.find((obj) => obj.accuracy_level === accuracyLevel)?.value ?? 0;
    };

    /**
     * @summary Set highest cutoff of all elements present, in case wavefunction
     * cutoff is present but no charge density cutoff, set charge density cutoff
     * to 4 or 8 times that of wavefunction cutoff
     * @return {Array<number>}: tuple of wavefunction and density cutoffs
     */
    get _cutoffsFromPseudos() {
        let ecutwfc = 0;
        let ecutrho = 0;
        const pseudos = this.methodData?.pseudo || [];

        pseudos.forEach((data) => {
            const ecutwfcStandard = this._cutoffFromPseudoData(data, "wavefunction");
            const ecutrhoStandard = this._cutoffFromPseudoData(data, "density");
            // set the highest cutoff of all elements
            ecutwfc = Math.max(ecutwfc, ecutwfcStandard);

            if (ecutrhoStandard > ecutrho) {
                ecutrho = ecutrhoStandard;
            } else {
                // if rho cutoff is not present, set it based on wfc cutoff
                // if it is ultrasoft pseudopotential set rho cutoff 8 times
                // that of wfc cutoff, otherwise 4 times that of wfc cutoff
                const rhoMultiplier = data?.type === "us" ? 8 : 4;
                ecutrho = Math.max(ecutrho, ecutwfc * rhoMultiplier);
            }
        });

        return [ecutwfc, ecutrho];
    }

    get defaultECUTWFC() {
        const [ecutwfc] = this._cutoffsFromPseudos;

        if (ecutwfc > 0) {
            return ecutwfc;
        }

        return this._cutoffConfigPerApplication.wavefunction || null;
    }

    get defaultECUTRHO() {
        const [, ecutrho] = this._cutoffsFromPseudos; // destructure and select second item

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

applicationContextMixin(PlanewaveCutoffsContextProvider.prototype);
methodDataContextMixin(PlanewaveCutoffsContextProvider.prototype);
