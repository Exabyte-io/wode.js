import { PERIODIC_TABLE } from "@exabyte-io/periodic-table.js";
import lodash from "lodash";
import _ from "underscore";
import s from "underscore.string";

import { jobContextMixin } from "../../../mixins/JobContextMixin";
import { materialContextMixin } from "../../../mixins/MaterialContextMixin";
import { methodDataContextMixin } from "../../../mixins/MethodDataContextMixin";
import { workflowContextMixin } from "../../../mixins/WorkflowContextMixin";
import ExecutableContextProvider from "../ExecutableContextProvider";

export default class NWChemTotalEnergyContextProvider extends ExecutableContextProvider {
    _material = undefined;

    constructor(config) {
        super(config);
        this.initMethodDataContextMixin();
        this.initWorkflowContextMixin();
        this.initJobContextMixin();
        this.initMaterialContextMixin();
    }

    get atomicPositionsWithoutConstraints() {
        return this.material.Basis.atomicPositions;
    }

    get atomicPositions() {
        const basis = this.material.Basis;
        basis.toCartesian();
        return basis.getAtomicPositionsWithConstraintsAsStrings();
    }

    get atomSymbols() {
        return this.material.Basis.uniqueElements;
    }

    get cartesianAtomicPositions() {
        return this.material.Basis.toCartesian !== undefined;
    }

    get ATOMIC_SPECIES() {
        return _.map(this.atomSymbols, (symbol) => {
            return NWChemTotalEnergyContextProvider.symbolToAtomicSpecies(symbol);
        }).join("\n");
    }

    /*
     * @NOTE: Overriding getData makes this provider "stateless", ie. delivering data from scratch each time and not
     *        considering the content of `this.data`, and `this.isEdited` field(s).
     */
    getData() {
        /*
        TODO: Create ability for user to define CHARGE, MULT, BASIS and FUNCTIONAL parameters.
         */
        const CHARGE = 0;
        const MULT = 1;
        const BASIS = "6-31G";
        const FUNCTIONAL = "B3LYP";

        return {
            CHARGE,
            MULT,
            BASIS,
            NAT: this.atomicPositions.length,
            NTYP: this.atomSymbols.length,
            ATOMIC_POSITIONS: this.atomicPositions.join("\n"),
            ATOMIC_POSITIONS_WITHOUT_CONSTRAINTS: this.atomicPositionsWithoutConstraints.join("\n"),
            ATOMIC_SPECIES: this.ATOMIC_SPECIES,
            FUNCTIONAL,
            CARTESIAN: this.cartesianAtomicPositions,
        };
    }

    static symbolToAtomicSpecies(symbol, pseudo) {
        const el = PERIODIC_TABLE[symbol];
        const filename = pseudo
            ? lodash.get(pseudo, "filename", s.strRightBack(pseudo.path || "", "/"))
            : "";
        return el ? s.sprintf("%s %f %s", symbol, el.atomic_mass, filename) : undefined;
    }
}

materialContextMixin(NWChemTotalEnergyContextProvider.prototype);
methodDataContextMixin(NWChemTotalEnergyContextProvider.prototype);
workflowContextMixin(NWChemTotalEnergyContextProvider.prototype);
jobContextMixin(NWChemTotalEnergyContextProvider.prototype);
