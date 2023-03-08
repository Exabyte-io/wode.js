import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@exabyte-io/code.js/dist/context";
import { math as codeJSMath } from "@exabyte-io/code.js/dist/math";
import { Made } from "@exabyte-io/made.js";
import lodash from "lodash";
import { mix } from "mixwith";
// TODO : pass appSettings to use defaultKPPRA

export class PointsGridFormDataProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this._divisor = config.divisor || 1; // KPPRA will be divided by this number

        this.dimensions = lodash.get(this.data, "dimensions") || this._defaultDimensions;
        this.shifts = lodash.get(this.data, "shifts") || this._defaultShifts;

        // init class fields from data (as constructed from context in parent)
        this.KPPRA = lodash.get(this.data, "KPPRA") || this._defaultKPPRA;
        this.preferKPPRA = lodash.get(this.data, "preferKPPRA", false);
    }

    // eslint-disable-next-line class-methods-use-this
    getDefaultShift() {
        return 0;
    }

    get _defaultDimensions() {
        return this._getGridFromKPPRA(this._defaultKPPRA).dimensions;
    }

    get _defaultShifts() {
        return Array(3).fill(this.getDefaultShift());
    }

    get _defaultKPPRA() {
        return Math.floor(5 / this._divisor);
    }

    get reciprocalVectorRatios() {
        const lattice = new Made.ReciprocalLattice(this.material.lattice);
        return lattice.reciprocalVectorRatios.map((r) =>
            Number(codeJSMath.numberToPrecision(r, 4)),
        );
    }

    get jsonSchema() {
        const kOrQ = this.name[0];
        const vector = {
            type: "array",
            items: {
                type: "number",
            },
            minItems: 3,
            maxItems: 3,
        };

        const vector_ = (defaultValue, isStringType = false) => {
            const isArray = Array.isArray(defaultValue);
            return {
                ...vector,
                items: {
                    type: isStringType ? "string" : "number",
                    ...(isArray ? {} : { default: defaultValue }),
                },
                ...(isArray ? { default: defaultValue } : {}),
            };
        };

        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            description: `3D grid with shifts. Default min value for ${kOrQ.toUpperCase()}PPRA (${kOrQ}pt per reciprocal atom) is ${
                this._defaultKPPRA
            }.`,
            type: "object",
            properties: {
                dimensions: vector_(this._defaultDimensions, this.isUsingJinjaVariables),
                shifts: vector_(this.getDefaultShift()),
                reciprocalVectorRatios: vector_(this.reciprocalVectorRatios),
                KPPRA: {
                    type: "integer",
                    minimum: 1,
                    default: this.KPPRA,
                },
                preferKPPRA: {
                    type: "boolean",
                    default: this.preferKPPRA,
                },
            },
            required: ["dimensions", "shifts"],
        };
    }

    _arraySubStyle(emptyValue = 0) {
        return {
            "ui:options": {
                addable: false,
                orderable: false,
                removable: false,
            },
            items: {
                "ui:disabled": this.preferKPPRA,
                // TODO: extract the actual current values from context
                "ui:placeholder": "1",
                "ui:emptyValue": emptyValue,
            },
        };
    }

    get uiSchema() {
        return {
            dimensions: this._arraySubStyle(1),
            shifts: this._arraySubStyle(0),
            KPPRA: {
                "ui:disabled": !this.preferKPPRA,
                "ui:emptyValue": this.KPPRA,
                "ui:placeholder": this.KPPRA.toString(), // make string to prevent prop type error
            },
            preferKPPRA: {
                ...this.fieldStyles("p-t-20"), // add padding top to level with other elements
                "ui:emptyValue": true,
                "ui:disabled": this.isUsingJinjaVariables,
            },
            reciprocalVectorRatios: {
                "ui:title": "reciprocal vector ratios",
                "ui:orderable": false,
                "ui:removable": false,
                "ui:readonly": true,
            },
        };
    }

    get _defaultData() {
        return {
            dimensions: this._defaultDimensions,
            shifts: this._defaultShifts,
            KPPRA: this._defaultKPPRA,
            preferKPPRA: false,
            reciprocalVectorRatios: this.reciprocalVectorRatios,
        };
    }

    get _defaultDataWithMaterial() {
        // if `data` is present and material is updated, prioritize `data` when `preferKPPRA` is not set
        return this.preferKPPRA
            ? this._getGridFromKPPRA(this.KPPRA)
            : this.data || this._defaultData;
    }

    get defaultData() {
        return this.material ? this._defaultDataWithMaterial : this._defaultData;
    }

    _getGridFromKPPRA(KPPRA) {
        const reciprocalLattice = new Made.ReciprocalLattice(this.material.lattice);
        const nAtoms = this.material ? this.material.Basis.nAtoms : 1;
        return {
            dimensions: reciprocalLattice.getDimensionsFromPoints(KPPRA / nAtoms),
            shifts: this._defaultShifts,
        };
    }

    _getKPPRAFromGrid(grid = this.defaultData) {
        const nAtoms = this.material ? this.material.Basis.nAtoms : 1;
        return grid.dimensions.reduce((a, b) => a * b) * nAtoms;
    }

    static _canTransform(data) {
        return (
            (data.preferKPPRA && data.KPPRA) ||
            (!data.preferKPPRA && data.dimensions.every((d) => typeof d === "number"))
        );
    }

    transformData(data) {
        if (!this.constructor._canTransform(data)) {
            return data;
        }
        // 1. check if KPPRA is preferred
        if (data.preferKPPRA) {
            // 2. KPPRA is preferred => recalculate grid; NOTE: `data.KPPRA` is undefined at first
            data.dimensions = this._getGridFromKPPRA(data.KPPRA).dimensions;
        } else {
            // 3. KPPRA is not preferred => grid was => recalculate KPPRA
            data.KPPRA = this._getKPPRAFromGrid(data);
        }
        return data;
    }
}
