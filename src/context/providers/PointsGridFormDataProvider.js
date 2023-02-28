import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@exabyte-io/code.js/dist/context";
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
        this.gridMetricType = lodash.get(this.data, "gridMetricType") || "KPPRA";
        this.gridMetricValue =
            lodash.get(this.data, "gridMetricValue") || this._getDefaultGridMetricValue("KPPRA");
        this.preferGridMetric = lodash.get(this.data, "preferGridMetric", false);
    }

    // eslint-disable-next-line class-methods-use-this
    getDefaultShift() {
        return 0;
    }

    get _defaultDimensions() {
        return this.calculateDimensions({
            gridMetricType: "KPPRA",
            gridMetricValue: this._getDefaultGridMetricValue("KPPRA"),
        });
    }

    get _defaultShifts() {
        return Array(3).fill(this.getDefaultShift());
    }

    _getDefaultGridMetricValue(metric) {
        switch (metric) {
            case "KPPRA":
                return Math.floor(5 / this._divisor);
            case "spacing":
                return 0.3;
            default:
                console.error("Metric type not recognized!");
                return 1;
        }
    }

    get _descriptionText() {
        let metric;
        if (this.gridMetricType === "KPPRA") {
            const kOrQ = this.name[0];
            metric = `${kOrQ.toUpperCase()}PPRA (${kOrQ}pt per reciprocal atom)`;
        } else if (this.gridMetricType === "spacing") {
            metric = "grid spacing";
        }
        return `Default min value for ${metric} is ${this._getDefaultGridMetricValue(
            this.gridMetricType,
        )}.`;
    }

    get reciprocalVectorRatios() {
        if (!this.material) return [1, 1, 1];
        const lattice = new Made.ReciprocalLattice(this.material.lattice);
        return lattice.reciprocalVectorRatios.map((r) => lodash.round(r, 3));
    }

    get jsonSchema() {
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
            description: `3D grid with shifts. ${this._descriptionText}`,
            type: "object",
            properties: {
                dimensions: vector_(this._defaultDimensions, this.isUsingJinjaVariables),
                shifts: vector_(this.getDefaultShift()),
                reciprocalVectorRatios: vector_(this.reciprocalVectorRatios),
                gridMetricType: {
                    type: "string",
                    enum: ["KPPRA", "spacing"],
                    default: "KPPRA",
                },
                gridMetricValue: {
                    type: "number",
                },
                preferGridMetric: {
                    type: "boolean",
                },
            },
            dependencies: {
                gridMetricType: {
                    oneOf: [
                        {
                            properties: {
                                gridMetricType: {
                                    enum: ["KPPRA"],
                                },
                                gridMetricValue: {
                                    type: "integer",
                                    minimum: 1,
                                    title: "Value",
                                    default: this.gridMetricValue,
                                },
                                preferGridMetric: {
                                    type: "boolean",
                                    title: "prefer KPPRA",
                                    default: this.preferGridMetric,
                                },
                            },
                        },
                        {
                            properties: {
                                gridMetricType: {
                                    enum: ["spacing"],
                                },
                                gridMetricValue: {
                                    type: "number",
                                    minimum: 0,
                                    title: "Value [1/Å]",
                                    default: this.gridMetricValue,
                                },
                                preferGridMetric: {
                                    type: "boolean",
                                    title: "prefer spacing",
                                    default: this.preferGridMetric,
                                },
                            },
                        },
                    ],
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
                "ui:disabled": this.preferGridMetric,
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
            gridMetricType: {
                "ui:title": "Grid Metric",
            },
            gridMetricValue: {
                "ui:disabled": !this.preferGridMetric,
                "ui:emptyValue": this.gridMetricValue,
                "ui:placeholder": this.gridMetricValue.toString(), // make string to prevent prop type error
            },
            preferGridMetric: {
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
            gridMetricType: "KPPRA",
            gridMetricValue: this._getDefaultGridMetricValue("KPPRA"),
            preferGridMetric: false,
            reciprocalVectorRatios: this.reciprocalVectorRatios,
        };
    }

    get _defaultDataWithMaterial() {
        const { gridMetricType, gridMetricValue } = this;
        // if `data` is present and material is updated, prioritize `data` when `preferGridMetric` is not set
        return this.preferGridMetric
            ? {
                  dimensions: this.calculateDimensions({ gridMetricType, gridMetricValue }),
                  shifts: this._defaultShifts,
              }
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

    _getDimensionsFromSpacing(spacing) {
        const reciprocalLattice = new Made.ReciprocalLattice(this.material.lattice);
        return reciprocalLattice.getDimensionsFromSpacing(spacing);
    }

    _getSpacingFromDimensions(dimensions) {
        const reciprocalLattice = new Made.ReciprocalLattice(this.material.lattice);
        return reciprocalLattice.getSpacingFromDimensions(dimensions);
    }

    static _canTransform(data) {
        return (
            (data.preferGridMetric && data.gridMetricType && data.gridMetricValue) ||
            (!data.preferGridMetric && data.dimensions.every((d) => typeof d === "number"))
        );
    }

    calculateDimensions({ gridMetricType, gridMetricValue }) {
        switch (gridMetricType) {
            case "KPPRA":
                return this._getGridFromKPPRA(gridMetricValue).dimensions;
            case "spacing":
                return this._getDimensionsFromSpacing(gridMetricValue);
            default:
                return [1, 1, 1];
        }
    }

    calculateGridMetric({ gridMetricType, dimensions }) {
        switch (gridMetricType) {
            case "KPPRA":
                return this._getKPPRAFromGrid({ dimensions });
            case "spacing":
                return lodash.round(this._getSpacingFromDimensions(dimensions), 3);
            default:
                return 1;
        }
    }

    transformData(data) {
        if (!this.constructor._canTransform(data)) {
            return data;
        }
        // dimensions are calculated from grid metric or vice versa
        if (data.preferGridMetric) {
            data.dimensions = this.calculateDimensions(data);
        } else {
            data.gridMetricValue = this.calculateGridMetric(data);
        }
        return data;
    }
}
