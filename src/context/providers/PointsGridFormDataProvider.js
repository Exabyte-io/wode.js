import JSONSchemaFormDataProvider from "@exabyte-io/ade.js/dist/js/context/JSONSchemaFormDataProvider";
import { units as UNITS } from "@mat3ra/code/dist/js/constants";
import { math as codeJSMath } from "@mat3ra/code/dist/js/math";
import { Made } from "@mat3ra/made";
import lodash from "lodash";

import { materialContextMixin } from "../mixins/MaterialContextMixin";
import { globalSettings } from "./settings";

export class PointsGridFormDataProvider extends JSONSchemaFormDataProvider {
    constructor(config) {
        super(config);
        this.initMaterialContextMixin();

        this._divisor = config.divisor || 1; // KPPRA will be divided by this number
        this.reciprocalLattice = new Made.ReciprocalLattice(this.material.lattice);

        this.dimensions = lodash.get(this.data, "dimensions") || this._defaultDimensions;
        this.shifts = lodash.get(this.data, "shifts") || this._defaultShifts;

        // init class fields from data (as constructed from context in parent)
        this.gridMetricType = lodash.get(this.data, "gridMetricType") || "KPPRA";
        this.gridMetricValue =
            lodash.get(this.data, "gridMetricValue") || this._getDefaultGridMetricValue("KPPRA");
        this.preferGridMetric = lodash.get(this.data, "preferGridMetric", false);

        this._metricDescription = {
            KPPRA: `${this.name[0].toUpperCase()}PPRA (${this.name[0]}pt per reciprocal atom)`, // KPPRA or QPPRA
            spacing: "grid spacing",
        };
        this.defaultClassNames = "col-xs-12 col-sm-6 col-md-3 col-lg-2";
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
                return Math.floor(globalSettings.defaultKPPRA / this._divisor);
            case "spacing":
                return 0.3;
            default:
                console.error("Metric type not recognized!");
                return 1;
        }
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

    get reciprocalVectorRatios() {
        return this.reciprocalLattice.reciprocalVectorRatios.map((r) =>
            Number(codeJSMath.numberToPrecision(r, 3)),
        );
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
            $schema: "http://json-schema.org/draft-07/schema#",
            description: `3D grid with shifts. Default min value for ${
                this._metricDescription[this.gridMetricType]
            } is ${this._getDefaultGridMetricValue(this.gridMetricType)}.`,
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

    get uiSchema() {
        const _arraySubStyle = (emptyValue = 0) => {
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
                    "ui:label": false,
                },
            };
        };

        return {
            dimensions: _arraySubStyle(1),
            shifts: _arraySubStyle(0),
            gridMetricType: {
                "ui:title": "Grid Metric",
            },
            gridMetricValue: {
                "ui:disabled": !this.preferGridMetric,
                "ui:emptyValue": this.gridMetricValue,
                "ui:placeholder": this.gridMetricValue.toString(), // make string to prevent prop type error
            },
            preferGridMetric: {
                "ui:emptyValue": true,
                "ui:disabled": this.isUsingJinjaVariables,
            },
            reciprocalVectorRatios: {
                "ui:title": "reciprocal vector ratios",
                "ui:orderable": false,
                "ui:removable": false,
                "ui:readonly": true,
                items: {
                    "ui:label": false,
                },
            },
        };
    }

    _getDimensionsFromKPPRA(KPPRA) {
        const nAtoms = this.material ? this.material.Basis.nAtoms : 1;
        return this.reciprocalLattice.getDimensionsFromPointsCount(KPPRA / nAtoms);
    }

    _getKPPRAFromDimensions(dimensions) {
        const nAtoms = this.material ? this.material.Basis.nAtoms : 1;
        return dimensions.reduce((a, b) => a * b) * nAtoms;
    }

    static _canTransform(data) {
        return (
            (data.preferGridMetric && data.gridMetricType && data.gridMetricValue) ||
            (!data.preferGridMetric && data.dimensions.every((d) => typeof d === "number"))
        );
    }

    calculateDimensions({ gridMetricType, gridMetricValue, units = UNITS.angstrom }) {
        switch (gridMetricType) {
            case "KPPRA":
                return this._getDimensionsFromKPPRA(gridMetricValue);
            case "spacing":
                return this.reciprocalLattice.getDimensionsFromSpacing(gridMetricValue, units);
            default:
                return [1, 1, 1];
        }
    }

    calculateGridMetric({ gridMetricType, dimensions, units = UNITS.angstrom }) {
        switch (gridMetricType) {
            case "KPPRA":
                return this._getKPPRAFromDimensions(dimensions);
            case "spacing":
                return lodash.round(
                    this.reciprocalLattice.getSpacingFromDimensions(dimensions, units),
                    3,
                );
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

materialContextMixin(PointsGridFormDataProvider.prototype);
