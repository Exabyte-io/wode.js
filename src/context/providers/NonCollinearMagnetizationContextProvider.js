import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@mat3ra/code/dist/js/context";
import { Made } from "@mat3ra/made";
import lodash from "lodash";
import { mix } from "mixwith";

export class NonCollinearMagnetizationContextProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.isConstrainedMagnetization = lodash.get(
            this.data,
            "isConstrainedMagnetization",
            false,
        );
    }

    get uniqueElementsWithLabels() {
        const elementsWithLabelsArray = this.material?.Basis?.elementsWithLabelsArray || [];
        return [...new Set(elementsWithLabelsArray)];
    }

    get defaultData() {
        const startingMagnetization = this.uniqueElementsWithLabels.map((element, index) => {
            return {
                index: index + 1,
                atomicSpecies: element,
                value: 0.0,
            };
        });

        const constrainedMagnetizationValues = this.uniqueElementsWithLabels.map(
            (element, index) => {
                return {
                    index: index + 1,
                    atomicSpecies: element,
                    angle1: 0.0,
                    angle2: 0.0,
                };
            },
        );

        return {
            startingMagnetization,
            isConstrainedMagnetization: false,
            constrainedMagnetization: {
                lambda: 0.0,
                lforcet: false,
                constrainedMagnetizationType: "atomic direction",
                values: constrainedMagnetizationValues,
            },
        };
    }

    // overwrite transformData from parent class
    transformData = (data) => {
        return data;
    };

    get uiSchemaStyled() {
        return {
            startingMagnetization: {
                items: {
                    atomicSpecies: {
                        "ui:classNames": "col-xs-3",
                        "ui:readonly": true,
                    },
                    value: {
                        "ui:classNames": "col-xs-6",
                    },
                },
                "ui:readonly": this.isConstrainedMagnetization,
                "ui:options": {
                    addable: false,
                    orderable: false,
                    removable: false,
                },
            },
            isConstrainedMagnetization: {},
            constrainedMagnetization: {
                values: {
                    items: {
                        atomicSpecies: {
                            "ui:classNames": "col-xs-3",
                            "ui:readonly": true,
                        },
                        angle1: {
                            "ui:classNames": "col-xs-3",
                        },
                        angle2: {
                            "ui:classNames": "col-xs-3",
                        },
                    },
                    "ui:options": {
                        addable: false,
                        orderable: false,
                        removable: false,
                    },
                },
                constrainedMagnetizationType: {
                    "ui:classNames": "col-xs-6",
                },
                lambda: {
                    "ui:classNames": "col-xs-3",
                },
                lforcet: {
                    "ui:classNames": "col-xs-12",
                    "ui:widget": "radio",
                    "ui:inline": true,
                },
            },
            "ui:readonly": !this.isConstrainedMagnetization,
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description:
                "Set initial parameters for non-collinear spin magnetic (SOC) calculation.",
            type: "object",
            properties: {
                startingMagnetization: {
                    type: "array",
                    minItems: this.uniqueElementsWithLabels.length,
                    maxItems: this.uniqueElementsWithLabels.length,
                    items: {
                        type: "object",
                        properties: {
                            atomicSpecies: {
                                type: "string",
                                title: "Atomic species",
                                enum: this.uniqueElementsWithLabels,
                                default: this.firstElement,
                            },
                            value: {
                                type: "number",
                                title: "Starting magnetization",
                                default: 0.0,
                                minimum: -1.0,
                                maximum: 1.0,
                            },
                        },
                    },
                },
                isConstrainedMagnetization: {
                    type: "boolean",
                    title: "Set constrained magnetization instead",
                    default: false,
                },
                constrainedMagnetization: {
                    type: "object",
                    properties: {
                        constrainedMagnetizationType: {
                            type: "string",
                            title: "Constrained magnetization",
                            enum: [
                                "none",
                                "total",
                                "atomic",
                                "total direction",
                                "atomic direction",
                            ],
                            default: "atomic direction",
                        },
                        lambda: {
                            type: "number",
                            title: "lambda",
                            default: 0.0,
                        },
                        lforcet: {
                            type: "boolean",
                            title: "lforcet",
                            oneOf: [
                                { const: true, title: "True" },
                                { const: false, title: "False" },
                            ],
                            default: false,
                        },
                        values: {
                            type: "array",
                            minItems: this.uniqueElementsWithLabels.length,
                            maxItems: this.uniqueElementsWithLabels.length,
                            items: {
                                type: "object",
                                properties: {
                                    atomicSpecies: {
                                        type: "string",
                                        title: "Atomic species",
                                        enum: this.uniqueElementsWithLabels,
                                        default: this.firstElement,
                                    },
                                    angle1: {
                                        type: "number",
                                        title: "Angle1 (deg)",
                                        default: 0.0,
                                    },
                                    angle2: {
                                        type: "number",
                                        title: "Angle2 (deg)",
                                        default: 0.0,
                                    },
                                },
                            },
                        },
                    },
                },
            },
        };
    }
}
