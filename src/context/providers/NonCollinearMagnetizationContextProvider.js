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
        this.isStartingMagnetization = lodash.get(this.data, "isStartingMagnetization", true);
        this.isConstrainedMagnetization = lodash.get(
            this.data,
            "isConstrainedMagnetization",
            false,
        );
        this.isExistingChargeDensity = lodash.get(this.data, "isExistingChargeDensity", false);
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

        const spinAngles = this.uniqueElementsWithLabels.map((element, index) => {
            return {
                index: index + 1,
                atomicSpecies: element,
                angle1: 0.0,
                angle2: 0.0,
            };
        });

        return {
            isStartingMagnetization: true,
            isConstrainedMagnetization: false,
            isExistingChargeDensity: false,
            startingMagnetization,
            constrainedMagnetization: {
                lambda: 0.0,
                constrainType: "atomic direction",
            },
            spinAngles,
        };
    }

    get uiSchemaStyled() {
        return {
            isStartingMagnetization: {},
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
                "ui:readonly": !this.isStartingMagnetization,
                "ui:options": {
                    addable: false,
                    orderable: false,
                    removable: false,
                },
            },
            isConstrainedMagnetization: {},
            constrainedMagnetization: {
                "ui:classNames": "col-xs-12",
                "ui:readonly": !this.isConstrainedMagnetization,
            },
            isExistingChargeDensity: {},
            spinAngles: {
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
                "ui:readonly": !this.isExistingChargeDensity,
                "ui:options": {
                    addable: false,
                    orderable: false,
                    removable: false,
                },
            },
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
                isStartingMagnetization: {
                    type: "boolean",
                    title: "Set starting magnetization",
                    default: true,
                },
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
                    title: "Set constrained magnetization",
                    default: false,
                },
                constrainedMagnetization: {
                    type: "object",
                    properties: {
                        constrainType: {
                            type: "string",
                            title: "Constrain type",
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
                    },
                },
                isExistingChargeDensity: {
                    type: "boolean",
                    title: "Start calculation from existing charge density",
                    default: false,
                },
                spinAngles: {
                    type: "array",
                    minItems: this.uniqueElementsWithLabels.length,
                    maxItems: this.uniqueElementsWithLabels.length,
                    items: {
                        type: "object",
                        properties: {
                            atomicSpecies: {
                                type: "string",
                                title: "Atomic species",
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
        };
    }
}
