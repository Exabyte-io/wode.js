import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";
import { mix } from "mixwith";

const defaultHubbardConfig = {
    hubbardUValue: 0.01,
};

export class HubbardContextProviderLegacy extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.uniqueElements = this.material?.Basis?.uniqueElements || [];
    }

    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
                atomicSpeciesIndex: 1,
            },
        ];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
            title: {
                "ui:classNames": "col-xs-12",
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicSpeciesIndex: { ...this.defaultFieldStyles, "ui:readonly": true },
                hubbardUValue: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-07/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U calculation.",
            type: "array",
            uniqueItems: true,
            minItems: 1,
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElements,
                    },
                    atomicSpeciesIndex: {
                        type: "integer",
                        title: "Species index",
                    },
                    hubbardUValue: {
                        type: "number",
                        title: "Hubbard U value",
                        default: defaultHubbardConfig.hubbardUValue,
                    },
                },
                dependencies: {
                    atomicSpecies: {
                        oneOf: this.uniqueElements.map((atom) => {
                            return {
                                properties: {
                                    atomicSpecies: {
                                        type: "string",
                                        enum: [atom],
                                        default: atom,
                                    },
                                    atomicSpeciesIndex: {
                                        type: "integer",
                                        default: this.uniqueElements.indexOf(atom) + 1,
                                    },
                                },
                            };
                        }),
                    },
                },
            },
        };
    }
}
