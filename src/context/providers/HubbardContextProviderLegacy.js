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

    /* eslint-disable class-methods-use-this */
    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
            },
        ];
    }

    atomicIndexFromSpecies = (atomicSpecies) => {
        return this.uniqueElements.indexOf(atomicSpecies) + 1;
    };

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
            title: {
                classNames: "col-xs-12",
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
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U calculation.",
            type: "array",
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
                        oneOf: this.uniqueElements.map((e) => {
                            return {
                                properties: {
                                    atomicSpecies: {
                                        enum: [e],
                                    },
                                    atomicSpeciesIndex: {
                                        type: "integer",
                                        enum: [this.uniqueElements.indexOf(e) + 1],
                                        default: this.uniqueElements.indexOf(e) + 1,
                                    },
                                },
                            };
                        }),
                    },
                },
            },
            minItems: 1,
        };
    }
}