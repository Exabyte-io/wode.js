import { HubbardUContextProvider } from "./HubbardUContextProvider";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "2p",
    atomicSpecies2: "",
    atomicOrbital2: "2p",
    hubbardVValue: 1.0,
};

export class HubbardVContextProvider extends HubbardUContextProvider {
    get defaultData() {
        return [{ ...defaultHubbardConfig, atomicSpecies2: this.uniqueElements }];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: true,
                removable: true,
            },
            title: {
                "ui:classNames": "col-xs-12",
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                atomicSpecies2: this.defaultFieldStyles,
                atomicOrbital2: this.defaultFieldStyles,
                hubbardVValue: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U+V calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElements,
                        default: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Atomic orbital",
                        enum: [
                            "2p",
                            "3s",
                            "3p",
                            "3d",
                            "4s",
                            "4p",
                            "4d",
                            "4f",
                            "5s",
                            "5p",
                            "5d",
                            "5f",
                            "6s",
                            "6p",
                            "6d",
                            "7s",
                            "7p",
                            "7d",
                        ],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    atomicSpecies2: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElements,
                        default: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
                    },
                    atomicOrbital2: {
                        type: "string",
                        title: "Atomic orbital",
                        enum: [
                            "2p",
                            "3s",
                            "3p",
                            "3d",
                            "4s",
                            "4p",
                            "4d",
                            "4f",
                            "5s",
                            "5p",
                            "5d",
                            "5f",
                            "6s",
                            "6p",
                            "6d",
                            "7s",
                            "7p",
                            "7d",
                        ],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    hubbardVValue: {
                        type: "number",
                        title: "Hubbard U (eV)",
                        default: defaultHubbardConfig.hubbardVValue,
                    },
                },
            },
            minItems: 1,
        };
    }
}
